/******************************************************************************************************/
/*                           HYREC: Hydrogen and Helium Recombination Code                            */
/*                      Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                      */
/*                                                                                                    */
/*         history.c: functions for recombination history                                             */
/*                                                                                                    */
/*         Version: January 2015                                                                      */
/*                                                                                                    */
/*         Revision history:                                                                          */
/*            - written November 2010                                                                 */
/*            - January 2011: changed various switches (notably for post-Saha expansions)             */
/*                             so that they remain valid for arbitrary cosmologies                    */
/*            - May 2012:   - added explicit dependence on fine structure constant and electron mass  */
/*                          - modified call of rec_build_history                                      */
/*                             and improved numerical radiative transfer equations                    */
/*                             so the Lyman-lines spectrum can be extracted                           */
/*                           - split some functions for more clarity                                  */
/*             - October 2012: added some wrapper functions for running CAMB with HyRec               */
/*                             (courtesy of Antony Lewis)                                             */
/*             - January 2015: - added DM annihilation and 21 cm routines.                            */
/*                             - changed cosmological parameters input form                           */
/*							   - possibility to change fine structure constant/ electron mass         */
/******************************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#include "hyrectools.h"
#include "helium.h"
#include "hydrogen.h"
#include "history.h"
#include "hyrec_21cm.h"

/*****************************************************************************
Setting derived cosmological parameters needed for recombination calculation
******************************************************************************/

void rec_set_derived_params(REC_COSMOPARAMS *param)
{
   /* Separated by AML to avoid duplication */

   double z, Pion, Tresc, RLya, four_betaB;

   param->nH0 = 11.223846333047 * param->obh2 * (1. - param->Y); /* number density of hydrogen today in m-3 */
   param->fHe = param->Y / (1 - param->Y) / 3.97153;             /* abundance of helium by number */
   /* these should depend on fsR and meR, strictly speaking; however, these are corrections to corrections */

   /* Total number of redshift steps */
   param->nz = (long)floor(2 + log((1. + ZSTART) / (1. + ZEND)) / DLNA);

   /* (Added May 2012)
      In order to save memory for the radiation field tables, compute right away
      the range of redshifts at which radiative transfer is actually followed,
      and the corresponding number of steps.
      Radiative transfer starts when Tr < TR_MAX and stops when the ionization probablity
      from n=2 (estimated with simple Peebles model) is less than PION_MAX. */

   param->izH0 = (long)floor(1 + log(kBoltz * param->T0 / square(param->fsR) / param->meR * (1. + ZSTART) / TR_MAX) / DLNA);
   param->zH0 = (1. + ZSTART) * exp(-param->izH0 * DLNA) - 1.;

   Pion = 1.;
   z = 900.;
   while (Pion > PION_MAX && z > ZEND)
   {
      Tresc = kBoltz * param->T0 * (1. + z) / param->fsR / param->fsR / param->meR; /* Tr rescaled for alpha, me */
      RLya = LYA_FACT(param->fsR, param->meR) * rec_HubbleConstant(param, z) / (1e-6 * param->nH0 * cube(1. + z));
      four_betaB = SAHA_FACT(param->fsR, param->meR) * Tresc * sqrt(Tresc) * exp(-0.25 * EI / Tresc) * alphaB_PPB(Tresc, param->fsR, param->meR);
      Pion = four_betaB / (3. * RLya + L2s_rescaled(param->fsR, param->meR) + four_betaB);
      z -= 10.;
   }

   param->nzrt = (long)floor(2 + log((1. + ZSTART) / (1. + z)) / DLNA) - param->izH0;
}

/********************************************************************************
Allocate and free memory for radiation field outputs (added December 2014)
********************************************************************************/

void allocate_rad_output(RAD_OUTPUTS *rad, long int nzrt)
{
   rad->Dfnu = create_2D_array(NVIRT, nzrt);
   rad->Dfminus_Ly[0] = create_1D_array(nzrt); /* Ly-alpha */
   rad->Dfminus_Ly[1] = create_1D_array(nzrt); /* Ly-beta  */
   rad->Dfminus_Ly[2] = create_1D_array(nzrt); /* Ly-gamma */
}

void free_rad_output(RAD_OUTPUTS *rad)
{
   free_2D_array(rad->Dfnu, NVIRT);
   free(rad->Dfminus_Ly[0]);
   free(rad->Dfminus_Ly[1]);
   free(rad->Dfminus_Ly[2]);
}

#ifdef CAMB

extern double dtauda_(double *);

/*************************************************************************************
Hubble expansion rate in sec^-1, obtained from CAMB.
*************************************************************************************/

double rec_HubbleConstant(REC_COSMOPARAMS *param, double z)
{
   double a;

   a = 1. / (1. + z);
   /* conversion from d tau/ da in Mpc to H(z) in 1/s */
   return 1. / (a * a) / dtauda_(&a) / 3.085678e22 * 2.99792458e8;
}

HRATEEFF rate_table;
TWO_PHOTON_PARAMS twog_params;

EFFECTIVE_EFF eff_ionization_annihilation;
EFFECTIVE_EFF eff_ionization_decay;
EFFECTIVE_EFF eff_ionization_blackhole;

EFFECTIVE_EFF eff_alpha_annihilation;
EFFECTIVE_EFF eff_alpha_decay;
EFFECTIVE_EFF eff_alpha_blackhole;

EFFECTIVE_EFF eff_thermal_annihilation;
EFFECTIVE_EFF eff_thermal_decay;
EFFECTIVE_EFF eff_thermal_blackhole;

REC_COSMOPARAMS param;
int firstTime = 0;
double *xe_output, *tm_output, *ts_output;
RAD_OUTPUTS rad;
double logstart;

void hyrec_init()
{

   char *buffer = (char *)malloc(1024);
   getcwd(buffer, 1024);
   chdir(HYRECPATH);

   /* Build effective rate table */
   read_rates(&rate_table);
   /* Read two-photon rate tables */
   read_twog_params(&twog_params);

   /* Read Effective Efficiency Tables */
   read_eff_eff(&eff_ionization_annihilation, ANNIHILATION_H_IONIZATION_ENERGY_AXIS_FILE, ANNIHILATION_H_IONIZATION_REDSHIFT_AXIS_FILE, ANNIHILATION_H_IONIZATION_EFF_TABLE_FILE, ANNIHILATION_ENERGY_AXIS_SIZE, ANNIHILATION_REDSHIFT_AXIS_SIZE);
   read_eff_eff(&eff_ionization_decay, DECAY_H_IONIZATION_ENERGY_AXIS_FILE, DECAY_H_IONIZATION_REDSHIFT_AXIS_FILE, DECAY_H_IONIZATION_EFF_TABLE_FILE, DECAY_ENERGY_AXIS_SIZE, DECAY_REDSHIFT_AXIS_SIZE);
   read_eff_eff(&eff_ionization_blackhole, BH_H_IONIZATION_ENERGY_AXIS_FILE, BH_H_IONIZATION_REDSHIFT_AXIS_FILE, BH_H_IONIZATION_EFF_TABLE_FILE, BH_ENERGY_AXIS_SIZE, BH_REDSHIFT_AXIS_SIZE);

   read_eff_eff(&eff_alpha_annihilation, ANNIHILATION_ALPHA_ENERGY_AXIS_FILE, ANNIHILATION_ALPHA_REDSHIFT_AXIS_FILE, ANNIHILATION_ALPHA_EFF_TABLE_FILE, ANNIHILATION_ENERGY_AXIS_SIZE, ANNIHILATION_REDSHIFT_AXIS_SIZE);
   read_eff_eff(&eff_alpha_decay, DECAY_ALPHA_ENERGY_AXIS_FILE, DECAY_ALPHA_REDSHIFT_AXIS_FILE, DECAY_ALPHA_EFF_TABLE_FILE, DECAY_ENERGY_AXIS_SIZE, DECAY_REDSHIFT_AXIS_SIZE);
   read_eff_eff(&eff_alpha_blackhole, BH_ALPHA_ENERGY_AXIS_FILE, BH_ALPHA_REDSHIFT_AXIS_FILE, BH_ALPHA_EFF_TABLE_FILE, BH_ENERGY_AXIS_SIZE, BH_REDSHIFT_AXIS_SIZE);

   read_eff_eff(&eff_thermal_annihilation, ANNIHILATION_HEATING_ENERGY_AXIS_FILE, ANNIHILATION_HEATING_REDSHIFT_AXIS_FILE, ANNIHILATION_HEATING_EFF_TABLE_FILE, ANNIHILATION_ENERGY_AXIS_SIZE, ANNIHILATION_REDSHIFT_AXIS_SIZE);
   read_eff_eff(&eff_thermal_decay, DECAY_HEATING_ENERGY_AXIS_FILE, DECAY_HEATING_REDSHIFT_AXIS_FILE, DECAY_HEATING_EFF_TABLE_FILE, DECAY_ENERGY_AXIS_SIZE, DECAY_REDSHIFT_AXIS_SIZE);
   read_eff_eff(&eff_thermal_blackhole, BH_HEATING_ENERGY_AXIS_FILE, BH_HEATING_REDSHIFT_AXIS_FILE, BH_HEATING_EFF_TABLE_FILE, BH_ENERGY_AXIS_SIZE, BH_REDSHIFT_AXIS_SIZE);

   chdir(buffer);
   free(buffer);
}

void rec_build_history_camb_(const double *OmegaC, const double *OmegaB, const double *OmegaN, const double *Omegav, const double *h0inp, const double *tcmb, const double *yp, const double *num_nu, const double *Pann, const double *DM_gamma, const double *DM_mass, const double *BH_mass, const double *DM_BH_frac, const double *fsR, const double *meR, const int *feedback)
{

   double h2 = *h0inp / 100.;
   h2 = h2 * h2;
   param.T0 = *tcmb;
   param.obh2 = *OmegaB * h2;
   param.omh2 = (*OmegaB + *OmegaC) * h2;
   param.okh2 = (1 - *OmegaC - *OmegaB - *Omegav - *OmegaN) * h2;
   param.odeh2 = *Omegav * h2;
   param.Y = *yp;
   param.Nnueff = *num_nu;
   param.fsR = *fsR;
   param.meR = *meR;               /*** Default: today's values = 1***/
   param.pann = *Pann * 1e-9;      /*** DM annihilation ***/
   param.dm_gamma = *DM_gamma;     /*** DM decay ***/
   param.dm_mass = *DM_mass * 1e9; /*** DM Mass ***/
   param.bh_mass = *BH_mass;       /*** BH Mass ***/
   param.dm_bh_frac = *DM_BH_frac; /*** DM BH Fraction ***/
   rec_set_derived_params(&param);
   if (*feedback > 0)
   {
      fprintf(stderr, " ********HYREC******************** \n");
      fprintf(stderr, " P_ann         = %e (in cm^3/s/eV) \n", param.pann);
      fprintf(stderr, " DM_Gamma      = %e (in 1/s) \n", param.dm_gamma);
      fprintf(stderr, " DM_Mass       = %e (in eV) \n", param.dm_mass);
      fprintf(stderr, " BH_Mass       = %e (in grams) \n", param.bh_mass);
      fprintf(stderr, " DM_BH_Frac    = %e \n", param.dm_bh_frac);
      fprintf(stderr, " m_e/m_e0      = %e \n", param.meR);
      fprintf(stderr, " alpha/alpha_0 = %e \n", param.fsR);
      fprintf(stderr, " ********HYREC******************** \n");
   }

   if (firstTime == 0)
   {
      hyrec_init();
      firstTime = 1;
      logstart = -log(1. + ZSTART);
      xe_output = create_1D_array(param.nz);
      tm_output = create_1D_array(param.nz);
   }
   allocate_rad_output(&rad, param.nzrt);

   rec_build_history(&param, &rate_table, &twog_params,
                     &eff_ionization_annihilation, &eff_ionization_decay, &eff_ionization_blackhole,
                     &eff_alpha_annihilation, &eff_alpha_decay, &eff_alpha_blackhole,
                     &eff_thermal_annihilation, &eff_thermal_decay, &eff_thermal_blackhole,
                     xe_output, tm_output, rad.Dfnu, rad.Dfminus_Ly);

   free_rad_output(&rad);
}

double hyrec_xe_(double *a)
{
   double loga = log(*a);
   if (loga < logstart)
      return xe_output[0];
   return rec_interp1d(logstart, DLNA, xe_output, param.nz, loga);
}

double hyrec_tm_(double *a)
{
   double loga = log(*a);
   if (loga < logstart)
      return tm_output[0];
   return rec_interp1d(logstart, DLNA, tm_output, param.nz, loga);
}

/* Added December 2014: 21 cm spin temperature, in Kelvins */
double hyrec_ts_(double *a)
{
   double loga = log(*a);
   if (loga < logstart)
      return param.T0 / (*a);
   double xe = rec_interp1d(logstart, DLNA, xe_output, param.nz, loga);
   double Tgas = rec_interp1d(logstart, DLNA, tm_output, param.nz, loga);
   return Tspin(xe, Tgas, param.T0 / (*a), 1e-6 * param.nH0 / cube(*a));
}

#else

/*************************************************************************************
Hubble expansion rate in sec^-1.
Note: neutrinos are assumed to be massless.
*************************************************************************************/

double rec_HubbleConstant(REC_COSMOPARAMS *param, double z)
{

   int j;
   double rho, rho_i[5], ainv;
   double ogh2;

   ainv = 1. + z; /* inverse scale factor */

   /* Matter */
   rho_i[0] = param->omh2 * ainv * ainv * ainv;

   /* Curvature */
   rho_i[1] = param->okh2 * ainv * ainv;

   /* Dark energy */
   rho_i[2] = param->odeh2;

   /* Get radiation density.
    * Coefficient based on 1 AU = 1.49597870691e11 m (JPL SSD) */
   ogh2 = 4.48162687719e-7 * param->T0 * param->T0 * param->T0 * param->T0;
   rho_i[3] = ogh2 * ainv * ainv * ainv * ainv;

   /* Neutrino density -- scaled from photons assuming lower temperature */
   rho_i[4] = 0.227107317660239 * rho_i[3] * param->Nnueff;

   /* Total density, including curvature contribution */
   rho = 0.;
   for (j = 0; j < 5; j++)
      rho += rho_i[j];

   /* Conversion to Hubble rate */
   return (3.2407792896393e-18 * sqrt(rho));
}
#endif

/*************************************************************************************************
Cosmological parameters Input/Output
*************************************************************************************************/

void rec_get_cosmoparam(FILE *fin, REC_COSMOPARAMS *param, REC_EXTRAS *extras)
{
   double h, Omega_b, Omega_m, Omega_L, Omega_K;
   char yorn[3];

   /* Cosmological parameters.
      January 2015: Changed the input form and removed equation of state params. */

   fscanf(fin, "%lg", &(param->T0));
   fscanf(fin, "%lg", &h);
   fscanf(fin, "%lg", &Omega_b);
   fscanf(fin, "%lg", &Omega_m);
   fscanf(fin, "%lg", &Omega_L);
   fscanf(fin, "%lg", &(param->Y));
   fscanf(fin, "%lg", &(param->Nnueff));

   Omega_K = 1. - (Omega_m + Omega_L);
   param->obh2 = Omega_b * h * h;
   param->omh2 = Omega_m * h * h;
   param->odeh2 = Omega_L * h * h;
   param->okh2 = Omega_K * h * h;

   /****** Added May 2012: explicit dependence on fine-structure constant and electron mass ******/
   /** fsR = alpha_fs(rec) / alpha_fs(today), meR = me(rec) / me(today) **/

   fscanf(fin, "%lg", &(param->fsR));
   fscanf(fin, "%lg", &(param->meR));

   /**** Added January 2015: dark matter annihilation parameter pann ****/
   fscanf(fin, "%lg", &(param->pann));
   /* Convert from cm^3/s/GeV to cm^3/s/eV */
   param->pann *= 1e-9;

   /**** Added February 2016: dark matter decay parameter dm_gamma ****/
   fscanf(fin, "%lg", &(param->dm_gamma));
   /* dm_gamma input in 1/s */

   /**** Added September 2016: dark matter mass dm_mass ****/
   fscanf(fin, "%lg", &(param->dm_mass));
   /*Convert from GeV to eV*/
   param->dm_mass *= 1e9;

   /**** Added September 2016: black hole mass bh_mass ****/
   fscanf(fin, "%lg", &(param->bh_mass));
   /* bh_mass input in grams */

   /**** Added September 2016: dark matter black hole dark fraction dm_bh_frac ****/
   fscanf(fin, "%lg", &(param->dm_bh_frac));
   /*dm_bh_frac = 0 for no black holes*/

   rec_set_derived_params(param);

   /**** Added January 2015: switches to print the Lyman-lines spectrum
         and/or 21cm temperatures ****/

   fscanf(fin, "%s", yorn);
   if (yorn[0] == 'y' || yorn[0] == 'Y')
      extras->print_spec = 1;
   else
      extras->print_spec = 0;

   fscanf(fin, "%s", extras->file_spec);
   fscanf(fin, "%le", &(extras->zmin_spec));
   fscanf(fin, "%le", &(extras->zmax_spec));
   fscanf(fin, "%d", &(extras->Nz_spec));

   fscanf(fin, "%s", yorn);
   if (yorn[0] == 'y' || yorn[0] == 'Y')
      extras->print_21cm = 1;
   else
      extras->print_21cm = 0;
   fscanf(fin, "%s", extras->file_21cm);
   fscanf(fin, "%le", &(extras->zmin_21cm));
   fscanf(fin, "%le", &(extras->zmax_21cm));
   fscanf(fin, "%d", &(extras->Nz_21cm));
}

/****************************************************************************************
Total volumic rate of energy realase due to DM annihilation, in eV/cm^3/s.
Can replace by any function (just keep the same name).
Added December 2014
Added DM decay February 2016
****************************************************************************************/
/*double dE_dtdV_func(double z, REC_COSMOPARAMS *param){
 return  (square(10537.4*(param->omh2-param->obh2)*cube(1.+z)) * param->pann)+(10537.4*(param->omh2-param->obh2)*cube(1.+z) * param->dm_gamma);  */
/* the prefactor is 3 H0^2/(8 Pi G) c^2 in eV/cm^3, H0 = 100km/s/Mpc */
/*}*/ /*Moved to bottom of file*/

/*****************************************************************************************
Matter temperature -- 1st order steady state, from Hirata 2008.
The input and output temperatures are in KELVIN.
Added December 2014: possibility for additional energy injection.
******************************************************************************************/

double rec_Tmss(double xe, double Tr, double H, double fHe,
                double fsR, double meR, double dE_dtdV, double nh)
{

   double coeff = fsR * fsR / meR / meR / meR * 4.91466895548409e-22 * Tr * Tr * Tr * Tr * xe / (1. + xe + fHe) / H;

   return Tr / (1. + 1. / coeff) + dE_dtdV / kBoltz / (1.5 * nh * (1. + xe + fHe)) / H / (1. + coeff);

   /* Coefficient = 8 sigma_T a_r / (3 m_e c) */
   /* Here Tr, Tm are the actual (not rescaled) temperatures */
}

/******************************************************************************************
Matter temperature evolution derivative. Input and output temperatures are in KELVIN.
Added May 2012: when Tm = Tr, return -Tr (needed by CLASS)
Added December 2014: possibility of additional energy injection.
******************************************************************************************/

double rec_dTmdlna(double xe, double Tm, double Tr, double H, double fHe,
                   double fsR, double meR, double dE_dtdV, double nh)
{
   return (/*Tr/Tm-1.<1e-10 ? -Tr : */
           -2. * Tm + fsR * fsR / meR / meR / meR * 4.91466895548409e-22 * Tr * Tr * Tr * Tr * xe / (1. + xe + fHe) * (Tr - Tm) / H + dE_dtdV / kBoltz / (1.5 * nh * (1. + xe + fHe)) / H);
   /* Coefficient = 8 sigma_T a_r / (3 m_e c) */
   /* Here Tr, Tm are the actual (not rescaled) temperatures */
}

/**********************************************************************************************
Second-order integrator for HeII abundance during HeII->HeI recombination.
Assumes hydrogen ionization is in Saha equilibrium (free electrons being provided by both H and He).
If xHeII is close enough to Saha equilibrium do a post-Saha expansion.
May 2012: removed unused z_prev, z_prev2 variables
***********************************************************************************************/

void rec_get_xe_next1_He(REC_COSMOPARAMS *param, double z_in, double *xHeII,
                         double *dxHeIIdlna_prev, double *dxHeIIdlna_prev2, int *post_saha)
{

   double H, xH1s, xH1s_p, xH1s_m, xHeIISaha, dxHeIISaha_dlna, DdxHeIIdlna_Dxe, dxHeIIdlna, z_out, Dxe;

   H = rec_HubbleConstant(param, z_in);
   xH1s = rec_saha_xH1s(*xHeII, param->nH0, param->T0, z_in, param->fsR, param->meR);
   dxHeIIdlna = rec_helium_dxHeIIdlna(xH1s, *xHeII, param->nH0, param->T0, param->fHe, H, z_in, param->fsR, param->meR);

   /* Post-Saha approximation during the early phase of HeII->HeI recombination */
   if (*post_saha == 1)
   {
      z_out = (1. + z_in) * exp(-DLNA) - 1.;
      H = rec_HubbleConstant(param, z_out);
      xHeIISaha = rec_saha_xHeII(param->nH0, param->T0, param->fHe, z_out, param->fsR, param->meR);

      dxHeIISaha_dlna = (1. + z_out) * (rec_saha_xHeII(param->nH0, param->T0, param->fHe, z_out - 0.5, param->fsR, param->meR) - rec_saha_xHeII(param->nH0, param->T0, param->fHe, z_out + 0.5, param->fsR, param->meR));

      Dxe = 0.01 * (param->fHe - xHeIISaha);
      xH1s_p = rec_saha_xH1s(xHeIISaha + Dxe, param->nH0, param->T0, z_out, param->fsR, param->meR);
      xH1s_m = rec_saha_xH1s(xHeIISaha - Dxe, param->nH0, param->T0, z_out, param->fsR, param->meR);

      DdxHeIIdlna_Dxe = (rec_helium_dxHeIIdlna(xH1s_p, xHeIISaha + Dxe, param->nH0, param->T0, param->fHe, H, z_out, param->fsR, param->meR) - rec_helium_dxHeIIdlna(xH1s_m, xHeIISaha - Dxe, param->nH0, param->T0, param->fHe, H, z_out, param->fsR, param->meR)) / 2. / Dxe;

      *xHeII = xHeIISaha + dxHeIISaha_dlna / DdxHeIIdlna_Dxe;

      /* Check that the post-Saha approximation is still valid. If not, switch it off for future iterations */
      if (fabs(*xHeII - xHeIISaha) > DXHEII_MAX)
         *post_saha = 0;
   }

   /* Otherwise integrate ODE */
   else
      *xHeII += DLNA * (1.25 * dxHeIIdlna - 0.25 * (*dxHeIIdlna_prev2));

   /* Update stored derivatives */
   *dxHeIIdlna_prev2 = *dxHeIIdlna_prev;
   *dxHeIIdlna_prev = dxHeIIdlna;
}

/****************************************************************************************************
Post-Saha value of neutral hydrogen abundance at the redshift z_out.
Tm = Tr is assumed.
xHeII_out is the given value of Helium abundance at z_out (may be non-zero if Helium is still recombining)
iz_out is the index corresponding to z_out.
****************************************************************************************************/

double rec_xH1s_postSaha(REC_COSMOPARAMS *param, unsigned iz_out, double z_out, double xHeII_out,
                         HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog_params,
                         double **Dfminus_hist, double *Dfminus_Ly_hist[], double **Dfnu_hist, int *post_saha,
                         EFFECTIVE_EFF *ionization_annihilation, EFFECTIVE_EFF *ionization_decay, EFFECTIVE_EFF *ionization_blackhole,
                         EFFECTIVE_EFF *alpha_annihilation, EFFECTIVE_EFF *alpha_decay, EFFECTIVE_EFF *alpha_blackhole)
{

   double ainv, xH1sSaha, xHIISaha, dxH1sSaha_dlna, dxH1sdlna_Saha, DdxH1sdlna_DxH1s, H, T, nH, Dxe, xH1s;
   double ionization_dE_dtdV = dE_dtdV_func(z_out, param, ionization_annihilation, ionization_decay, ionization_blackhole);
   double alpha_dE_dtdV = dE_dtdV_func(z_out, param, alpha_annihilation, alpha_decay, alpha_blackhole);

   xH1sSaha = rec_saha_xH1s(xHeII_out, param->nH0, param->T0, z_out, param->fsR, param->meR);
   xHIISaha = 1. - xH1sSaha;
   H = rec_HubbleConstant(param, z_out);
   T = kBoltz * param->T0 * (ainv = 1. + z_out); /* Convert to eV for hydrogen rec functions */
   nH = 1e-6 * param->nH0 * ainv * ainv * ainv;  /* Convert to cm^-3 for hydrogen rec functions */

   dxH1sSaha_dlna = (1. + z_out) * (rec_saha_xH1s(xHeII_out, param->nH0, param->T0, z_out - 0.5, param->fsR, param->meR) - rec_saha_xH1s(xHeII_out, param->nH0, param->T0, z_out + 0.5, param->fsR, param->meR));
   /* (partial xHII)/(partial lna). Use xH1s = 1-xHII for better accuracy. */
   if (xHeII_out != 0.)
   {
      Dxe = 0.01 * (param->fHe - xHeII_out);
      dxH1sSaha_dlna += (rec_saha_xH1s(xHeII_out + Dxe, param->nH0, param->T0, z_out, param->fsR, param->meR) - rec_saha_xH1s(xHeII_out - Dxe, param->nH0, param->T0, z_out, param->fsR, param->meR)) / 2. / Dxe * rec_helium_dxHeIIdlna(xH1sSaha, xHeII_out, param->nH0, param->T0, param->fHe, H, z_out, param->fsR, param->meR);
      /* (partial xHII)/(partial xHeII).dxHeII/dlna */
   }
   dxH1sdlna_Saha = -rec_dxHIIdlna(MODEL, xHIISaha + xHeII_out, xHIISaha, nH, H, T, T, rate_table, twog_params,
                                   Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist, param->zH0, iz_out - param->izH0, z_out,
                                   param->fsR, param->meR, ionization_dE_dtdV, alpha_dE_dtdV);
   Dxe = 0.01 * xH1sSaha;
   DdxH1sdlna_DxH1s = (rec_dxHIIdlna(MODEL, xHIISaha + Dxe + xHeII_out, xHIISaha + Dxe, nH, H, T, T, rate_table, twog_params,
                                     Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist, param->zH0, iz_out - param->izH0, z_out,
                                     param->fsR, param->meR, ionization_dE_dtdV, alpha_dE_dtdV) -
                       rec_dxHIIdlna(MODEL, xHIISaha - Dxe + xHeII_out, xHIISaha - Dxe, nH, H, T, T, rate_table, twog_params,
                                     Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist, param->zH0, iz_out - param->izH0, z_out,
                                     param->fsR, param->meR, ionization_dE_dtdV, alpha_dE_dtdV)) /
                      2. / Dxe;

   xH1s = xH1sSaha + (dxH1sSaha_dlna - dxH1sdlna_Saha) / DdxH1sdlna_DxH1s;

   /* Check that we are still close enough to Saha equilibrium. If not, switch post-saha expansion off */
   if (fabs(xH1s - xH1sSaha) > DXHII_MAX)
      *post_saha = 0;

   return xH1s;
}

/*****************************************************************************************************
Second-order integrator used to evolve simultaneously H(1s) and HeII.
When hydrogen is close to Saha equilibrium but there is still a significant amount of HeII,
use a post-Saha expansion for hydrogen. The other way around is not treated (simply integrate H and He until
there is almost no HeII left, then integrate H only)
******************************************************************************************************/

void get_rec_next2_HHe(REC_COSMOPARAMS *param, unsigned iz_in, double z_in, double Tm_in, double *xH1s, double *xHeII,
                       HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog_params, double **Dfminus_hist, double *Dfminus_Ly_hist[], double **Dfnu_hist,
                       double *dxHIIdlna_prev, double *dxHeIIdlna_prev, double *dxHIIdlna_prev2, double *dxHeIIdlna_prev2, int *post_saha,
                       EFFECTIVE_EFF *ionization_annihilation, EFFECTIVE_EFF *ionization_decay, EFFECTIVE_EFF *ionization_blackhole,
                       EFFECTIVE_EFF *alpha_annihilation, EFFECTIVE_EFF *alpha_decay, EFFECTIVE_EFF *alpha_blackhole)
{

   double H, dxHeIIdlna, dxHIIdlna, TR, TM, nH, zout, ainv, xH1s_in, xHeII_in, xe_in;
   double ionization_dE_dtdV = dE_dtdV_func(z_in, param, ionization_annihilation, ionization_decay, ionization_blackhole);
   double alpha_dE_dtdV = dE_dtdV_func(z_in, param, alpha_annihilation, alpha_decay, alpha_blackhole);

   xH1s_in = *xH1s;
   xHeII_in = *xHeII;
   xe_in = xHeII_in + (1. - xH1s_in);

   /* Evolve HeII by solving ODE */
   H = rec_HubbleConstant(param, z_in);
   dxHeIIdlna = rec_helium_dxHeIIdlna(xH1s_in, xHeII_in, param->nH0, param->T0, param->fHe, H, z_in, param->fsR, param->meR);
   *xHeII += DLNA * (1.25 * dxHeIIdlna - 0.25 * (*dxHeIIdlna_prev2));

   /* Compute Hydrogen derivative at input time. Even if using the post-Saha expansion, needed for getting the correct radiation field at z_in */
   TR = kBoltz * param->T0 * (ainv = 1. + z_in);
   TM = kBoltz * Tm_in;
   nH = 1e-6 * param->nH0 * ainv * ainv * ainv;
   dxHIIdlna = rec_dxHIIdlna(MODEL, xe_in, (1. - xH1s_in), nH, H, TM, TR, rate_table, twog_params,
                             Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist, param->zH0, iz_in - param->izH0, z_in,
                             param->fsR, param->meR, ionization_dE_dtdV, alpha_dE_dtdV);

   /* If Hydrogen is still close to Saha equilibrium do a post-Saha expansion for Hydrogen */
   if (*post_saha == 1)
   {
      zout = (1. + z_in) * exp(-DLNA) - 1.; /* Redshift for the output */
      *xH1s = rec_xH1s_postSaha(param, iz_in + 1, zout, *xHeII, rate_table, twog_params,
                                Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist, post_saha,
                                ionization_annihilation, ionization_decay, ionization_blackhole,
                                alpha_annihilation, alpha_decay, alpha_blackhole);
   }

   /* Otherwise solve HII ODE */
   else
      *xH1s -= DLNA * (1.25 * dxHIIdlna - 0.25 * (*dxHIIdlna_prev2));

   /* Update derivatives */
   *dxHIIdlna_prev2 = *dxHIIdlna_prev;
   *dxHIIdlna_prev = dxHIIdlna;
   *dxHeIIdlna_prev2 = *dxHeIIdlna_prev;
   *dxHeIIdlna_prev = dxHeIIdlna;
}

/*********************************************************************************************************
Second-order integrator to evolve hydrogen only, assuming helium has entirely recombined.
Tm is given as an input (to avoid computing it twice) and fixed to quasi-equilibrium value with Tr.
**********************************************************************************************************/

void rec_get_xe_next1_H(REC_COSMOPARAMS *param, double z_in, double xe_in, double Tm_in, double *xe_out,
                        HRATEEFF *rate_table, unsigned iz, TWO_PHOTON_PARAMS *twog_params,
                        double **Dfminus_hist, double *Dfminus_Ly_hist[], double **Dfnu_hist,
                        double *dxedlna_prev, double *dxedlna_prev2, int *post_saha,
                        EFFECTIVE_EFF *ionization_annihilation, EFFECTIVE_EFF *ionization_decay, EFFECTIVE_EFF *ionization_blackhole,
                        EFFECTIVE_EFF *alpha_annihilation, EFFECTIVE_EFF *alpha_decay, EFFECTIVE_EFF *alpha_blackhole)
{

   double dxedlna, TR, nH, ainv, H, TM, zout;
   int model;
   double ionization_dE_dtdV = dE_dtdV_func(z_in, param, ionization_annihilation, ionization_decay, ionization_blackhole);
   double alpha_dE_dtdV = dE_dtdV_func(z_in, param, alpha_annihilation, alpha_decay, alpha_blackhole);

   TR = kBoltz * param->T0 * (ainv = 1. + z_in);
   nH = 1e-6 * param->nH0 * ainv * ainv * ainv;
   H = rec_HubbleConstant(param, z_in);
   TM = kBoltz * Tm_in;

   /* Switch off radiative transfer calculation if needed (param->nzrt and izH0 are evaluated in rec_get_cosmoparam) */
   model = (iz - param->izH0 < param->nzrt || MODEL != FULL) ? MODEL : EMLA2s2p;

   dxedlna = rec_dxHIIdlna(model, xe_in, xe_in, nH, H, TM, TR, rate_table, twog_params,
                           Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist, param->zH0, iz - param->izH0, z_in,
                           param->fsR, param->meR, ionization_dE_dtdV, alpha_dE_dtdV);

   /* If close to Saha equilibrium (with xHeII = 0), do a post-Saha expansion */
   if (*post_saha == 1)
   {
      zout = (1. + z_in) * exp(-DLNA) - 1.; /* Redshift for the output */
      *xe_out = 1. - rec_xH1s_postSaha(param, iz + 1, zout, 0., rate_table, twog_params,
                                       Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist, post_saha,
                                       ionization_annihilation, ionization_decay, ionization_blackhole,
                                       alpha_annihilation, alpha_decay, alpha_blackhole);
   }

   /* Otherwise evolve ODE */
   else
      *xe_out = xe_in + DLNA * (1.25 * dxedlna - 0.25 * (*dxedlna_prev2));

   /* Update previous derivatives */
   *dxedlna_prev2 = *dxedlna_prev;
   *dxedlna_prev = dxedlna;
}

/**********************************************************************************************
Second-order integrator for evolving xe and Tm simultaneously
Used for Hydrogen recombination only
May 2012: added a switch so Peebles model can be used at low redshift.
***********************************************************************************************/

void rec_get_xe_next2_HTm(int func_select, REC_COSMOPARAMS *param, double z_in, double xe_in, double Tm_in, double *xe_out, double *Tm_out,
                          HRATEEFF *rate_table, unsigned iz, TWO_PHOTON_PARAMS *twog_params, double **Dfminus_hist, double *Dfminus_Ly_hist[],
                          double **Dfnu_hist, double *dxedlna_prev, double *dTmdlna_prev, double *dxedlna_prev2, double *dTmdlna_prev2,
                          EFFECTIVE_EFF *ionization_annihilation, EFFECTIVE_EFF *ionization_decay, EFFECTIVE_EFF *ionization_blackhole,
                          EFFECTIVE_EFF *alpha_annihilation, EFFECTIVE_EFF *alpha_decay, EFFECTIVE_EFF *alpha_blackhole,
                          EFFECTIVE_EFF *thermal_annihilation, EFFECTIVE_EFF *thermal_decay, EFFECTIVE_EFF *thermal_blackhole)
{

   double dxedlna, dTmdlna, TR, nH, ainv, H, TM;
   int model;
   double ionization_dE_dtdV = dE_dtdV_func(z_in, param, ionization_annihilation, ionization_decay, ionization_blackhole);
   double alpha_dE_dtdV = dE_dtdV_func(z_in, param, alpha_annihilation, alpha_decay, alpha_blackhole);
   double thermal_dE_dtdV = dE_dtdV_func(z_in, param, thermal_annihilation, thermal_decay, thermal_blackhole);

   TR = kBoltz * param->T0 * (ainv = 1. + z_in);
   nH = 1e-6 * param->nH0 * ainv * ainv * ainv;
   H = rec_HubbleConstant(param, z_in);
   TM = kBoltz * Tm_in;

   /* Switch off radiative transfer calculation if needed */
   model = (iz - param->izH0 < param->nzrt || func_select != FULL) ? func_select : EMLA2s2p;

   dxedlna = rec_dxHIIdlna(model, xe_in, xe_in, nH, H, TM, TR, rate_table, twog_params,
                           Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist, param->zH0, iz - param->izH0, z_in,
                           param->fsR, param->meR, ionization_dE_dtdV, alpha_dE_dtdV);

   dTmdlna = rec_dTmdlna(xe_in, Tm_in, TR / kBoltz, H, param->fHe,
                         param->fsR, param->meR, thermal_dE_dtdV, nH);

   *xe_out = xe_in + DLNA * (1.25 * dxedlna - 0.25 * (*dxedlna_prev2));
   *Tm_out = Tm_in + DLNA * (1.25 * dTmdlna - 0.25 * (*dTmdlna_prev2));

   *dxedlna_prev2 = *dxedlna_prev;
   *dTmdlna_prev2 = *dTmdlna_prev;
   *dxedlna_prev = dxedlna;
   *dTmdlna_prev = dTmdlna;
}

/****************************************************************************************************
Builds a recombination history
Added May 2012: The radiation field was added as an input so it can be extracted if desired
****************************************************************************************************/

void rec_build_history(REC_COSMOPARAMS *param, HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog_params,
                       EFFECTIVE_EFF *eff_ionization_annihilation, EFFECTIVE_EFF *eff_ionization_decay, EFFECTIVE_EFF *eff_ionization_blackhole,
                       EFFECTIVE_EFF *eff_alpha_annihilation, EFFECTIVE_EFF *eff_alpha_decay, EFFECTIVE_EFF *eff_alpha_blackhole,
                       EFFECTIVE_EFF *eff_thermal_annihilation, EFFECTIVE_EFF *eff_thermal_decay, EFFECTIVE_EFF *eff_thermal_blackhole,
                       double *xe_output, double *Tm_output, double **Dfnu_hist, double *Dfminus_Ly_hist[3])
{

   long iz;
   double z, dxHIIdlna_prev, dxHIIdlna_prev2, dTmdlna_prev, dTmdlna_prev2, dxHeIIdlna_prev, dxHeIIdlna_prev2;
   double Delta_xe, xHeII, xH1s;
   double **Dfminus_hist;
   int post_saha;

   Dfminus_hist = create_2D_array(NVIRT, param->nzrt);

   /* Make sure the input spectrum is initialized at zero */
   for (iz = 0; iz < param->nzrt; iz++)
      Dfminus_Ly_hist[0][iz] = Dfminus_Ly_hist[1][iz] = Dfminus_Ly_hist[2][iz] = 0;

   z = ZSTART;

   /********* He III -> II Saha phase. Tm = Tr. Stop when xHeIII = 1e-8 *********/
   Delta_xe = param->fHe; /* Delta_xe = xHeIII here */

   for (iz = 0; iz < param->nz && Delta_xe > 1e-8; iz++)
   {
      z = (1. + ZSTART) * exp(-DLNA * iz) - 1.;
      xe_output[iz] = rec_xesaha_HeII_III(param->nH0, param->T0, param->fHe, z, &Delta_xe, param->fsR, param->meR);
      Tm_output[iz] = param->T0 * (1. + z);
   }

   /******** He II -> I recombination.
             Hydrogen in Saha equilibrium with the free electrons.
             Tm fixed to steady state.
             Integrate until TR is low enough that can start integrating hydrogen recombination
             (this occurs at index izH0 computed in rec_get_cosmoparam).
             Start with post-Saha expansion.
    ********/

   dxHeIIdlna_prev2 = (xe_output[iz - 2] - xe_output[iz - 4]) / 2. / DLNA;
   dxHeIIdlna_prev = (xe_output[iz - 1] - xe_output[iz - 3]) / 2. / DLNA;

   xHeII = rec_saha_xHeII(param->nH0, param->T0, param->fHe, z, param->fsR, param->meR);
   post_saha = 1; /* Start with post-saha expansion */

   for (; iz < param->izH0 + 1; iz++)
   {
      rec_get_xe_next1_He(param, z, &xHeII, &dxHeIIdlna_prev, &dxHeIIdlna_prev2, &post_saha);
      z = (1. + ZSTART) * exp(-DLNA * iz) - 1.;
      xH1s = rec_saha_xH1s(xHeII, param->nH0, param->T0, z, param->fsR, param->meR);
      xe_output[iz] = (1. - xH1s) + xHeII;
      Tm_output[iz] = rec_Tmss(xe_output[iz], param->T0 * (1. + z), rec_HubbleConstant(param, z), param->fHe,
                               param->fsR, param->meR, dE_dtdV_func(z, param, eff_thermal_annihilation, eff_thermal_decay, eff_thermal_blackhole), 1e-6 * param->nH0 * cube(1. + z));
   }

   /******** H II -> I and He II -> I simultaneous recombination (rarely needed but just in case)
             Tm fixed to steady state.
             Integrate H and He simultaneously until xHeII < XHEII_MIN
             Start with post-saha expansion for hydrogen
    ********/

   dxHIIdlna_prev2 = (xe_output[iz - 2] - xe_output[iz - 4]) / 2. / DLNA - dxHeIIdlna_prev2;
   dxHIIdlna_prev = (xe_output[iz - 1] - xe_output[iz - 3]) / 2. / DLNA - dxHeIIdlna_prev;
   post_saha = 1;

   for (; iz < param->nz && xHeII > XHEII_MIN; iz++)
   {
      get_rec_next2_HHe(param, iz - 1, z, Tm_output[iz - 1], &xH1s, &xHeII, rate_table, twog_params, Dfminus_hist, Dfminus_Ly_hist,
                        Dfnu_hist, &dxHIIdlna_prev, &dxHeIIdlna_prev, &dxHIIdlna_prev2, &dxHeIIdlna_prev2, &post_saha,
                        eff_ionization_annihilation, eff_ionization_decay, eff_ionization_blackhole,
                        eff_alpha_annihilation, eff_alpha_decay, eff_alpha_blackhole);
      xe_output[iz] = (1. - xH1s) + xHeII;
      z = (1. + ZSTART) * exp(-DLNA * iz) - 1.;
      Tm_output[iz] = rec_Tmss(xe_output[iz], param->T0 * (1. + z), rec_HubbleConstant(param, z), param->fHe,
                               param->fsR, param->meR, dE_dtdV_func(z, param, eff_thermal_annihilation, eff_thermal_decay, eff_thermal_blackhole), 1e-6 * param->nH0 * cube(1. + z));
   }

   /******** H recombination. Helium assumed entirely neutral.
             Tm fixed to steady-state until its relative difference from Tr is DLNT_MAX
   ********/

   for (; iz < param->nz && 1. - Tm_output[iz - 1] / param->T0 / (1. + z) < DLNT_MAX; iz++)
   {
      rec_get_xe_next1_H(param, z, xe_output[iz - 1], Tm_output[iz - 1], xe_output + iz, rate_table, iz - 1, twog_params,
                         Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist, &dxHIIdlna_prev, &dxHIIdlna_prev2, &post_saha,
                         eff_ionization_annihilation, eff_ionization_decay, eff_ionization_blackhole,
                         eff_alpha_annihilation, eff_alpha_decay, eff_alpha_blackhole);
      z = (1. + ZSTART) * exp(-DLNA * iz) - 1.;
      Tm_output[iz] = rec_Tmss(xe_output[iz], param->T0 * (1. + z), rec_HubbleConstant(param, z), param->fHe,
                               param->fsR, param->meR, dE_dtdV_func(z, param, eff_thermal_annihilation, eff_thermal_decay, eff_thermal_blackhole), 1e-6 * param->nH0 * cube(1. + z));
   }

   /******** Evolve xe and Tm simultaneously until the lower bounds of integration tables are reached.
             Note that the radiative transfer calculation is switched off automatically in the functions
             rec_get_xe_next1_H and rec_get_xe_next2_HTm when it is no longer relevant.
   ********/

   dTmdlna_prev2 = (Tm_output[iz - 2] - Tm_output[iz - 4]) / 2. / DLNA;
   dTmdlna_prev = (Tm_output[iz - 1] - Tm_output[iz - 3]) / 2. / DLNA;

   for (; iz < param->nz && kBoltz * param->T0 * (1. + z) / param->fsR / param->fsR / param->meR > TR_MIN && Tm_output[iz - 1] / param->T0 / (1. + z) > TM_TR_MIN; iz++)
   {
      rec_get_xe_next2_HTm(MODEL, param, z, xe_output[iz - 1], Tm_output[iz - 1], xe_output + iz, Tm_output + iz,
                           rate_table, iz - 1, twog_params, Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist,
                           &dxHIIdlna_prev, &dTmdlna_prev, &dxHIIdlna_prev2, &dTmdlna_prev2,
                           eff_ionization_annihilation, eff_ionization_decay, eff_ionization_blackhole,
                           eff_alpha_annihilation, eff_alpha_decay, eff_alpha_blackhole,
                           eff_thermal_annihilation, eff_thermal_decay, eff_thermal_blackhole);
      z = (1. + ZSTART) * exp(-DLNA * iz) - 1.;
   }

   /***** For low redshifts (z < 20 or so) use Peeble's model (Tm is evolved with xe).
           The precise model does not metter much here as
           1) the free electron fraction is basically zero (~1e-4) in any case and
           2) the universe is going to be reionized around that epoch
    *****/

   for (; iz < param->nz; iz++)
   {
      rec_get_xe_next2_HTm(PEEBLES, param, z, xe_output[iz - 1], Tm_output[iz - 1], xe_output + iz, Tm_output + iz,
                           rate_table, iz - 1, twog_params, Dfminus_hist, Dfminus_Ly_hist, Dfnu_hist,
                           &dxHIIdlna_prev, &dTmdlna_prev, &dxHIIdlna_prev2, &dTmdlna_prev2,
                           eff_ionization_annihilation, eff_ionization_decay, eff_ionization_blackhole,
                           eff_alpha_annihilation, eff_alpha_decay, eff_alpha_blackhole,
                           eff_thermal_annihilation, eff_thermal_decay, eff_thermal_blackhole);
      z = (1. + ZSTART) * exp(-DLNA * iz) - 1.;
   }

   /* Cleanup */
   free_2D_array(Dfminus_hist, NVIRT);
}
/******************************************************************************************************************************/

/** Question: has the effect of DM annihilation on Helium recombination been properly accounted for?
    Since dE/dtdV ~ (1+z)^6 it could be that it has a large effect on Helium recombination?
    Also, is Hydrogen still in Saha equilibrium at high-z and does this affect the absorption of hydrogen on Helium lyman alpha line?
**/

/****************************************************************************************/
/**************************************************************************
Reads effective efficiency table data.
***************************************************************************/

void read_eff_eff(EFFECTIVE_EFF *effeff, char *ENERGY_AXIS_FILE, char *REDSHIFT_AXIS_FILE, char *EFFICIENCY_FILE, unsigned ENERGY_AXIS_SIZE, unsigned REDSHIFT_AXIS_SIZE)
{

   unsigned i, j;

   FILE *FEeff = fopen(ENERGY_AXIS_FILE, "r");
   effeff->energy_axis = create_1D_array(ENERGY_AXIS_SIZE);
   for (i = 0; i < ENERGY_AXIS_SIZE; i++)
   {
      fscanf(FEeff, "%le", &(effeff->energy_axis[i]));
   }
   fclose(FEeff);

   FILE *FReff = fopen(REDSHIFT_AXIS_FILE, "r");
   effeff->redshift_axis = create_1D_array(REDSHIFT_AXIS_SIZE);
   for (i = 0; i < REDSHIFT_AXIS_SIZE; i++)
   {
      fscanf(FReff, "%le", &(effeff->redshift_axis[i]));
   }
   fclose(FReff);

   FILE *FEFFeff = fopen(EFFICIENCY_FILE, "r");
   effeff->efficiency = create_2D_array(REDSHIFT_AXIS_SIZE, ENERGY_AXIS_SIZE);
   for (i = 0; i < REDSHIFT_AXIS_SIZE; i++)
   {
      for (j = 0; j < ENERGY_AXIS_SIZE; j++)
      {
         fscanf(FEFFeff, "%le", &(effeff->efficiency[i][j]));
      }
   }
   fclose(FEFFeff);
}

/**************************************************************************
Interpolation for Effective Efficiency.
***************************************************************************/

double interp_eff_eff(EFFECTIVE_EFF *effeff, double energy, double redshift, unsigned ENERGY_AXIS_SIZE, unsigned REDSHIFT_AXIS_SIZE)
{

   /*Find Energy Near Points*/
   unsigned count, step, iteration, index, bin_start_location, first_point, energy_pt1, energy_pt2, redshift_pt1, redshift_pt2;
   double MinE, MaxE, MinZp, MaxZp;
   MinZp = effeff->redshift_axis[0];
   MaxZp = effeff->redshift_axis[REDSHIFT_AXIS_SIZE - 1];
   MinE = effeff->energy_axis[0];
   MaxE = effeff->energy_axis[ENERGY_AXIS_SIZE - 1];

   // double energy_value, redshift_value;
   count = ENERGY_AXIS_SIZE;
   bin_start_location = 0;
   while (count > 0)
   {
      iteration = bin_start_location;
      step = count / 2;
      index = iteration + step;
      if (effeff->energy_axis[index] < energy)
      {
         bin_start_location = index + 1;
         count -= step + 1;
      }
      else
      {
         count = step;
      }
   }
   if (bin_start_location == 0)
   {
      energy_pt1 = 0;
      energy_pt2 = 1;
   }
   else if (bin_start_location == ENERGY_AXIS_SIZE)
   {
      energy_pt1 = bin_start_location - 2;
      energy_pt2 = bin_start_location - 1;
   }
   else
   {
      energy_pt1 = bin_start_location - 1;
      energy_pt2 = bin_start_location;
   }

   /*Find Redshift Near Points*/
   count = REDSHIFT_AXIS_SIZE;
   bin_start_location = 0;
   while (count > 0)
   {
      iteration = bin_start_location;
      step = count / 2;
      index = iteration + step;
      if (effeff->redshift_axis[index] < redshift)
      {
         bin_start_location = index + 1;
         count -= step + 1;
      }
      else
      {
         count = step;
      }
   }
   if (bin_start_location == 0)
   {
      redshift_pt1 = 0;
      redshift_pt2 = 1;
   }
   else if (bin_start_location == REDSHIFT_AXIS_SIZE)
   {
      redshift_pt1 = bin_start_location - 2;
      redshift_pt2 = bin_start_location - 1;
   }
   else
   {
      redshift_pt1 = bin_start_location - 1;
      redshift_pt2 = bin_start_location;
   }

   /*Interpolates the Effective Efficiency*/
   // Use linear method in log axis
   double x, x1, x2, y, y1, y2, f11, f12, f21, f22, f, f1, f2;

   x = log10(energy);
   y = log10(redshift);
   x1 = log10(effeff->energy_axis[energy_pt1]);
   x2 = log10(effeff->energy_axis[energy_pt2]);
   y1 = log10(effeff->redshift_axis[redshift_pt1]);
   y2 = log10(effeff->redshift_axis[redshift_pt2]);
   f11 = effeff->efficiency[redshift_pt1][energy_pt1];
   f12 = effeff->efficiency[redshift_pt1][energy_pt2];
   f21 = effeff->efficiency[redshift_pt2][energy_pt1];
   f22 = effeff->efficiency[redshift_pt2][energy_pt2];

   f1 = (f21 - f11) * (y - y1) / (y2 - y1) + f11;
   f2 = (f22 - f12) * (y - y1) / (y2 - y1) + f12;
   f = (f2 - f1) * (x - x1) / (x2 - x1) + f1;

   if (energy < MinE || MaxE < energy || redshift < MinZp || MaxZp < redshift)
   {
      return 0;
   }
   else
   {
      return f;
   }
}

/****************************************************************************************
Total volumic rate of energy realase due to DM annihilation, in eV/cm^3/s.
Can replace by any function (just keep the same name).
Added December 2014
Added DM decay February 2016
****************************************************************************************/
double dE_dtdV_func(double z, REC_COSMOPARAMS *param, EFFECTIVE_EFF *annihilation, EFFECTIVE_EFF *decay, EFFECTIVE_EFF *blackhole)
{

   double dE_dtdV_annihilation, dE_dtdV_decay, dE_dt_dV_blackhole, energy;
   // The energy axis is for DM mass
   energy = param->dm_mass;
   if (param->pann > 0)
   {
      double eff_ann;
      // energy = param->dm_mass;
      eff_ann = interp_eff_eff(annihilation, energy, 1. + z, ANNIHILATION_ENERGY_AXIS_SIZE, ANNIHILATION_REDSHIFT_AXIS_SIZE);
      dE_dtdV_annihilation = eff_ann * (square(10537.4 * (1. - param->dm_bh_frac) * (param->omh2 - param->obh2) * cube(1. + z)) * param->pann);
      /* the prefactor of 10537.4 is 3 H0^2/(8 Pi G) c^2 in eV/cm^3, H0 = 100km/s/Mpc */
   }
   else
   {
      dE_dtdV_annihilation = 0;
   }

   if (param->dm_gamma > 0)
   {
      double eff_decay;
      // energy = param->dm_mass/2;
      eff_decay = interp_eff_eff(decay, energy, 1. + z, DECAY_ENERGY_AXIS_SIZE, DECAY_REDSHIFT_AXIS_SIZE);
      dE_dtdV_decay = eff_decay * (10537.4 * (1. - param->dm_bh_frac) * (param->omh2 - param->obh2) * cube(1. + z) * param->dm_gamma);
      /* the prefactor of 10537.4 is 3 H0^2/(8 Pi G) c^2 in eV/cm^3, H0 = 100km/s/Mpc */
   }
   else
   {
      dE_dtdV_decay = 0;
   }

   if (param->dm_bh_frac > 0)
   {
      double eff_blackhole;
      energy = 1.06e22 / param->bh_mass; /* Blackhole Temperature */
      eff_blackhole = interp_eff_eff(blackhole, energy, 1. + z, BH_ENERGY_AXIS_SIZE, BH_REDSHIFT_AXIS_SIZE);
      dE_dt_dV_blackhole = eff_blackhole * (5.626976744186047e+29 / cube(param->bh_mass) * param->dm_bh_frac * (param->omh2 - param->obh2) * cube(1. + z));
   }
   else
   {
      dE_dt_dV_blackhole = 0;
   }

   return (dE_dtdV_annihilation + dE_dtdV_decay + dE_dt_dV_blackhole);
}

/***********************************************************************************************
Free the memory for effective efficiency tables. Written Sept 2016
***********************************************************************************************/
void free_eff_eff(EFFECTIVE_EFF *effeff, unsigned REDSHIFT_AXIS_SIZE)
{
   free(effeff->energy_axis);
   free(effeff->redshift_axis);
   free_2D_array(effeff->efficiency, REDSHIFT_AXIS_SIZE);
}

/*************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                 */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              */
/*                                                                                               */
/*         history.h: functions for recombination history                                        */
/*                                                                                               */
/*************************************************************************************************/

/**** Structure for cosmological parameters.
      Include information on starting and ending redshit and timestep
      Added May 2012: the period of Hydrogen recombination is evaluated right away so
      tables of radiation field can be reduced to only what is needed and memory is saved.
****/

typedef struct
{
      double T0;               /* CMB temperature today in K*/
      double obh2, omh2, okh2; /* cosmological parameters */
      double odeh2, w0, wa;    /* dark energy parameters */
      double Y;                /* primordial helium abundance */
      double Nnueff;           /* effective number of neutrinos */

      double fsR, meR; /* alpha_fs/alpha_fs(today) and me/me(today) (Added April 2012)*/

      /* Secondary parameters, to avoid recalculating every time */
      double nH0; /* density of hydrogen today in m^{-3} */
      double fHe; /* Helium fraction by number */

      long nz;    /* total number of redshift steps */
      long izH0;  /* index when H recombination starts to be considered */
      double zH0; /* Redshift at which H recombination starts (zH0 = z[izH0]) */
      long nzrt;  /* number of redshift steps while radiative transfer is computed */

      /*-------- Dark Matter and PBH params --------*/
      double Mdm;   // DM mass in GeV
      double Pann;  /* DM annihilation parameter in the smooth background and in haloes */
                    /* Units of pann are cm^3/s/GeV */
      double Gamma; // DM decay width in s^-1
      double Mbh;  /* primordial black hole mass, if PBH_Model==1, Mbh is in solar mass, otherwise Mbh is in gram */
      double fbh;  /* Fraction of DM made of primordial black holes */

      double PBH_Model;  /* Chose PBH energy injection model
                         1 - Accretion
                         2 - Hawking Radiation, with hadronisation emissions, detailed in arxiv 2108.13256
                         3 - Hawking Radiation, with NO hadronisation emissions, detailed in arxiv 2108.13256
                         4 - Hawking Radiation accounting only for electron+positron+gamma emissions, applicapable for Mbh = [10^15, 10^17] g
                         */
      double PBH_Spin;   /* Reduced Kerr spin of PBH, see arxiv 2108.13256 for definition
                         Currently only applicable to Hawking radiation for following values: [0 0.25 0.5 0.75 0.999 0.9999]
                         will automatically find closest matching spin if entered value os not among above list
                         */
      double DM_Channel; /* DM decay/annihilation channel
                     1 - Photon
                     2 - Electron
                     3 - Higgs
                     4 - Muon
                     5 - Tau
                     6 - Q
                     7 - Charm
                     8 - Bottom
                     9 - Top
                     10 - W
                     11 - Z
                     12 - Gluon
                     */

} REC_COSMOPARAMS;

void rec_get_cosmoparam(FILE *fin, FILE *fout, REC_COSMOPARAMS *param);
double rec_HubbleConstant(REC_COSMOPARAMS *param, double z);
double rec_Tmss(double xe, double Tr, double H, double fHe, double fsR, double meR, REC_COSMOPARAMS *param);
double rec_dTmdlna(double xe, double Tm, double Tr, double H, double fHe, double fsR, double meR, REC_COSMOPARAMS *param);
void rec_get_xe_next1_He(REC_COSMOPARAMS *param, double z_in, double *xHeII,
                         double *dxHeIIdlna_prev, double *dxHeIIdlna_prev2, int *post_saha);
double rec_xH1s_postSaha(REC_COSMOPARAMS *param, unsigned iz_out, double z_out, double xHeII_out,
                         HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog_params,
                         double **Dfminus_hist, double *Dfminus_Ly_hist[], double **Dfnu_hist, int *post_saha);
void get_rec_next2_HHe(REC_COSMOPARAMS *param, unsigned iz_in, double z_in, double Tm_in, double *xH1s, double *xHeII,
                       HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog_params, double **Dfminus_hist, double *Dfminus_Ly_hist[], double **Dfnu_hist,
                       double *dxHIIdlna_prev, double *dxHeIIdlna_prev, double *dxHIIdlna_prev2, double *dxHeIIdlna_prev2, int *post_saha);
void rec_get_xe_next1_H(REC_COSMOPARAMS *param, double z_in, double xe_in, double Tm_in, double *xe_out,
                        HRATEEFF *rate_table, unsigned iz, TWO_PHOTON_PARAMS *twog_params, double **Dfminus_hist, double *Dfminus_Ly_hist[],
                        double **Dfnu_hist, double *dxedlna_prev, double *dxedlna_prev2, int *post_saha);
void rec_get_xe_next2_HTm(int func_select, REC_COSMOPARAMS *param, double z_in, double xe_in, double Tm_in, double *xe_out, double *Tm_out,
                          HRATEEFF *rate_table, unsigned iz, TWO_PHOTON_PARAMS *twog_params, double **Dfminus_hist, double *Dfminus_Ly_hist[],
                          double **Dfnu_hist, double *dxedlna_prev, double *dTmdlna_prev, double *dxedlna_prev2, double *dTmdlna_prev2);
void rec_build_history(REC_COSMOPARAMS *param, HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog_params,
                       double *xe_output, double *Tm_output, double **Dfnu_hist, double *Dfminus_Ly_hist[3]);

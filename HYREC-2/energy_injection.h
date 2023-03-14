/* Structure with all energy injection parameters */
/* If adding a new energy injection process
   make sure to add relevant parameters here */

#ifndef __ENERGY_INJECTION__
#define __ENERGY_INJECTION__

typedef struct
{

  double odmh2; /* Omega_dm h^2 */
  
  double Mdm; // DM mass in GeV
  double pann;  /* DM annihilation parameter in the smooth background and in haloes */
                /* Units of pann are cm^3/s/GeV */
  double decay; // DM decay width in s^-1
  double Mpbh; /* primordial black hole mass, if PBH_Model==1, Mbh is in solar mass, otherwise Mbh is in gram */
  double fpbh; /* Fraction of DM made of primordial black holes */
 
  double PBH_Model; /* Chose PBH energy injection model
                    1 - Accretion
                    2 - Hawking Radiation, with hadronisation emissions, detailed in arxiv 2108.13256
                    3 - Hawking Radiation, with NO hadronisation emissions, detailed in arxiv 2108.13256
                    4 - Hawking Radiation accounting only for electron+positron+gamma emissions, applicapable for Mbh = [10^15, 10^17] g
                    */
  double PBH_Spin; /* Reduced Kerr spin of PBH, see arxiv 2108.13256 for definition
                   Currently only applicable to Hawking radiation for following values: [0 0.25 0.5 0.75 0.999 0.9999]
                   will automatically find closest matching spin if entered value os not among above list
                   */
  double ion; // dEdVdt_HIon / nH / EI
  double exclya; // dEdVdt_LyA / nH / E21

  double DM_Channel;/* DM decay/annihilation channel
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

} INJ_PARAMS;

double dEdVdt_deposited(double z, INJ_PARAMS *params, int dep_channel);
void Check_Error(double xe, double T);
void Validate_Inputs(INJ_PARAMS *params);

#endif

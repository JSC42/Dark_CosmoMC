#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Tables.h"

int Sign(double x)
{
  // Find sign of a number
  int result;
  if (x > 0.0)
  {
    result = 1;
  }
  else if (x < 0.0)
  {
    result = -1;
  }
  else
  {
    // maybe x is too close to 0?
    result = 0;
  }
  return result;
}

int Same_Sign(double x, double y)
{
  int sx, sy, r;
  sx = Sign(x);
  sy = Sign(y);
  if (sx == sy)
  {
    r = 1;
  }
  else
  {
    r = 0;
  }
  return r;
}

int Find_Index_1D(double *Tab, double x, int nx)
{
  // Find closest left element index
  // return -1 if not in range
  double x1, x2, x3;
  int id1, id2, id3, Stop, s1, s2, s3, idx, count;
  Stop = 0;
  id1 = 0;
  id3 = nx - 1;

  count = 0;
  while (Stop == 0)
  {
    count = count + 1;
    id2 = (int)round((((double)id1) + ((double)id3)) / 2.0);
    if (id3 == id1 + 1)
    {
      idx = id1;
      // printf("Stopping, id1 = %d, id2 = %d, id3 = %d\n", id1, id2, id3);
      Stop = 1;
    }

    x1 = Tab[id1];
    x2 = Tab[id2];
    x3 = Tab[id3];

    s1 = Sign(x - x1);
    s2 = Sign(x - x2);
    s3 = Sign(x - x3);
    if (s1 != s2)
    {
      id3 = id2;
    }
    else if (s2 != s3)
    {
      id1 = id2;
    }
    if ((s1 == s3) || (count > 50))
    {
      idx = -1;
      Stop = 1;
    }
  }
  return idx;
}

double Read_Tab_2D(int yid, int xid, int ny, int nx, double *Tab)
{
  // A way to read 2D table, I can't get pointer to work
  int id;
  if ((xid > nx - 1) || (yid > ny - 1))
  {
    printf("Error: array index overflow.\n");
    exit(1);
  }
  id = yid * nx + xid;
  return Tab[id];
}

double Interp_2D(double m, double z, double *m_axis, double *z_axis, double *Tab, int nm, int nz, int Overflow_Handle)
{
  /* Interpolate (in log) from a 2D table
   -- inputs --
   m : m target
   z : z target
   m_axis : m axis pointer
   z_axis : z axis pointer
   Tab : Data Table, index : Tab[z_id][m_id]
   nm : m axis size
   nz : z axis size
   Overflow_Handle : decide what to do when m or z is not in range
                     0 : Raise error and exit
                     1 : return 0 when z is not in range, if m is also not in range then raise error
  */
  int mid1, mid2, zid1, zid2;
  double lm, lm1, lm2, lz, lz1, lz2, f11, f12, f21, f22, f1, f2, F1, F2, f;
  mid1 = Find_Index_1D(m_axis, m, nm);
  zid1 = Find_Index_1D(z_axis, z, nz);
  // printf("debug: m = %E,  z = %E\n", m, z);
  
  if ((mid1 == -1) || (zid1 == -1))
  {
    // Overflow detected
    if (Overflow_Handle == 0)
    {
      printf("Error: m or z not in range, exitting, debugging info:.\n");
      printf("m1 = %E, m = %E, m2 = %E\n", m_axis[mid1], m, m_axis[mid1 + 1]);
      printf("z1 = %E, z = %E, z2 = %E\n", z_axis[zid1], z, z_axis[zid1 + 1]);
      printf("min(m) = %E, m = %E, max(m) = %E\n", m_axis[0], m, m_axis[nm - 1]);
      printf("min(z) = %E, z = %E, max(z) = %E\n", z_axis[0], z, z_axis[nz - 1]);
      exit(1);
    }
    else if ((zid1 == -1) && (Overflow_Handle == 1))
    {
      if (mid1 == -1)
      {
        printf("Error: m not in range, exitting, debugging info:.\n");
        printf("m1 = %E, m = %E, m2 = %E\n", m_axis[mid1], m, m_axis[mid1 + 1]);
        printf("z1 = %E, z = %E, z2 = %E\n", z_axis[zid1], z, z_axis[zid1 + 1]);
        printf("min(m) = %E, m = %E, max(m) = %E\n", m_axis[0], m, m_axis[nm - 1]);
        printf("min(z) = %E, z = %E, max(z) = %E\n", z_axis[0], z, z_axis[nz - 1]);
        exit(1);
      }
      else
      {
        // Our z axis is bounded
        return 0.0;
      }
    }
  }

  mid2 = mid1 + 1;
  zid2 = zid1 + 1;

  lm = log10(m);
  lm1 = log10(m_axis[mid1]);
  lm2 = log10(m_axis[mid2]);

  lz = log10(z);
  lz1 = log10(z_axis[zid1]);
  lz2 = log10(z_axis[zid2]);

  f11 = Read_Tab_2D(zid1, mid1, nz, nm, Tab);
  f12 = Read_Tab_2D(zid1, mid2, nz, nm, Tab);
  f21 = Read_Tab_2D(zid2, mid1, nz, nm, Tab);
  f22 = Read_Tab_2D(zid2, mid2, nz, nm, Tab);

  // fix m1
  f1 = f11;
  f2 = f21;
  F1 = (f2 - f1) / (lz2 - lz1) * (lz - lz1) + f1;
  // fix m2
  f1 = f12;
  f2 = f22;
  F2 = (f2 - f1) / (lz2 - lz1) * (lz - lz1) + f1;

  f = (F2 - F1) / (lm2 - lm1) * (lm - lm1) + F1;

  return f;
}

double Interp_EFF_DM_Annihilation(double Mdm, double z, int Dep_Channel, int Reaction_Channel)
{
  /* Interpolate deposition efficiencies for DM annihilation
  -- inputs --
  Mdm : DM mass in GeV
  z: redshfit
  Dep_Channel: Deposition Channel
          1 - HIon
          3 - LyA
          4 - Heat
  Reaction_Channel: Product that DM decay/annihilate into
          1 - Photon
          2 - Electron
          3 - Higgs
          4 - Muon
          5 - Tau
          6 - Q
          7 - CHarm
          8 - Bottom
          9 - Top
          10 - W
          11 - Z
          12 - Gluon
  */

  double Mdm_eV, x_in, EFF;
  Mdm_eV = Mdm * 1.0E9;
  // Validate inputs
  if ((Dep_Channel != 1) && (Dep_Channel != 3) && (Dep_Channel != 4))
  {
    printf("Error: wrong choice of deposition channels.\n");
    exit(1);
  }

  // ---- Gamma ----
  if (Reaction_Channel == 1)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Phot_HIon_Ann_HMG, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Phot_LyA_Ann_HMG, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Phot_Heat_Ann_HMG, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Elec ----
  else if (Reaction_Channel == 2)
  {
    x_in = Mdm_eV - Electron_Mass_eV;
    if (x_in < 0.0)
    {
      printf("Error: DM mass too small to annihilate into electron+positron.\n");
      exit(1);
    }
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Elec_HIon_Ann_HMG, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Elec_LyA_Ann_HMG, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Elec_Heat_Ann_HMG, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Higgs ----
  else if (Reaction_Channel == 3)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Higgs_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Higgs_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Higgs_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Muon ----
  else if (Reaction_Channel == 4)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Muon_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Muon_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Muon_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Tau ----
  else if (Reaction_Channel == 5)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Tau_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Tau_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Tau_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Q ----
  else if (Reaction_Channel == 6)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Q_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Q_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Q_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Charm ----
  else if (Reaction_Channel == 7)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Charm_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Charm_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Charm_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Bottom ----
  else if (Reaction_Channel == 8)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Bottom_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Bottom_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Bottom_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Top ----
  else if (Reaction_Channel == 9)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Top_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Top_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Top_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- W ----
  else if (Reaction_Channel == 10)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_W_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_W_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_W_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Z ----
  else if (Reaction_Channel == 11)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Z_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Z_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Z_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Gluon ----
  else if (Reaction_Channel == 12)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Gluon_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Gluon_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Gluon_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  else
  {
    printf("Unknown DM decay/annihilation channel.\n");
    exit(1);
  }
  return EFF;
}

double Interp_EFF_DM_Decay(double Mdm, double z, int Dep_Channel, int Reaction_Channel)
{
  /* Interpolate deposition efficiencies for DM decay
  -- inputs --
  Mdm : DM mass in GeV
  z: redshfit
  Dep_Channel: Deposition Channel
          1 - HIon
          3 - LyA
          4 - Heat
  Reaction_Channel: Product that DM decay/annihilate into
          1 - Photon
          2 - Electron
          3 - Higgs
          4 - Muon
          5 - Tau
          6 - Q
          7 - CHarm
          8 - Bottom
          9 - Top
          10 - W
          11 - Z
          12 - Gluon
  */

  double Mdm_eV, x_in, EFF;
  Mdm_eV = Mdm * 1.0E9;
  // Validate inputs
  if ((Dep_Channel != 1) && (Dep_Channel != 3) && (Dep_Channel != 4))
  {
    printf("Error: wrong choice of deposition channels\n");
    exit(1);
  }

  // ---- Gamma ----
  if (Reaction_Channel == 1)
  {
    x_in = Mdm_eV / 2.0;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Phot_HIon_Decay, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Phot_LyA_Decay, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Phot_Heat_Decay, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Electron ----
  else if (Reaction_Channel == 2)
  {
    x_in = Mdm_eV / 2.0 - Electron_Mass_eV;
    if (x_in < 0.0)
    {
      printf("Error: DM mass too small to decay into electron+positron.\n");
      exit(1);
    }

    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Elec_HIon_Decay, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Elec_LyA_Decay, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Elec_Heat_Decay, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Higgs ----
  else if (Reaction_Channel == 3)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Higgs_HIon_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Higgs_LyA_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Higgs_Heat_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Muon ----
  else if (Reaction_Channel == 4)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Muon_HIon_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Muon_LyA_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Muon_Heat_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Tau ----
  else if (Reaction_Channel == 5)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Tau_HIon_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Tau_LyA_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Tau_Heat_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Q ----
  else if (Reaction_Channel == 6)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Q_HIon_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Q_LyA_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Q_Heat_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Charm ----
  else if (Reaction_Channel == 7)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Charm_HIon_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Charm_LyA_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Charm_Heat_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Bottom ----
  else if (Reaction_Channel == 8)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Bottom_HIon_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Bottom_LyA_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Bottom_Heat_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Top ----
  else if (Reaction_Channel == 9)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Top_HIon_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Top_LyA_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Top_Heat_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- W ----
  else if (Reaction_Channel == 10)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_W_HIon_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_W_LyA_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_W_Heat_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Z ----
  else if (Reaction_Channel == 11)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Z_HIon_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Z_LyA_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Z_Heat_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Gluon ----
  else if (Reaction_Channel == 12)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Gluon_HIon_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Gluon_LyA_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Decay, Redshift_Axis, EFF_Gluon_Heat_Decay, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  else
  {
    printf("Unknown DM decay/annihilation channel.\n");
    exit(1);
  }

  return EFF;
}

int Convert_to_Int(double x)
{
  int ix;
  ix = (int)round(x);
  return ix;
}

double Interp_EFF_Hawking(double Mbh, double z, double PBH_Spin, double PBH_Model, int Dep_Channel)
{
  /* Interpolate PBH Hawking Radiation deposition efficiencies
  inputs:
  Mbh : PBH mass in gram
  z : redshift
  PBH_Spin : PBH_Spin
  PBH_Model : PBH_Model defined in energy_injection.h
  Dep_Channel: Deposition Channel
          1 - HIon
          3 - LyA
          4 - Heat
  */
  double EFF, Tbh;
  int Spin_id, model_id;
  model_id = Convert_to_Int(PBH_Model);
  Tbh = 1.06E6 * (1.0E16 / Mbh); // BH temperature in eV
  if (PBH_Spin < 0.2)
  {
    Spin_id = 1;
  }
  else if ((PBH_Spin > 0.24) && (PBH_Spin < 0.26))
  {
    Spin_id = 2;
  }
  else if ((PBH_Spin > 0.49) && (PBH_Spin < 0.51))
  {
    Spin_id = 3;
  }
  else if ((PBH_Spin > 0.74) && (PBH_Spin < 0.76))
  {
    Spin_id = 4;
  }
  else if ((PBH_Spin > 0.991) && (PBH_Spin < 0.9995))
  {
    Spin_id = 5;
  }
  else if ((PBH_Spin > 0.9996) && (PBH_Spin < 0.99995))
  {
    Spin_id = 6;
  }

  if (Spin_id == 1)
  {
    if (model_id == 2)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K1, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K1, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K1, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
    else if (model_id == 3)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K1B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K1B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K1B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
  }
  else if (Spin_id == 2)
  {
    if (model_id == 2)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K2, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K2, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K2, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
    else if (model_id == 3)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K2B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K2B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K2B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
  }
  else if (Spin_id == 3)
  {
    if (model_id == 2)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K3, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K3, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K3, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
    else if (model_id == 3)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K3B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K3B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K3B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
  }
  else if (Spin_id == 4)
  {
    if (model_id == 2)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K4, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K4, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K4, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
    else if (model_id == 3)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K4B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K4B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K4B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
  }
  else if (Spin_id == 5)
  {
    if (model_id == 2)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K5, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K5, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K5, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
    else if (model_id == 3)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K5B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K5B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K5B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
  }
  else if (Spin_id == 6)
  {
    if (model_id == 2)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K6, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K6, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K6, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
    else if (model_id == 3)
    {
      if (Dep_Channel == 1)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_HIon_K6B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 3)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_LyA_K6B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
      else if (Dep_Channel == 4)
      {
        EFF = Interp_2D(Tbh, 1 + z, Kerr_PBH_Temperature_Axis, Redshift_Axis, EFF_Hawking_Heat_K6B, Kerr_PBH_Temperature_Size, Redshift_Size, 1);
      }
    }
  }

  return EFF;
}

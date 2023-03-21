#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "EFF_Tables.h"
#include "Useful_Functions.h"

/* ---- main ---- */

/* Check mdot_nacked with Fig.3 of 0709.0524 */
double Tm_Fit(double z)
{
  double adec, zdec, zp, b, r, a;

  zdec = 132;
  b = 1.72;
  adec = 1 / (1 + zdec);
  zp = 1 + z;
  a = 1 / zp;
  r = 2730 * (zp / 1000) * adec / pow((pow(a, b) + pow(adec, b)), 1 / b);

  return r;
}

double Veff(double z)
{
  double Cs, Vr, Tm, r;

  Cs = 5700 * sqrt(Tm_Fit(z) / 2730);
  Vr = 30000 * fmin(1, (1 + z) / 1000);
  r = sqrt(pow(Cs, 2) + pow(Vr, 2));
  return r;
}

double viscosity(double m, double z, double xe)
{
  // Gas viscosity withot DM accretion, Eq.3 of 2003.12589
  // m : bh mass in msun
  double r, zp, veff, T, Cs;
  zp = 1 + z;
  T = Tm_Fit(z);

  r = 0.257 + 1.45 * (xe / 0.01) * pow(zp / 1000, 2.5);
  r *= (m / 1E4) * pow(zp / 1000, 1.5);
  veff = Veff(z);
  r *= pow(veff / 5740, -3);

  return r;
}

double Lambda(double m, double z, double xe, double b)
{
  // BH accretion efficiency
  double r, xcr, T;
  T = Tm_Fit(z);
  xcr = (sqrt(1 + b) - 1) / b;
  r = xcr * xcr * exp(4.5 / (3 + pow(b, 0.75)));
  return r;
}

double mdot_naked(double m, double z, double xe)
{
  /* Dimensionless accretion rate for a nacked PBH
   */
  double mdot, r, l, vr, zp, T, b;
  T = Tm_Fit(z);
  vr = Veff(z);
  b = viscosity(m, z, xe);
  l = Lambda(m, z, xe, b);
  zp = 1 + z;
  mdot = 0.4 * l * pow(zp / 1000, 3) * m * pow(vr / 1000, -3);

  return mdot;
}

double K_factor(double m, double z, double xe)
{
  // BH K factor, see 2003.12589
  double mh, zp, veff, r;
  zp = 1 + z;
  mh = 3 * m * 1000 / zp;
  veff = Veff(z);
  r = 0.22 * (zp / 1000) * pow(mh, 2 / 3) * pow(veff / 1000, -2);

  return r;
}


double mdot_clothed(double m, double z, double xe)
{
  /* Accretion rate for PBH clothed in DM particle
   */
  double k, r, zp, a, bh, b0, lh, y, veff, mdot, mh, p, T, l0;
  a = 2.25; // halo density profile, 2108.11130
  T = Tm_Fit(z);
  zp = 1 + z;
  k = K_factor(m, z, xe);
  p = 2 - a;
  mh = 3 * m * 1000 / zp;
  if (k > 2)
  {
    printf("hello\n");
    r = mdot_naked(mh, z, xe);
    
    return r;
  }

  b0 = viscosity(m, z, xe);
  bh = b0 * pow(k, p / (1 - p));
  y = pow(1 + 10 * bh, 0.1);
  y *= exp(2 - k) * pow(k / 2, 2);
  l0 = Lambda(m, z, xe, bh);
  lh = l0*pow(y, p / (1 - p));

  veff = Veff(z);
  
  // dimensionless
  // mdot = 0.016 * lh * (zp / 1000) * m * pow(veff / 5740, -3);
  mdot = 0.016 * l0 * (zp / 1000) * m * pow(veff / 5740, -3);

  return mdot;
}

int main()
{
  double vz[100];
  int nz = 100;
  int id;
  double m, xe, z, mdot;

  xe = 1;
  m = 1E5;

  Fill_Linear_Array(1, 3000, vz, nz, 1);

  for (id = 0; id < nz; id++)
  {
    z = vz[id];
    mdot = mdot_clothed(m, z, xe);
    printf("%f  %E\n", 1 + z, mdot);
  }

  return 0;
}

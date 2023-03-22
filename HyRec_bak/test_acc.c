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
    mdot = mdot_naked(m, z, xe);
    printf("%f  %E\n", 1 + z, mdot);
  }

  return 0;
}

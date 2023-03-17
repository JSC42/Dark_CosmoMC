#include "Functions.h"
int main()
{
  double Ek, zp, r;
  scanf("%lf\n", &Ek);
  scanf("%lf", &zp);
  r = Interp_2D_v2(Ek, zp, Kinetic_Energy_Axis, Redshift_Axis, EFF_Phot_LyA_Ann_HMG, Kinetic_Energy_Axis_Size, Redshift_Size);
  printf("%E\n",r);
  
  return 0;
  
}
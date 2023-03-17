#include "DarkSide.h"
// Functions.h renamed into DarkSide.h

int main()
{
  double Mdm, z;
  double r;
  scanf("%lf", &Mdm);
  scanf("%lf", &z);
  r = Interp_EFF_DM_Annihilation(Mdm, z, 1, 1);

  printf("%f\n", r);
  return 0;
  
}
#include "Functions.h"

int main()
{
  int mid, zid;
  double r;
  scanf("%d", &mid);
  scanf("%d", &zid);
  r = Read_Tab_2D(zid, mid, Redshift_Size, Kinetic_Energy_Axis_Size, EFF_Phot_HIon_Ann_HMG);
  //double Read_Tab_2D(int yid, int xid, int ny, int nx, double *Tab)
  printf("%f\n", r);
  return 0;
  
}
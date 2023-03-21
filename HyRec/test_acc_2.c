#include "Accreting_PBH.h"
int main()
{
  double m0, z_out[Mbh_Evo_nz], m_out[Mbh_Evo_nz];
  // m0 = 30;
  scanf("%lf",&m0);
  Get_Mass_Evolution(m0, z_out, m_out);
  return 0;
}
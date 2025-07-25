#include "globes/globes.h"

#define MAX_INTERACTIONS 5


#include <complex.h>

struct WEFT {
  double complex epsilon_CC[MAX_INTERACTIONS][2][3][MAX_FLAVORS][MAX_FLAVORS]; // Interaction index("L", "R", "S", "P", "T"),up-like index, down-like index,  charged lepton, neutrino
  double complex epsilon_NC[2][3][3][MAX_FLAVORS][MAX_FLAVORS]; // Interaction index("L", "R"),matterfermions, matterfermions,  laptonflavors, laptonflavors
};

struct SMEFT {
  double complex c3_phi_l[3][3],c_phi_u[3][3],c_phi_d[3][3],c_phi_e[3][3],c_phi_q[3][3],c_phi_l[3][3], c3_phi_q[3][3],c3_phi_e[3][3],c3_phi_le[3][3],c_phi_ud[3][3],c_ll[3][3][3][3], c3_lq[3][3][3][3], c_lq[3][3][3][3],c_ld[3][3][3][3],c_lu[3][3][3][3],c_le[3][3][3][3],c_ledq[3][3][3][3],c_lequ[3][3][3][3],c3_lequ[3][3][3][3],c3_ledq[3][3][3][3], c_HWB, c_HD;

};

void Running_and_matching(double running_smeft[3][3], double running_weft_mb[5][5],
                          double running_weft_low_scale[5][5], struct SMEFT High_scale,
                          struct WEFT *Low_scale);

void running_filling(double running_smeft[3][3], double running_weft_mb[5][5],
                     double running_weft_low_scale[5][5]);

double f(double T3, double Q, struct SMEFT smeft_param);




#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with sympy 1.9                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_4853822566850710957) {
   out_4853822566850710957[0] = delta_x[0] + nom_x[0];
   out_4853822566850710957[1] = delta_x[1] + nom_x[1];
   out_4853822566850710957[2] = delta_x[2] + nom_x[2];
   out_4853822566850710957[3] = delta_x[3] + nom_x[3];
   out_4853822566850710957[4] = delta_x[4] + nom_x[4];
   out_4853822566850710957[5] = delta_x[5] + nom_x[5];
   out_4853822566850710957[6] = delta_x[6] + nom_x[6];
   out_4853822566850710957[7] = delta_x[7] + nom_x[7];
   out_4853822566850710957[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_7910005515458276405) {
   out_7910005515458276405[0] = -nom_x[0] + true_x[0];
   out_7910005515458276405[1] = -nom_x[1] + true_x[1];
   out_7910005515458276405[2] = -nom_x[2] + true_x[2];
   out_7910005515458276405[3] = -nom_x[3] + true_x[3];
   out_7910005515458276405[4] = -nom_x[4] + true_x[4];
   out_7910005515458276405[5] = -nom_x[5] + true_x[5];
   out_7910005515458276405[6] = -nom_x[6] + true_x[6];
   out_7910005515458276405[7] = -nom_x[7] + true_x[7];
   out_7910005515458276405[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_4981524804303897279) {
   out_4981524804303897279[0] = 1.0;
   out_4981524804303897279[1] = 0;
   out_4981524804303897279[2] = 0;
   out_4981524804303897279[3] = 0;
   out_4981524804303897279[4] = 0;
   out_4981524804303897279[5] = 0;
   out_4981524804303897279[6] = 0;
   out_4981524804303897279[7] = 0;
   out_4981524804303897279[8] = 0;
   out_4981524804303897279[9] = 0;
   out_4981524804303897279[10] = 1.0;
   out_4981524804303897279[11] = 0;
   out_4981524804303897279[12] = 0;
   out_4981524804303897279[13] = 0;
   out_4981524804303897279[14] = 0;
   out_4981524804303897279[15] = 0;
   out_4981524804303897279[16] = 0;
   out_4981524804303897279[17] = 0;
   out_4981524804303897279[18] = 0;
   out_4981524804303897279[19] = 0;
   out_4981524804303897279[20] = 1.0;
   out_4981524804303897279[21] = 0;
   out_4981524804303897279[22] = 0;
   out_4981524804303897279[23] = 0;
   out_4981524804303897279[24] = 0;
   out_4981524804303897279[25] = 0;
   out_4981524804303897279[26] = 0;
   out_4981524804303897279[27] = 0;
   out_4981524804303897279[28] = 0;
   out_4981524804303897279[29] = 0;
   out_4981524804303897279[30] = 1.0;
   out_4981524804303897279[31] = 0;
   out_4981524804303897279[32] = 0;
   out_4981524804303897279[33] = 0;
   out_4981524804303897279[34] = 0;
   out_4981524804303897279[35] = 0;
   out_4981524804303897279[36] = 0;
   out_4981524804303897279[37] = 0;
   out_4981524804303897279[38] = 0;
   out_4981524804303897279[39] = 0;
   out_4981524804303897279[40] = 1.0;
   out_4981524804303897279[41] = 0;
   out_4981524804303897279[42] = 0;
   out_4981524804303897279[43] = 0;
   out_4981524804303897279[44] = 0;
   out_4981524804303897279[45] = 0;
   out_4981524804303897279[46] = 0;
   out_4981524804303897279[47] = 0;
   out_4981524804303897279[48] = 0;
   out_4981524804303897279[49] = 0;
   out_4981524804303897279[50] = 1.0;
   out_4981524804303897279[51] = 0;
   out_4981524804303897279[52] = 0;
   out_4981524804303897279[53] = 0;
   out_4981524804303897279[54] = 0;
   out_4981524804303897279[55] = 0;
   out_4981524804303897279[56] = 0;
   out_4981524804303897279[57] = 0;
   out_4981524804303897279[58] = 0;
   out_4981524804303897279[59] = 0;
   out_4981524804303897279[60] = 1.0;
   out_4981524804303897279[61] = 0;
   out_4981524804303897279[62] = 0;
   out_4981524804303897279[63] = 0;
   out_4981524804303897279[64] = 0;
   out_4981524804303897279[65] = 0;
   out_4981524804303897279[66] = 0;
   out_4981524804303897279[67] = 0;
   out_4981524804303897279[68] = 0;
   out_4981524804303897279[69] = 0;
   out_4981524804303897279[70] = 1.0;
   out_4981524804303897279[71] = 0;
   out_4981524804303897279[72] = 0;
   out_4981524804303897279[73] = 0;
   out_4981524804303897279[74] = 0;
   out_4981524804303897279[75] = 0;
   out_4981524804303897279[76] = 0;
   out_4981524804303897279[77] = 0;
   out_4981524804303897279[78] = 0;
   out_4981524804303897279[79] = 0;
   out_4981524804303897279[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_3473102566445950304) {
   out_3473102566445950304[0] = state[0];
   out_3473102566445950304[1] = state[1];
   out_3473102566445950304[2] = state[2];
   out_3473102566445950304[3] = state[3];
   out_3473102566445950304[4] = state[4];
   out_3473102566445950304[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_3473102566445950304[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_3473102566445950304[7] = state[7];
   out_3473102566445950304[8] = state[8];
}
void F_fun(double *state, double dt, double *out_1612317163615163753) {
   out_1612317163615163753[0] = 1;
   out_1612317163615163753[1] = 0;
   out_1612317163615163753[2] = 0;
   out_1612317163615163753[3] = 0;
   out_1612317163615163753[4] = 0;
   out_1612317163615163753[5] = 0;
   out_1612317163615163753[6] = 0;
   out_1612317163615163753[7] = 0;
   out_1612317163615163753[8] = 0;
   out_1612317163615163753[9] = 0;
   out_1612317163615163753[10] = 1;
   out_1612317163615163753[11] = 0;
   out_1612317163615163753[12] = 0;
   out_1612317163615163753[13] = 0;
   out_1612317163615163753[14] = 0;
   out_1612317163615163753[15] = 0;
   out_1612317163615163753[16] = 0;
   out_1612317163615163753[17] = 0;
   out_1612317163615163753[18] = 0;
   out_1612317163615163753[19] = 0;
   out_1612317163615163753[20] = 1;
   out_1612317163615163753[21] = 0;
   out_1612317163615163753[22] = 0;
   out_1612317163615163753[23] = 0;
   out_1612317163615163753[24] = 0;
   out_1612317163615163753[25] = 0;
   out_1612317163615163753[26] = 0;
   out_1612317163615163753[27] = 0;
   out_1612317163615163753[28] = 0;
   out_1612317163615163753[29] = 0;
   out_1612317163615163753[30] = 1;
   out_1612317163615163753[31] = 0;
   out_1612317163615163753[32] = 0;
   out_1612317163615163753[33] = 0;
   out_1612317163615163753[34] = 0;
   out_1612317163615163753[35] = 0;
   out_1612317163615163753[36] = 0;
   out_1612317163615163753[37] = 0;
   out_1612317163615163753[38] = 0;
   out_1612317163615163753[39] = 0;
   out_1612317163615163753[40] = 1;
   out_1612317163615163753[41] = 0;
   out_1612317163615163753[42] = 0;
   out_1612317163615163753[43] = 0;
   out_1612317163615163753[44] = 0;
   out_1612317163615163753[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_1612317163615163753[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_1612317163615163753[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1612317163615163753[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1612317163615163753[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_1612317163615163753[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_1612317163615163753[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_1612317163615163753[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_1612317163615163753[53] = -9.8000000000000007*dt;
   out_1612317163615163753[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_1612317163615163753[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_1612317163615163753[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1612317163615163753[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1612317163615163753[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_1612317163615163753[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_1612317163615163753[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_1612317163615163753[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1612317163615163753[62] = 0;
   out_1612317163615163753[63] = 0;
   out_1612317163615163753[64] = 0;
   out_1612317163615163753[65] = 0;
   out_1612317163615163753[66] = 0;
   out_1612317163615163753[67] = 0;
   out_1612317163615163753[68] = 0;
   out_1612317163615163753[69] = 0;
   out_1612317163615163753[70] = 1;
   out_1612317163615163753[71] = 0;
   out_1612317163615163753[72] = 0;
   out_1612317163615163753[73] = 0;
   out_1612317163615163753[74] = 0;
   out_1612317163615163753[75] = 0;
   out_1612317163615163753[76] = 0;
   out_1612317163615163753[77] = 0;
   out_1612317163615163753[78] = 0;
   out_1612317163615163753[79] = 0;
   out_1612317163615163753[80] = 1;
}
void h_25(double *state, double *unused, double *out_4217284194084537253) {
   out_4217284194084537253[0] = state[6];
}
void H_25(double *state, double *unused, double *out_7353830765589797144) {
   out_7353830765589797144[0] = 0;
   out_7353830765589797144[1] = 0;
   out_7353830765589797144[2] = 0;
   out_7353830765589797144[3] = 0;
   out_7353830765589797144[4] = 0;
   out_7353830765589797144[5] = 0;
   out_7353830765589797144[6] = 1;
   out_7353830765589797144[7] = 0;
   out_7353830765589797144[8] = 0;
}
void h_24(double *state, double *unused, double *out_7320050265028110894) {
   out_7320050265028110894[0] = state[4];
   out_7320050265028110894[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4517341501528236371) {
   out_4517341501528236371[0] = 0;
   out_4517341501528236371[1] = 0;
   out_4517341501528236371[2] = 0;
   out_4517341501528236371[3] = 0;
   out_4517341501528236371[4] = 1;
   out_4517341501528236371[5] = 0;
   out_4517341501528236371[6] = 0;
   out_4517341501528236371[7] = 0;
   out_4517341501528236371[8] = 0;
   out_4517341501528236371[9] = 0;
   out_4517341501528236371[10] = 0;
   out_4517341501528236371[11] = 0;
   out_4517341501528236371[12] = 0;
   out_4517341501528236371[13] = 0;
   out_4517341501528236371[14] = 1;
   out_4517341501528236371[15] = 0;
   out_4517341501528236371[16] = 0;
   out_4517341501528236371[17] = 0;
}
void h_30(double *state, double *unused, double *out_7458624653274663427) {
   out_7458624653274663427[0] = state[4];
}
void H_30(double *state, double *unused, double *out_4176222966628137717) {
   out_4176222966628137717[0] = 0;
   out_4176222966628137717[1] = 0;
   out_4176222966628137717[2] = 0;
   out_4176222966628137717[3] = 0;
   out_4176222966628137717[4] = 1;
   out_4176222966628137717[5] = 0;
   out_4176222966628137717[6] = 0;
   out_4176222966628137717[7] = 0;
   out_4176222966628137717[8] = 0;
}
void h_26(double *state, double *unused, double *out_2456951385613817782) {
   out_2456951385613817782[0] = state[7];
}
void H_26(double *state, double *unused, double *out_3612327446715740920) {
   out_3612327446715740920[0] = 0;
   out_3612327446715740920[1] = 0;
   out_3612327446715740920[2] = 0;
   out_3612327446715740920[3] = 0;
   out_3612327446715740920[4] = 0;
   out_3612327446715740920[5] = 0;
   out_3612327446715740920[6] = 0;
   out_3612327446715740920[7] = 1;
   out_3612327446715740920[8] = 0;
}
void h_27(double *state, double *unused, double *out_5160794046986199688) {
   out_5160794046986199688[0] = state[3];
}
void H_27(double *state, double *unused, double *out_6350986278428562628) {
   out_6350986278428562628[0] = 0;
   out_6350986278428562628[1] = 0;
   out_6350986278428562628[2] = 0;
   out_6350986278428562628[3] = 1;
   out_6350986278428562628[4] = 0;
   out_6350986278428562628[5] = 0;
   out_6350986278428562628[6] = 0;
   out_6350986278428562628[7] = 0;
   out_6350986278428562628[8] = 0;
}
void h_29(double *state, double *unused, double *out_1059583135860814281) {
   out_1059583135860814281[0] = state[1];
}
void H_29(double *state, double *unused, double *out_3665991622313745533) {
   out_3665991622313745533[0] = 0;
   out_3665991622313745533[1] = 1;
   out_3665991622313745533[2] = 0;
   out_3665991622313745533[3] = 0;
   out_3665991622313745533[4] = 0;
   out_3665991622313745533[5] = 0;
   out_3665991622313745533[6] = 0;
   out_3665991622313745533[7] = 0;
   out_3665991622313745533[8] = 0;
}
void h_28(double *state, double *unused, double *out_4557478162309722171) {
   out_4557478162309722171[0] = state[0];
}
void H_28(double *state, double *unused, double *out_8748390639383276107) {
   out_8748390639383276107[0] = 1;
   out_8748390639383276107[1] = 0;
   out_8748390639383276107[2] = 0;
   out_8748390639383276107[3] = 0;
   out_8748390639383276107[4] = 0;
   out_8748390639383276107[5] = 0;
   out_8748390639383276107[6] = 0;
   out_8748390639383276107[7] = 0;
   out_8748390639383276107[8] = 0;
}
void h_31(double *state, double *unused, double *out_3675655300799986106) {
   out_3675655300799986106[0] = state[8];
}
void H_31(double *state, double *unused, double *out_7384476727466757572) {
   out_7384476727466757572[0] = 0;
   out_7384476727466757572[1] = 0;
   out_7384476727466757572[2] = 0;
   out_7384476727466757572[3] = 0;
   out_7384476727466757572[4] = 0;
   out_7384476727466757572[5] = 0;
   out_7384476727466757572[6] = 0;
   out_7384476727466757572[7] = 0;
   out_7384476727466757572[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_4853822566850710957) {
  err_fun(nom_x, delta_x, out_4853822566850710957);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_7910005515458276405) {
  inv_err_fun(nom_x, true_x, out_7910005515458276405);
}
void car_H_mod_fun(double *state, double *out_4981524804303897279) {
  H_mod_fun(state, out_4981524804303897279);
}
void car_f_fun(double *state, double dt, double *out_3473102566445950304) {
  f_fun(state,  dt, out_3473102566445950304);
}
void car_F_fun(double *state, double dt, double *out_1612317163615163753) {
  F_fun(state,  dt, out_1612317163615163753);
}
void car_h_25(double *state, double *unused, double *out_4217284194084537253) {
  h_25(state, unused, out_4217284194084537253);
}
void car_H_25(double *state, double *unused, double *out_7353830765589797144) {
  H_25(state, unused, out_7353830765589797144);
}
void car_h_24(double *state, double *unused, double *out_7320050265028110894) {
  h_24(state, unused, out_7320050265028110894);
}
void car_H_24(double *state, double *unused, double *out_4517341501528236371) {
  H_24(state, unused, out_4517341501528236371);
}
void car_h_30(double *state, double *unused, double *out_7458624653274663427) {
  h_30(state, unused, out_7458624653274663427);
}
void car_H_30(double *state, double *unused, double *out_4176222966628137717) {
  H_30(state, unused, out_4176222966628137717);
}
void car_h_26(double *state, double *unused, double *out_2456951385613817782) {
  h_26(state, unused, out_2456951385613817782);
}
void car_H_26(double *state, double *unused, double *out_3612327446715740920) {
  H_26(state, unused, out_3612327446715740920);
}
void car_h_27(double *state, double *unused, double *out_5160794046986199688) {
  h_27(state, unused, out_5160794046986199688);
}
void car_H_27(double *state, double *unused, double *out_6350986278428562628) {
  H_27(state, unused, out_6350986278428562628);
}
void car_h_29(double *state, double *unused, double *out_1059583135860814281) {
  h_29(state, unused, out_1059583135860814281);
}
void car_H_29(double *state, double *unused, double *out_3665991622313745533) {
  H_29(state, unused, out_3665991622313745533);
}
void car_h_28(double *state, double *unused, double *out_4557478162309722171) {
  h_28(state, unused, out_4557478162309722171);
}
void car_H_28(double *state, double *unused, double *out_8748390639383276107) {
  H_28(state, unused, out_8748390639383276107);
}
void car_h_31(double *state, double *unused, double *out_3675655300799986106) {
  h_31(state, unused, out_3675655300799986106);
}
void car_H_31(double *state, double *unused, double *out_7384476727466757572) {
  H_31(state, unused, out_7384476727466757572);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);

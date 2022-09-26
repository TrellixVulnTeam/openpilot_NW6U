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
void err_fun(double *nom_x, double *delta_x, double *out_3495569113160768895) {
   out_3495569113160768895[0] = delta_x[0] + nom_x[0];
   out_3495569113160768895[1] = delta_x[1] + nom_x[1];
   out_3495569113160768895[2] = delta_x[2] + nom_x[2];
   out_3495569113160768895[3] = delta_x[3] + nom_x[3];
   out_3495569113160768895[4] = delta_x[4] + nom_x[4];
   out_3495569113160768895[5] = delta_x[5] + nom_x[5];
   out_3495569113160768895[6] = delta_x[6] + nom_x[6];
   out_3495569113160768895[7] = delta_x[7] + nom_x[7];
   out_3495569113160768895[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_3332561543766221488) {
   out_3332561543766221488[0] = -nom_x[0] + true_x[0];
   out_3332561543766221488[1] = -nom_x[1] + true_x[1];
   out_3332561543766221488[2] = -nom_x[2] + true_x[2];
   out_3332561543766221488[3] = -nom_x[3] + true_x[3];
   out_3332561543766221488[4] = -nom_x[4] + true_x[4];
   out_3332561543766221488[5] = -nom_x[5] + true_x[5];
   out_3332561543766221488[6] = -nom_x[6] + true_x[6];
   out_3332561543766221488[7] = -nom_x[7] + true_x[7];
   out_3332561543766221488[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_1059310438201751306) {
   out_1059310438201751306[0] = 1.0;
   out_1059310438201751306[1] = 0;
   out_1059310438201751306[2] = 0;
   out_1059310438201751306[3] = 0;
   out_1059310438201751306[4] = 0;
   out_1059310438201751306[5] = 0;
   out_1059310438201751306[6] = 0;
   out_1059310438201751306[7] = 0;
   out_1059310438201751306[8] = 0;
   out_1059310438201751306[9] = 0;
   out_1059310438201751306[10] = 1.0;
   out_1059310438201751306[11] = 0;
   out_1059310438201751306[12] = 0;
   out_1059310438201751306[13] = 0;
   out_1059310438201751306[14] = 0;
   out_1059310438201751306[15] = 0;
   out_1059310438201751306[16] = 0;
   out_1059310438201751306[17] = 0;
   out_1059310438201751306[18] = 0;
   out_1059310438201751306[19] = 0;
   out_1059310438201751306[20] = 1.0;
   out_1059310438201751306[21] = 0;
   out_1059310438201751306[22] = 0;
   out_1059310438201751306[23] = 0;
   out_1059310438201751306[24] = 0;
   out_1059310438201751306[25] = 0;
   out_1059310438201751306[26] = 0;
   out_1059310438201751306[27] = 0;
   out_1059310438201751306[28] = 0;
   out_1059310438201751306[29] = 0;
   out_1059310438201751306[30] = 1.0;
   out_1059310438201751306[31] = 0;
   out_1059310438201751306[32] = 0;
   out_1059310438201751306[33] = 0;
   out_1059310438201751306[34] = 0;
   out_1059310438201751306[35] = 0;
   out_1059310438201751306[36] = 0;
   out_1059310438201751306[37] = 0;
   out_1059310438201751306[38] = 0;
   out_1059310438201751306[39] = 0;
   out_1059310438201751306[40] = 1.0;
   out_1059310438201751306[41] = 0;
   out_1059310438201751306[42] = 0;
   out_1059310438201751306[43] = 0;
   out_1059310438201751306[44] = 0;
   out_1059310438201751306[45] = 0;
   out_1059310438201751306[46] = 0;
   out_1059310438201751306[47] = 0;
   out_1059310438201751306[48] = 0;
   out_1059310438201751306[49] = 0;
   out_1059310438201751306[50] = 1.0;
   out_1059310438201751306[51] = 0;
   out_1059310438201751306[52] = 0;
   out_1059310438201751306[53] = 0;
   out_1059310438201751306[54] = 0;
   out_1059310438201751306[55] = 0;
   out_1059310438201751306[56] = 0;
   out_1059310438201751306[57] = 0;
   out_1059310438201751306[58] = 0;
   out_1059310438201751306[59] = 0;
   out_1059310438201751306[60] = 1.0;
   out_1059310438201751306[61] = 0;
   out_1059310438201751306[62] = 0;
   out_1059310438201751306[63] = 0;
   out_1059310438201751306[64] = 0;
   out_1059310438201751306[65] = 0;
   out_1059310438201751306[66] = 0;
   out_1059310438201751306[67] = 0;
   out_1059310438201751306[68] = 0;
   out_1059310438201751306[69] = 0;
   out_1059310438201751306[70] = 1.0;
   out_1059310438201751306[71] = 0;
   out_1059310438201751306[72] = 0;
   out_1059310438201751306[73] = 0;
   out_1059310438201751306[74] = 0;
   out_1059310438201751306[75] = 0;
   out_1059310438201751306[76] = 0;
   out_1059310438201751306[77] = 0;
   out_1059310438201751306[78] = 0;
   out_1059310438201751306[79] = 0;
   out_1059310438201751306[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_7662000099711208314) {
   out_7662000099711208314[0] = state[0];
   out_7662000099711208314[1] = state[1];
   out_7662000099711208314[2] = state[2];
   out_7662000099711208314[3] = state[3];
   out_7662000099711208314[4] = state[4];
   out_7662000099711208314[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_7662000099711208314[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_7662000099711208314[7] = state[7];
   out_7662000099711208314[8] = state[8];
}
void F_fun(double *state, double dt, double *out_1291851874522726762) {
   out_1291851874522726762[0] = 1;
   out_1291851874522726762[1] = 0;
   out_1291851874522726762[2] = 0;
   out_1291851874522726762[3] = 0;
   out_1291851874522726762[4] = 0;
   out_1291851874522726762[5] = 0;
   out_1291851874522726762[6] = 0;
   out_1291851874522726762[7] = 0;
   out_1291851874522726762[8] = 0;
   out_1291851874522726762[9] = 0;
   out_1291851874522726762[10] = 1;
   out_1291851874522726762[11] = 0;
   out_1291851874522726762[12] = 0;
   out_1291851874522726762[13] = 0;
   out_1291851874522726762[14] = 0;
   out_1291851874522726762[15] = 0;
   out_1291851874522726762[16] = 0;
   out_1291851874522726762[17] = 0;
   out_1291851874522726762[18] = 0;
   out_1291851874522726762[19] = 0;
   out_1291851874522726762[20] = 1;
   out_1291851874522726762[21] = 0;
   out_1291851874522726762[22] = 0;
   out_1291851874522726762[23] = 0;
   out_1291851874522726762[24] = 0;
   out_1291851874522726762[25] = 0;
   out_1291851874522726762[26] = 0;
   out_1291851874522726762[27] = 0;
   out_1291851874522726762[28] = 0;
   out_1291851874522726762[29] = 0;
   out_1291851874522726762[30] = 1;
   out_1291851874522726762[31] = 0;
   out_1291851874522726762[32] = 0;
   out_1291851874522726762[33] = 0;
   out_1291851874522726762[34] = 0;
   out_1291851874522726762[35] = 0;
   out_1291851874522726762[36] = 0;
   out_1291851874522726762[37] = 0;
   out_1291851874522726762[38] = 0;
   out_1291851874522726762[39] = 0;
   out_1291851874522726762[40] = 1;
   out_1291851874522726762[41] = 0;
   out_1291851874522726762[42] = 0;
   out_1291851874522726762[43] = 0;
   out_1291851874522726762[44] = 0;
   out_1291851874522726762[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_1291851874522726762[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_1291851874522726762[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1291851874522726762[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_1291851874522726762[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_1291851874522726762[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_1291851874522726762[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_1291851874522726762[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_1291851874522726762[53] = -9.8000000000000007*dt;
   out_1291851874522726762[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_1291851874522726762[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_1291851874522726762[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1291851874522726762[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1291851874522726762[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_1291851874522726762[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_1291851874522726762[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_1291851874522726762[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_1291851874522726762[62] = 0;
   out_1291851874522726762[63] = 0;
   out_1291851874522726762[64] = 0;
   out_1291851874522726762[65] = 0;
   out_1291851874522726762[66] = 0;
   out_1291851874522726762[67] = 0;
   out_1291851874522726762[68] = 0;
   out_1291851874522726762[69] = 0;
   out_1291851874522726762[70] = 1;
   out_1291851874522726762[71] = 0;
   out_1291851874522726762[72] = 0;
   out_1291851874522726762[73] = 0;
   out_1291851874522726762[74] = 0;
   out_1291851874522726762[75] = 0;
   out_1291851874522726762[76] = 0;
   out_1291851874522726762[77] = 0;
   out_1291851874522726762[78] = 0;
   out_1291851874522726762[79] = 0;
   out_1291851874522726762[80] = 1;
}
void h_25(double *state, double *unused, double *out_7422228400916633105) {
   out_7422228400916633105[0] = state[6];
}
void H_25(double *state, double *unused, double *out_140850095994413611) {
   out_140850095994413611[0] = 0;
   out_140850095994413611[1] = 0;
   out_140850095994413611[2] = 0;
   out_140850095994413611[3] = 0;
   out_140850095994413611[4] = 0;
   out_140850095994413611[5] = 0;
   out_140850095994413611[6] = 1;
   out_140850095994413611[7] = 0;
   out_140850095994413611[8] = 0;
}
void h_24(double *state, double *unused, double *out_4331290668325946325) {
   out_4331290668325946325[0] = state[4];
   out_4331290668325946325[1] = state[5];
}
void H_24(double *state, double *unused, double *out_7277541537682513282) {
   out_7277541537682513282[0] = 0;
   out_7277541537682513282[1] = 0;
   out_7277541537682513282[2] = 0;
   out_7277541537682513282[3] = 0;
   out_7277541537682513282[4] = 1;
   out_7277541537682513282[5] = 0;
   out_7277541537682513282[6] = 0;
   out_7277541537682513282[7] = 0;
   out_7277541537682513282[8] = 0;
   out_7277541537682513282[9] = 0;
   out_7277541537682513282[10] = 0;
   out_7277541537682513282[11] = 0;
   out_7277541537682513282[12] = 0;
   out_7277541537682513282[13] = 0;
   out_7277541537682513282[14] = 1;
   out_7277541537682513282[15] = 0;
   out_7277541537682513282[16] = 0;
   out_7277541537682513282[17] = 0;
}
void h_30(double *state, double *unused, double *out_7697422463201138994) {
   out_7697422463201138994[0] = state[4];
}
void H_30(double *state, double *unused, double *out_270189043137653681) {
   out_270189043137653681[0] = 0;
   out_270189043137653681[1] = 0;
   out_270189043137653681[2] = 0;
   out_270189043137653681[3] = 0;
   out_270189043137653681[4] = 1;
   out_270189043137653681[5] = 0;
   out_270189043137653681[6] = 0;
   out_270189043137653681[7] = 0;
   out_270189043137653681[8] = 0;
}
void h_26(double *state, double *unused, double *out_1594903027467944604) {
   out_1594903027467944604[0] = state[7];
}
void H_26(double *state, double *unused, double *out_3882353414868469835) {
   out_3882353414868469835[0] = 0;
   out_3882353414868469835[1] = 0;
   out_3882353414868469835[2] = 0;
   out_3882353414868469835[3] = 0;
   out_3882353414868469835[4] = 0;
   out_3882353414868469835[5] = 0;
   out_3882353414868469835[6] = 0;
   out_3882353414868469835[7] = 1;
   out_3882353414868469835[8] = 0;
}
void h_27(double *state, double *unused, double *out_6022842405132072866) {
   out_6022842405132072866[0] = state[3];
}
void H_27(double *state, double *unused, double *out_2444952354938078592) {
   out_2444952354938078592[0] = 0;
   out_2444952354938078592[1] = 0;
   out_2444952354938078592[2] = 0;
   out_2444952354938078592[3] = 1;
   out_2444952354938078592[4] = 0;
   out_2444952354938078592[5] = 0;
   out_2444952354938078592[6] = 0;
   out_2444952354938078592[7] = 0;
   out_2444952354938078592[8] = 0;
}
void h_29(double *state, double *unused, double *out_5747648342847566977) {
   out_5747648342847566977[0] = state[1];
}
void H_29(double *state, double *unused, double *out_240042301176738503) {
   out_240042301176738503[0] = 0;
   out_240042301176738503[1] = 1;
   out_240042301176738503[2] = 0;
   out_240042301176738503[3] = 0;
   out_240042301176738503[4] = 0;
   out_240042301176738503[5] = 0;
   out_240042301176738503[6] = 0;
   out_240042301176738503[7] = 0;
   out_240042301176738503[8] = 0;
}
void h_28(double *state, double *unused, double *out_1021169137471227221) {
   out_1021169137471227221[0] = state[0];
}
void H_28(double *state, double *unused, double *out_2194684810242303374) {
   out_2194684810242303374[0] = 1;
   out_2194684810242303374[1] = 0;
   out_2194684810242303374[2] = 0;
   out_2194684810242303374[3] = 0;
   out_2194684810242303374[4] = 0;
   out_2194684810242303374[5] = 0;
   out_2194684810242303374[6] = 0;
   out_2194684810242303374[7] = 0;
   out_2194684810242303374[8] = 0;
}
void h_31(double *state, double *unused, double *out_8872859793545516272) {
   out_8872859793545516272[0] = state[8];
}
void H_31(double *state, double *unused, double *out_110204134117453183) {
   out_110204134117453183[0] = 0;
   out_110204134117453183[1] = 0;
   out_110204134117453183[2] = 0;
   out_110204134117453183[3] = 0;
   out_110204134117453183[4] = 0;
   out_110204134117453183[5] = 0;
   out_110204134117453183[6] = 0;
   out_110204134117453183[7] = 0;
   out_110204134117453183[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_3495569113160768895) {
  err_fun(nom_x, delta_x, out_3495569113160768895);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_3332561543766221488) {
  inv_err_fun(nom_x, true_x, out_3332561543766221488);
}
void car_H_mod_fun(double *state, double *out_1059310438201751306) {
  H_mod_fun(state, out_1059310438201751306);
}
void car_f_fun(double *state, double dt, double *out_7662000099711208314) {
  f_fun(state,  dt, out_7662000099711208314);
}
void car_F_fun(double *state, double dt, double *out_1291851874522726762) {
  F_fun(state,  dt, out_1291851874522726762);
}
void car_h_25(double *state, double *unused, double *out_7422228400916633105) {
  h_25(state, unused, out_7422228400916633105);
}
void car_H_25(double *state, double *unused, double *out_140850095994413611) {
  H_25(state, unused, out_140850095994413611);
}
void car_h_24(double *state, double *unused, double *out_4331290668325946325) {
  h_24(state, unused, out_4331290668325946325);
}
void car_H_24(double *state, double *unused, double *out_7277541537682513282) {
  H_24(state, unused, out_7277541537682513282);
}
void car_h_30(double *state, double *unused, double *out_7697422463201138994) {
  h_30(state, unused, out_7697422463201138994);
}
void car_H_30(double *state, double *unused, double *out_270189043137653681) {
  H_30(state, unused, out_270189043137653681);
}
void car_h_26(double *state, double *unused, double *out_1594903027467944604) {
  h_26(state, unused, out_1594903027467944604);
}
void car_H_26(double *state, double *unused, double *out_3882353414868469835) {
  H_26(state, unused, out_3882353414868469835);
}
void car_h_27(double *state, double *unused, double *out_6022842405132072866) {
  h_27(state, unused, out_6022842405132072866);
}
void car_H_27(double *state, double *unused, double *out_2444952354938078592) {
  H_27(state, unused, out_2444952354938078592);
}
void car_h_29(double *state, double *unused, double *out_5747648342847566977) {
  h_29(state, unused, out_5747648342847566977);
}
void car_H_29(double *state, double *unused, double *out_240042301176738503) {
  H_29(state, unused, out_240042301176738503);
}
void car_h_28(double *state, double *unused, double *out_1021169137471227221) {
  h_28(state, unused, out_1021169137471227221);
}
void car_H_28(double *state, double *unused, double *out_2194684810242303374) {
  H_28(state, unused, out_2194684810242303374);
}
void car_h_31(double *state, double *unused, double *out_8872859793545516272) {
  h_31(state, unused, out_8872859793545516272);
}
void car_H_31(double *state, double *unused, double *out_110204134117453183) {
  H_31(state, unused, out_110204134117453183);
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

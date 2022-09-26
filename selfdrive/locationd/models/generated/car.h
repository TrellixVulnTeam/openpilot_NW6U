#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_3495569113160768895);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_3332561543766221488);
void car_H_mod_fun(double *state, double *out_1059310438201751306);
void car_f_fun(double *state, double dt, double *out_7662000099711208314);
void car_F_fun(double *state, double dt, double *out_1291851874522726762);
void car_h_25(double *state, double *unused, double *out_7422228400916633105);
void car_H_25(double *state, double *unused, double *out_140850095994413611);
void car_h_24(double *state, double *unused, double *out_4331290668325946325);
void car_H_24(double *state, double *unused, double *out_7277541537682513282);
void car_h_30(double *state, double *unused, double *out_7697422463201138994);
void car_H_30(double *state, double *unused, double *out_270189043137653681);
void car_h_26(double *state, double *unused, double *out_1594903027467944604);
void car_H_26(double *state, double *unused, double *out_3882353414868469835);
void car_h_27(double *state, double *unused, double *out_6022842405132072866);
void car_H_27(double *state, double *unused, double *out_2444952354938078592);
void car_h_29(double *state, double *unused, double *out_5747648342847566977);
void car_H_29(double *state, double *unused, double *out_240042301176738503);
void car_h_28(double *state, double *unused, double *out_1021169137471227221);
void car_H_28(double *state, double *unused, double *out_2194684810242303374);
void car_h_31(double *state, double *unused, double *out_8872859793545516272);
void car_H_31(double *state, double *unused, double *out_110204134117453183);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}
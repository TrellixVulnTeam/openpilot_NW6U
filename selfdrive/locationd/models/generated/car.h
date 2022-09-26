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
void car_err_fun(double *nom_x, double *delta_x, double *out_4853822566850710957);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_7910005515458276405);
void car_H_mod_fun(double *state, double *out_4981524804303897279);
void car_f_fun(double *state, double dt, double *out_3473102566445950304);
void car_F_fun(double *state, double dt, double *out_1612317163615163753);
void car_h_25(double *state, double *unused, double *out_4217284194084537253);
void car_H_25(double *state, double *unused, double *out_7353830765589797144);
void car_h_24(double *state, double *unused, double *out_7320050265028110894);
void car_H_24(double *state, double *unused, double *out_4517341501528236371);
void car_h_30(double *state, double *unused, double *out_7458624653274663427);
void car_H_30(double *state, double *unused, double *out_4176222966628137717);
void car_h_26(double *state, double *unused, double *out_2456951385613817782);
void car_H_26(double *state, double *unused, double *out_3612327446715740920);
void car_h_27(double *state, double *unused, double *out_5160794046986199688);
void car_H_27(double *state, double *unused, double *out_6350986278428562628);
void car_h_29(double *state, double *unused, double *out_1059583135860814281);
void car_H_29(double *state, double *unused, double *out_3665991622313745533);
void car_h_28(double *state, double *unused, double *out_4557478162309722171);
void car_H_28(double *state, double *unused, double *out_8748390639383276107);
void car_h_31(double *state, double *unused, double *out_3675655300799986106);
void car_H_31(double *state, double *unused, double *out_7384476727466757572);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}
#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_1220311343949854590);
void live_err_fun(double *nom_x, double *delta_x, double *out_3153226439935894210);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_8463461953737375796);
void live_H_mod_fun(double *state, double *out_21562585824775244);
void live_f_fun(double *state, double dt, double *out_8333998904903447252);
void live_F_fun(double *state, double dt, double *out_7893806645490710208);
void live_h_4(double *state, double *unused, double *out_7403799576841385348);
void live_H_4(double *state, double *unused, double *out_3340629146868555404);
void live_h_9(double *state, double *unused, double *out_3237950739396852411);
void live_H_9(double *state, double *unused, double *out_451767594588476062);
void live_h_10(double *state, double *unused, double *out_5171938564839258654);
void live_H_10(double *state, double *unused, double *out_56454275011937580);
void live_h_12(double *state, double *unused, double *out_4968959824595442666);
void live_H_12(double *state, double *unused, double *out_2719530121820961737);
void live_h_31(double *state, double *unused, double *out_4693443347889562258);
void live_H_31(double *state, double *unused, double *out_26032910504051972);
void live_h_32(double *state, double *unused, double *out_6872470075943617214);
void live_H_32(double *state, double *unused, double *out_7826837534583484090);
void live_h_13(double *state, double *unused, double *out_1840851040610273547);
void live_H_13(double *state, double *unused, double *out_1182368437291281484);
void live_h_14(double *state, double *unused, double *out_3237950739396852411);
void live_H_14(double *state, double *unused, double *out_451767594588476062);
void live_h_33(double *state, double *unused, double *out_903316753777834409);
void live_H_33(double *state, double *unused, double *out_3176589915142909576);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}
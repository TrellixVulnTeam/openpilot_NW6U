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
void live_H(double *in_vec, double *out_7309746982530501858);
void live_err_fun(double *nom_x, double *delta_x, double *out_5779196364896740181);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_4961782207135168051);
void live_H_mod_fun(double *state, double *out_7146165817346037743);
void live_f_fun(double *state, double dt, double *out_2262751468361150474);
void live_F_fun(double *state, double dt, double *out_7749511526562044438);
void live_h_4(double *state, double *unused, double *out_9007062715361016091);
void live_H_4(double *state, double *unused, double *out_7471456318502086351);
void live_h_9(double *state, double *unused, double *out_9115699719530867729);
void live_H_9(double *state, double *unused, double *out_184237383237638881);
void live_h_10(double *state, double *unused, double *out_2755636301544987966);
void live_H_10(double *state, double *unused, double *out_7945717824255188420);
void live_h_12(double *state, double *unused, double *out_1783803183341343539);
void live_H_12(double *state, double *unused, double *out_4594029378164732269);
void live_h_31(double *state, double *unused, double *out_4590535298572725627);
void live_H_31(double *state, double *unused, double *out_293563121854889153);
void live_h_32(double *state, double *unused, double *out_4586094836937765954);
void live_H_32(double *state, double *unused, double *out_8323493498801388045);
void live_h_13(double *state, double *unused, double *out_6216853385747818522);
void live_H_13(double *state, double *unused, double *out_3463561509396612115);
void live_h_14(double *state, double *unused, double *out_9115699719530867729);
void live_H_14(double *state, double *unused, double *out_184237383237638881);
void live_h_33(double *state, double *unused, double *out_1564424883571661231);
void live_H_33(double *state, double *unused, double *out_954237256490621371);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}
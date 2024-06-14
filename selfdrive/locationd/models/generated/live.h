#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_35(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_1525046797953799792);
void live_err_fun(double *nom_x, double *delta_x, double *out_6787273198622262185);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_1163963483344490897);
void live_H_mod_fun(double *state, double *out_2148830550237504684);
void live_f_fun(double *state, double dt, double *out_9175569816144614590);
void live_F_fun(double *state, double dt, double *out_6638881892447955891);
void live_h_4(double *state, double *unused, double *out_4748581981128319294);
void live_H_4(double *state, double *unused, double *out_8350786998460166632);
void live_h_9(double *state, double *unused, double *out_3299619827359850763);
void live_H_9(double *state, double *unused, double *out_2808738139984937514);
void live_h_10(double *state, double *unused, double *out_6759185635745913961);
void live_H_10(double *state, double *unused, double *out_8931560196458308531);
void live_h_12(double *state, double *unused, double *out_3977376211679548972);
void live_H_12(double *state, double *unused, double *out_8971886023507760299);
void live_h_35(double *state, double *unused, double *out_6053256877303408526);
void live_H_35(double *state, double *unused, double *out_6729295017876777608);
void live_h_32(double *state, double *unused, double *out_8660421758678747156);
void live_H_32(double *state, double *unused, double *out_6598130232936037139);
void live_h_13(double *state, double *unused, double *out_1407584935076879677);
void live_H_13(double *state, double *unused, double *out_2297943131464073829);
void live_h_14(double *state, double *unused, double *out_3299619827359850763);
void live_H_14(double *state, double *unused, double *out_2808738139984937514);
void live_h_33(double *state, double *unused, double *out_2880611375828088826);
void live_H_33(double *state, double *unused, double *out_3578738013237920004);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}
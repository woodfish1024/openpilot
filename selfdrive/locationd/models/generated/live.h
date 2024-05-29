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
void live_H(double *in_vec, double *out_5397508461351933175);
void live_err_fun(double *nom_x, double *delta_x, double *out_6616526752652423851);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_489500073009557101);
void live_H_mod_fun(double *state, double *out_2299434493560280780);
void live_f_fun(double *state, double dt, double *out_8025762527095415664);
void live_F_fun(double *state, double dt, double *out_1558757361660199893);
void live_h_4(double *state, double *unused, double *out_477891168818797768);
void live_H_4(double *state, double *unused, double *out_7679088972560206961);
void live_h_9(double *state, double *unused, double *out_7897025827656141437);
void live_H_9(double *state, double *unused, double *out_7878793548869265313);
void live_h_10(double *state, double *unused, double *out_6864998424677469431);
void live_H_10(double *state, double *unused, double *out_6005479821192351653);
void live_h_12(double *state, double *unused, double *out_1092987823989953485);
void live_H_12(double *state, double *unused, double *out_8300187997607800628);
void live_h_35(double *state, double *unused, double *out_4516492124335460556);
void live_H_35(double *state, double *unused, double *out_7400993043776737279);
void live_h_32(double *state, double *unused, double *out_4951607991013811668);
void live_H_32(double *state, double *unused, double *out_5446340617065732104);
void live_h_13(double *state, double *unused, double *out_2532856805632969062);
void live_H_13(double *state, double *unused, double *out_4026559930567762347);
void live_h_14(double *state, double *unused, double *out_7897025827656141437);
void live_H_14(double *state, double *unused, double *out_7878793548869265313);
void live_h_33(double *state, double *unused, double *out_2979582365054270509);
void live_H_33(double *state, double *unused, double *out_4250436039137879675);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}
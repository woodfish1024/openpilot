#pragma once
#include "rednose/helpers/ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_5209869411261452488);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_1703827342772737256);
void car_H_mod_fun(double *state, double *out_7039949086143635320);
void car_f_fun(double *state, double dt, double *out_3112150461183044984);
void car_F_fun(double *state, double dt, double *out_3917161440273399707);
void car_h_25(double *state, double *unused, double *out_2376580926527384629);
void car_H_25(double *state, double *unused, double *out_7719553694858263201);
void car_h_24(double *state, double *unused, double *out_4215534704440492819);
void car_H_24(double *state, double *unused, double *out_2850738829830556349);
void car_h_30(double *state, double *unused, double *out_2101386864242878740);
void car_H_30(double *state, double *unused, double *out_3191857364730655003);
void car_h_26(double *state, double *unused, double *out_4297654653170970406);
void car_H_26(double *state, double *unused, double *out_3978050375984206977);
void car_h_27(double *state, double *unused, double *out_1773265041850907112);
void car_H_27(double *state, double *unused, double *out_1017094052930230092);
void car_h_29(double *state, double *unused, double *out_2900286403417966905);
void car_H_29(double *state, double *unused, double *out_3702088709045047187);
void car_h_28(double *state, double *unused, double *out_2716774894752569547);
void car_H_28(double *state, double *unused, double *out_5665718980610373438);
void car_h_31(double *state, double *unused, double *out_5324306916882869590);
void car_H_31(double *state, double *unused, double *out_3351842273750855501);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}
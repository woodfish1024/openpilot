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
void car_err_fun(double *nom_x, double *delta_x, double *out_2817899192381410803);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_7343685783341145042);
void car_H_mod_fun(double *state, double *out_4599531654129897066);
void car_f_fun(double *state, double dt, double *out_442343507582822689);
void car_F_fun(double *state, double dt, double *out_9136000802478021924);
void car_h_25(double *state, double *unused, double *out_7170330515616841330);
void car_H_25(double *state, double *unused, double *out_1204733297232033761);
void car_h_24(double *state, double *unused, double *out_2924262862368423500);
void car_H_24(double *state, double *unused, double *out_2353192301570289375);
void car_h_30(double *state, double *unused, double *out_1382334063405462395);
void car_H_30(double *state, double *unused, double *out_1313599661275214866);
void car_h_26(double *state, double *unused, double *out_1654451988038870290);
void car_H_26(double *state, double *unused, double *out_4946236616106089985);
void car_h_27(double *state, double *unused, double *out_6567014630940363813);
void car_H_27(double *state, double *unused, double *out_3537193732459158083);
void car_h_29(double *state, double *unused, double *out_6724201299291492958);
void car_H_29(double *state, double *unused, double *out_1823831005589607050);
void car_h_28(double *state, double *unused, double *out_8842943391324855118);
void car_H_28(double *state, double *unused, double *out_3258568011479923524);
void car_h_31(double *state, double *unused, double *out_3072027217337469466);
void car_H_31(double *state, double *unused, double *out_5572444718339441461);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}
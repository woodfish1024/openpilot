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
 *                       Code generated with SymPy 1.12                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_5209869411261452488) {
   out_5209869411261452488[0] = delta_x[0] + nom_x[0];
   out_5209869411261452488[1] = delta_x[1] + nom_x[1];
   out_5209869411261452488[2] = delta_x[2] + nom_x[2];
   out_5209869411261452488[3] = delta_x[3] + nom_x[3];
   out_5209869411261452488[4] = delta_x[4] + nom_x[4];
   out_5209869411261452488[5] = delta_x[5] + nom_x[5];
   out_5209869411261452488[6] = delta_x[6] + nom_x[6];
   out_5209869411261452488[7] = delta_x[7] + nom_x[7];
   out_5209869411261452488[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_1703827342772737256) {
   out_1703827342772737256[0] = -nom_x[0] + true_x[0];
   out_1703827342772737256[1] = -nom_x[1] + true_x[1];
   out_1703827342772737256[2] = -nom_x[2] + true_x[2];
   out_1703827342772737256[3] = -nom_x[3] + true_x[3];
   out_1703827342772737256[4] = -nom_x[4] + true_x[4];
   out_1703827342772737256[5] = -nom_x[5] + true_x[5];
   out_1703827342772737256[6] = -nom_x[6] + true_x[6];
   out_1703827342772737256[7] = -nom_x[7] + true_x[7];
   out_1703827342772737256[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_7039949086143635320) {
   out_7039949086143635320[0] = 1.0;
   out_7039949086143635320[1] = 0;
   out_7039949086143635320[2] = 0;
   out_7039949086143635320[3] = 0;
   out_7039949086143635320[4] = 0;
   out_7039949086143635320[5] = 0;
   out_7039949086143635320[6] = 0;
   out_7039949086143635320[7] = 0;
   out_7039949086143635320[8] = 0;
   out_7039949086143635320[9] = 0;
   out_7039949086143635320[10] = 1.0;
   out_7039949086143635320[11] = 0;
   out_7039949086143635320[12] = 0;
   out_7039949086143635320[13] = 0;
   out_7039949086143635320[14] = 0;
   out_7039949086143635320[15] = 0;
   out_7039949086143635320[16] = 0;
   out_7039949086143635320[17] = 0;
   out_7039949086143635320[18] = 0;
   out_7039949086143635320[19] = 0;
   out_7039949086143635320[20] = 1.0;
   out_7039949086143635320[21] = 0;
   out_7039949086143635320[22] = 0;
   out_7039949086143635320[23] = 0;
   out_7039949086143635320[24] = 0;
   out_7039949086143635320[25] = 0;
   out_7039949086143635320[26] = 0;
   out_7039949086143635320[27] = 0;
   out_7039949086143635320[28] = 0;
   out_7039949086143635320[29] = 0;
   out_7039949086143635320[30] = 1.0;
   out_7039949086143635320[31] = 0;
   out_7039949086143635320[32] = 0;
   out_7039949086143635320[33] = 0;
   out_7039949086143635320[34] = 0;
   out_7039949086143635320[35] = 0;
   out_7039949086143635320[36] = 0;
   out_7039949086143635320[37] = 0;
   out_7039949086143635320[38] = 0;
   out_7039949086143635320[39] = 0;
   out_7039949086143635320[40] = 1.0;
   out_7039949086143635320[41] = 0;
   out_7039949086143635320[42] = 0;
   out_7039949086143635320[43] = 0;
   out_7039949086143635320[44] = 0;
   out_7039949086143635320[45] = 0;
   out_7039949086143635320[46] = 0;
   out_7039949086143635320[47] = 0;
   out_7039949086143635320[48] = 0;
   out_7039949086143635320[49] = 0;
   out_7039949086143635320[50] = 1.0;
   out_7039949086143635320[51] = 0;
   out_7039949086143635320[52] = 0;
   out_7039949086143635320[53] = 0;
   out_7039949086143635320[54] = 0;
   out_7039949086143635320[55] = 0;
   out_7039949086143635320[56] = 0;
   out_7039949086143635320[57] = 0;
   out_7039949086143635320[58] = 0;
   out_7039949086143635320[59] = 0;
   out_7039949086143635320[60] = 1.0;
   out_7039949086143635320[61] = 0;
   out_7039949086143635320[62] = 0;
   out_7039949086143635320[63] = 0;
   out_7039949086143635320[64] = 0;
   out_7039949086143635320[65] = 0;
   out_7039949086143635320[66] = 0;
   out_7039949086143635320[67] = 0;
   out_7039949086143635320[68] = 0;
   out_7039949086143635320[69] = 0;
   out_7039949086143635320[70] = 1.0;
   out_7039949086143635320[71] = 0;
   out_7039949086143635320[72] = 0;
   out_7039949086143635320[73] = 0;
   out_7039949086143635320[74] = 0;
   out_7039949086143635320[75] = 0;
   out_7039949086143635320[76] = 0;
   out_7039949086143635320[77] = 0;
   out_7039949086143635320[78] = 0;
   out_7039949086143635320[79] = 0;
   out_7039949086143635320[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_3112150461183044984) {
   out_3112150461183044984[0] = state[0];
   out_3112150461183044984[1] = state[1];
   out_3112150461183044984[2] = state[2];
   out_3112150461183044984[3] = state[3];
   out_3112150461183044984[4] = state[4];
   out_3112150461183044984[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_3112150461183044984[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_3112150461183044984[7] = state[7];
   out_3112150461183044984[8] = state[8];
}
void F_fun(double *state, double dt, double *out_3917161440273399707) {
   out_3917161440273399707[0] = 1;
   out_3917161440273399707[1] = 0;
   out_3917161440273399707[2] = 0;
   out_3917161440273399707[3] = 0;
   out_3917161440273399707[4] = 0;
   out_3917161440273399707[5] = 0;
   out_3917161440273399707[6] = 0;
   out_3917161440273399707[7] = 0;
   out_3917161440273399707[8] = 0;
   out_3917161440273399707[9] = 0;
   out_3917161440273399707[10] = 1;
   out_3917161440273399707[11] = 0;
   out_3917161440273399707[12] = 0;
   out_3917161440273399707[13] = 0;
   out_3917161440273399707[14] = 0;
   out_3917161440273399707[15] = 0;
   out_3917161440273399707[16] = 0;
   out_3917161440273399707[17] = 0;
   out_3917161440273399707[18] = 0;
   out_3917161440273399707[19] = 0;
   out_3917161440273399707[20] = 1;
   out_3917161440273399707[21] = 0;
   out_3917161440273399707[22] = 0;
   out_3917161440273399707[23] = 0;
   out_3917161440273399707[24] = 0;
   out_3917161440273399707[25] = 0;
   out_3917161440273399707[26] = 0;
   out_3917161440273399707[27] = 0;
   out_3917161440273399707[28] = 0;
   out_3917161440273399707[29] = 0;
   out_3917161440273399707[30] = 1;
   out_3917161440273399707[31] = 0;
   out_3917161440273399707[32] = 0;
   out_3917161440273399707[33] = 0;
   out_3917161440273399707[34] = 0;
   out_3917161440273399707[35] = 0;
   out_3917161440273399707[36] = 0;
   out_3917161440273399707[37] = 0;
   out_3917161440273399707[38] = 0;
   out_3917161440273399707[39] = 0;
   out_3917161440273399707[40] = 1;
   out_3917161440273399707[41] = 0;
   out_3917161440273399707[42] = 0;
   out_3917161440273399707[43] = 0;
   out_3917161440273399707[44] = 0;
   out_3917161440273399707[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_3917161440273399707[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_3917161440273399707[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3917161440273399707[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_3917161440273399707[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_3917161440273399707[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_3917161440273399707[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_3917161440273399707[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_3917161440273399707[53] = -9.8000000000000007*dt;
   out_3917161440273399707[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_3917161440273399707[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_3917161440273399707[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3917161440273399707[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3917161440273399707[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_3917161440273399707[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_3917161440273399707[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_3917161440273399707[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_3917161440273399707[62] = 0;
   out_3917161440273399707[63] = 0;
   out_3917161440273399707[64] = 0;
   out_3917161440273399707[65] = 0;
   out_3917161440273399707[66] = 0;
   out_3917161440273399707[67] = 0;
   out_3917161440273399707[68] = 0;
   out_3917161440273399707[69] = 0;
   out_3917161440273399707[70] = 1;
   out_3917161440273399707[71] = 0;
   out_3917161440273399707[72] = 0;
   out_3917161440273399707[73] = 0;
   out_3917161440273399707[74] = 0;
   out_3917161440273399707[75] = 0;
   out_3917161440273399707[76] = 0;
   out_3917161440273399707[77] = 0;
   out_3917161440273399707[78] = 0;
   out_3917161440273399707[79] = 0;
   out_3917161440273399707[80] = 1;
}
void h_25(double *state, double *unused, double *out_2376580926527384629) {
   out_2376580926527384629[0] = state[6];
}
void H_25(double *state, double *unused, double *out_7719553694858263201) {
   out_7719553694858263201[0] = 0;
   out_7719553694858263201[1] = 0;
   out_7719553694858263201[2] = 0;
   out_7719553694858263201[3] = 0;
   out_7719553694858263201[4] = 0;
   out_7719553694858263201[5] = 0;
   out_7719553694858263201[6] = 1;
   out_7719553694858263201[7] = 0;
   out_7719553694858263201[8] = 0;
}
void h_24(double *state, double *unused, double *out_4215534704440492819) {
   out_4215534704440492819[0] = state[4];
   out_4215534704440492819[1] = state[5];
}
void H_24(double *state, double *unused, double *out_2850738829830556349) {
   out_2850738829830556349[0] = 0;
   out_2850738829830556349[1] = 0;
   out_2850738829830556349[2] = 0;
   out_2850738829830556349[3] = 0;
   out_2850738829830556349[4] = 1;
   out_2850738829830556349[5] = 0;
   out_2850738829830556349[6] = 0;
   out_2850738829830556349[7] = 0;
   out_2850738829830556349[8] = 0;
   out_2850738829830556349[9] = 0;
   out_2850738829830556349[10] = 0;
   out_2850738829830556349[11] = 0;
   out_2850738829830556349[12] = 0;
   out_2850738829830556349[13] = 0;
   out_2850738829830556349[14] = 1;
   out_2850738829830556349[15] = 0;
   out_2850738829830556349[16] = 0;
   out_2850738829830556349[17] = 0;
}
void h_30(double *state, double *unused, double *out_2101386864242878740) {
   out_2101386864242878740[0] = state[4];
}
void H_30(double *state, double *unused, double *out_3191857364730655003) {
   out_3191857364730655003[0] = 0;
   out_3191857364730655003[1] = 0;
   out_3191857364730655003[2] = 0;
   out_3191857364730655003[3] = 0;
   out_3191857364730655003[4] = 1;
   out_3191857364730655003[5] = 0;
   out_3191857364730655003[6] = 0;
   out_3191857364730655003[7] = 0;
   out_3191857364730655003[8] = 0;
}
void h_26(double *state, double *unused, double *out_4297654653170970406) {
   out_4297654653170970406[0] = state[7];
}
void H_26(double *state, double *unused, double *out_3978050375984206977) {
   out_3978050375984206977[0] = 0;
   out_3978050375984206977[1] = 0;
   out_3978050375984206977[2] = 0;
   out_3978050375984206977[3] = 0;
   out_3978050375984206977[4] = 0;
   out_3978050375984206977[5] = 0;
   out_3978050375984206977[6] = 0;
   out_3978050375984206977[7] = 1;
   out_3978050375984206977[8] = 0;
}
void h_27(double *state, double *unused, double *out_1773265041850907112) {
   out_1773265041850907112[0] = state[3];
}
void H_27(double *state, double *unused, double *out_1017094052930230092) {
   out_1017094052930230092[0] = 0;
   out_1017094052930230092[1] = 0;
   out_1017094052930230092[2] = 0;
   out_1017094052930230092[3] = 1;
   out_1017094052930230092[4] = 0;
   out_1017094052930230092[5] = 0;
   out_1017094052930230092[6] = 0;
   out_1017094052930230092[7] = 0;
   out_1017094052930230092[8] = 0;
}
void h_29(double *state, double *unused, double *out_2900286403417966905) {
   out_2900286403417966905[0] = state[1];
}
void H_29(double *state, double *unused, double *out_3702088709045047187) {
   out_3702088709045047187[0] = 0;
   out_3702088709045047187[1] = 1;
   out_3702088709045047187[2] = 0;
   out_3702088709045047187[3] = 0;
   out_3702088709045047187[4] = 0;
   out_3702088709045047187[5] = 0;
   out_3702088709045047187[6] = 0;
   out_3702088709045047187[7] = 0;
   out_3702088709045047187[8] = 0;
}
void h_28(double *state, double *unused, double *out_2716774894752569547) {
   out_2716774894752569547[0] = state[0];
}
void H_28(double *state, double *unused, double *out_5665718980610373438) {
   out_5665718980610373438[0] = 1;
   out_5665718980610373438[1] = 0;
   out_5665718980610373438[2] = 0;
   out_5665718980610373438[3] = 0;
   out_5665718980610373438[4] = 0;
   out_5665718980610373438[5] = 0;
   out_5665718980610373438[6] = 0;
   out_5665718980610373438[7] = 0;
   out_5665718980610373438[8] = 0;
}
void h_31(double *state, double *unused, double *out_5324306916882869590) {
   out_5324306916882869590[0] = state[8];
}
void H_31(double *state, double *unused, double *out_3351842273750855501) {
   out_3351842273750855501[0] = 0;
   out_3351842273750855501[1] = 0;
   out_3351842273750855501[2] = 0;
   out_3351842273750855501[3] = 0;
   out_3351842273750855501[4] = 0;
   out_3351842273750855501[5] = 0;
   out_3351842273750855501[6] = 0;
   out_3351842273750855501[7] = 0;
   out_3351842273750855501[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_5209869411261452488) {
  err_fun(nom_x, delta_x, out_5209869411261452488);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_1703827342772737256) {
  inv_err_fun(nom_x, true_x, out_1703827342772737256);
}
void car_H_mod_fun(double *state, double *out_7039949086143635320) {
  H_mod_fun(state, out_7039949086143635320);
}
void car_f_fun(double *state, double dt, double *out_3112150461183044984) {
  f_fun(state,  dt, out_3112150461183044984);
}
void car_F_fun(double *state, double dt, double *out_3917161440273399707) {
  F_fun(state,  dt, out_3917161440273399707);
}
void car_h_25(double *state, double *unused, double *out_2376580926527384629) {
  h_25(state, unused, out_2376580926527384629);
}
void car_H_25(double *state, double *unused, double *out_7719553694858263201) {
  H_25(state, unused, out_7719553694858263201);
}
void car_h_24(double *state, double *unused, double *out_4215534704440492819) {
  h_24(state, unused, out_4215534704440492819);
}
void car_H_24(double *state, double *unused, double *out_2850738829830556349) {
  H_24(state, unused, out_2850738829830556349);
}
void car_h_30(double *state, double *unused, double *out_2101386864242878740) {
  h_30(state, unused, out_2101386864242878740);
}
void car_H_30(double *state, double *unused, double *out_3191857364730655003) {
  H_30(state, unused, out_3191857364730655003);
}
void car_h_26(double *state, double *unused, double *out_4297654653170970406) {
  h_26(state, unused, out_4297654653170970406);
}
void car_H_26(double *state, double *unused, double *out_3978050375984206977) {
  H_26(state, unused, out_3978050375984206977);
}
void car_h_27(double *state, double *unused, double *out_1773265041850907112) {
  h_27(state, unused, out_1773265041850907112);
}
void car_H_27(double *state, double *unused, double *out_1017094052930230092) {
  H_27(state, unused, out_1017094052930230092);
}
void car_h_29(double *state, double *unused, double *out_2900286403417966905) {
  h_29(state, unused, out_2900286403417966905);
}
void car_H_29(double *state, double *unused, double *out_3702088709045047187) {
  H_29(state, unused, out_3702088709045047187);
}
void car_h_28(double *state, double *unused, double *out_2716774894752569547) {
  h_28(state, unused, out_2716774894752569547);
}
void car_H_28(double *state, double *unused, double *out_5665718980610373438) {
  H_28(state, unused, out_5665718980610373438);
}
void car_h_31(double *state, double *unused, double *out_5324306916882869590) {
  h_31(state, unused, out_5324306916882869590);
}
void car_H_31(double *state, double *unused, double *out_3351842273750855501) {
  H_31(state, unused, out_3351842273750855501);
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

ekf_lib_init(car)

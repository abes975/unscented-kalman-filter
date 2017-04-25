#include <iostream>
#include <cfloat>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}


VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  if (!estimations.size() || (estimations.size() != ground_truth.size()))
    return rmse;

  for (int i = 0; i < estimations.size(); ++i) {
      VectorXd res = estimations.at(i) - ground_truth.at(i);
      res = res.array() * res.array();
      rmse += res;
  }
  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}

VectorXd Tools::PolarToCartesian(const VectorXd& polar) {
    VectorXd cartesian(4);

    double rho = polar(0);
    double phi = polar(1);
    // velocity
    double rho_v = polar(2);

    double px = rho * cos(phi);
    double py = rho * sin(phi);
    double vx = rho_v * cos(phi);
    double vy = rho_v * sin(phi);

    cartesian << px, py, vx, vy;
    return cartesian;
}

VectorXd Tools::CartesianToPolar(const VectorXd& cartesian) {
    float px = cartesian(0);
    float py = cartesian(1);
    float vx = cartesian(2);
    float vy = cartesian(3);
    double rho = sqrt(px*px+py*py);

    VectorXd polar(3);
    polar << 0,0,0;

    if (!px) {
        px = FLT_EPSILON;
        rho = sqrt(px*px+py*py);
    }
    if (!py) {
        py = FLT_EPSILON;
        rho = sqrt(px*px+py*py);
    }

    if (!rho) {
        rho = FLT_EPSILON;
    }

    double phi = atan2(py, px);
    double rho_v = (px*vx+py*vy)/rho;
    polar << rho,
             phi,
             rho_v;

    return polar;
}

double Tools::normalize_angle(double to_normalize) {
  while (to_normalize > M_PI)
    to_normalize -= 2. * M_PI;
  while (to_normalize < -M_PI)
    to_normalize += 2. * M_PI;
  return to_normalize;
}

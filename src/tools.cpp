#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    
    VectorXd sq_diff = VectorXd(4); // px, py, vx, vy
    sq_diff.fill(0.0);
    
    // border coniditions
    if(estimations.size() != ground_truth.size()) {
      return sq_diff;
    }
    
    for(int i=0; i < estimations.size(); i++) {
      VectorXd diff = estimations[i] - ground_truth[i];
      sq_diff += diff.square();
    }

    // find mean
    VectorXd sq_diff_mean = sq_diff/estimations.size();

    return sq_diff_mean.sqrt();
}

VectorXd convert_to_cartesian(const VectorXd& polar) {
  VectorXd cartesian = VectorXd(2);
  double px = -polar(0)*sin(polar(1));
  double py = polar(0)*cos(polar(1));
  cartesian << px, py;
  return cartesian;
}

VectorXd convert_to_polar(const VectorXd& cartesian) {
  // 
  double px = cartesian(0);
  double py = cartesian(1);
  double v = cartesian(2);
  double yaw = cartesian(3);
  
  double rho = sqrt(px*px+ py*py);
  double phi = atan2(py, px);
  while(phi < -M_PI && phi > M_PI) phi < 0 ? phi += 2*M_PI : phi -= 2*M_PI;

  double rhodot = (px*cos(yaw)*v + py*sin(yaw)*v)/rho;

  VectorXd polar = VectorXd(3);
  polar << rho, phi, rhod;
  return polar;
}
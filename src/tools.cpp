#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

double Tools::percentile_check(const vector<double> &dist, double percentile_value) {
  ssize_t total = dist.size();
  ssize_t samples = 0;
  for(auto value : dist) {
    if(value > percentile_value) {
      samples++;
    } 
  }

  return (double(samples)/total);
}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    
    VectorXd sq_diff = VectorXd(4); // px, py, vx, vy
    sq_diff.fill(0.0);

    ssize_t length = estimations.size();
    if(length == 0) {
      return sq_diff;
    }
    
    // border coniditions
    if(length != ground_truth.size()) {
      return sq_diff;
    }
    
    for(int i=0; i < length; i++) {
      VectorXd diff = estimations[i] - ground_truth[i];
      diff = diff.array()*diff.array();
      sq_diff += diff;
    }

    // find mean
    VectorXd sq_diff_mean = sq_diff/length;
    VectorXd rmse = sq_diff_mean.array().sqrt();

    return rmse;
}

VectorXd Tools::convert_to_cartesian(const VectorXd& polar) {
  VectorXd cartesian = VectorXd(2);
  double px = -polar(0)*sin(polar(1));
  double py = polar(0)*cos(polar(1));
  cartesian << px, py;
  return cartesian;
}

VectorXd Tools::convert_to_polar(const VectorXd& cartesian) {
  // 
  double px = cartesian(0);
  double py = cartesian(1);
  double v = cartesian(2);
  double yaw = cartesian(3);

  double rho = sqrt(px*px+ py*py);
  double phi = 0;
  double rhod = 0.0;
  VectorXd polar = VectorXd(3);
  if (fabs(px) < 0.0001 && fabs(py) < 0.0001) {
    polar << rho, phi, rhod;
    return polar;
  }
  phi = atan2(py, px);
  if (fabs(rho) < 0.0001) {
    polar << rho, phi, rhod;
    return polar;
  }
  rhod = (px*cos(yaw)*v + py*sin(yaw)*v)/rho;

  polar << rho, phi, rhod;
  return polar;
}

double Tools::phi_range(double phi) {
  double b = phi;
  while(phi < -M_PI) phi += 2*M_PI;
  while(phi > M_PI) phi -= 2*M_PI;
  if (fabs(b - phi) > 0.001) {
    std::cout << "b: " << b << " a: " << phi << std::endl;
  }
  return phi;
}
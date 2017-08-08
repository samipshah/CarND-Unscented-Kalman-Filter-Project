#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

bool debug = true;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5; // may be less than 1 for a bicyclist

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.1; 

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  n_x_ = 5;

  n_aug_ = 7; // added process noise elements

  lambda_ = 3 - n_aug_;

  Xsig_pred_ = MatrixXd(5, 2*n_aug_+1);
  Xsig_pred_.fill(0.0);

  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_/(lambda_+n_aug_);
  for(int i=1;i<(2*n_aug_+1);i++) {
      weights_(i) = .5/(lambda_+n_aug_);
  }

  R_ = MatrixXd(3,3);
  R_.fill(0.0);
  R_(0,0) = std_radr_*std_radr_;
  R_(1,1) = std_radphi_*std_radphi_;
  R_(2,2) = std_radrd_*std_radrd_;

  L_ = MatrixXd(2,2);
  L_.fill(0.0);
  L_(0,0) = std_laspx_*std_laspx_;
  L_(1,1) = std_laspy_*std_laspy_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  // initialize
  if(!is_initialized_) {
    // initialize
    time_us_ = meas_package.timestamp_;
    VectorXd cartesian = VectorXd(2);
    if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // Laser
      cartesian << meas_package.raw_measurements_.head(2);
    } else {
      // Radar
      cartesian = tools_.convert_to_cartesian(meas_package.raw_measurements_.head(3));
    }
    
    x_.fill(0.0);
    x_.head(2) = cartesian;
    x_(2) = 3; // approx 3m/s velocity
    x_(3) = 0; 
    x_(4) = 0.5; 

    // initialize P_
    P_.fill(0.0);
    for (int i=0; i<x_.size(); i++) {
      P_(i,i) = 1;
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  // prediction
  double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0; // convert to seconds
  time_us_ = meas_package.timestamp_;
  Prediction(delta_t);

  // check sensor type
  if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
    if(!use_laser_) {
      return;
    }
    // update
    UpdateLidar(meas_package);

  } else {
    if(!use_radar_) {
      return;
    }
    // update  
    UpdateRadar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // use x_ and find augmented sigma points using augmented process covariance.
  MatrixXd P_aug = MatrixXd(n_x_ + 2, n_x_ + 2);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  MatrixXd sqrt_P = P_aug.llt().matrixL();
  double c = sqrt(lambda_ + n_aug_);

  MatrixXd Xsig_aug  = MatrixXd(7, 2*n_aug_ + 1);
  VectorXd x_aug = VectorXd(n_x_ + 2);
  x_aug.fill(0.0);
  x_aug.head(5) = x_;
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i+1)      = x_aug + c * sqrt_P.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - c * sqrt_P.col(i);
  }

  // predict all sigma points to next state px, py using velocity, delta_t, and yaw and yaw rate
  // and additive process noise 
  double sq_delta_t = delta_t*delta_t;
  VectorXd process_noise = VectorXd(5);
  VectorXd state_transition = VectorXd(5);
  
  //predict sigma points
  for(int i=0; i<(2*n_aug_+1); i++) {
      process_noise.fill(0.0);
      // noise calculation 
      double vk = Xsig_aug.col(i)(2);
      double yawk = Xsig_aug.col(i)(3);
      double yawdk = Xsig_aug.col(i)(4);
      double vak = Xsig_aug.col(i)(5);
      double vyawddk = Xsig_aug.col(i)(6);
      if (vak != 0.0) {
          process_noise(2) = delta_t*vak;
          process_noise(0) = .5*sq_delta_t*cos(yawk)*vak;
          process_noise(1) = .5*sq_delta_t*sin(yawk)*vak;
      }
      
      if (vyawddk != 0.0) {
          process_noise(3) = .5*sq_delta_t*vyawddk;
          process_noise(4) = delta_t*vyawddk;
      }

      // state transition vector
      state_transition.fill(0.0);
      if (fabs(yawdk) > 0.001) {
        state_transition(0) = (vk/yawdk)*(sin(yawk+yawdk*delta_t)-sin(yawk));
        state_transition(1) = (vk/yawdk)*(-cos(yawk+yawdk*delta_t)+cos(yawk));
        state_transition(3) = yawdk*delta_t;
      } else {
        state_transition(0) = vk*cos(yawk)*delta_t;
        state_transition(1) = vk*sin(yawk)*delta_t;
      }

      Xsig_pred_.col(i) = Xsig_aug.col(i).head(5) + state_transition + process_noise;
      Xsig_pred_.col(i)(4) = tools_.phi_range(Xsig_pred_.col(i)(4));
  }

  // calculate weights and find mean from these predicted sigma points for the one next
  // state vector
  
  // predict state by finding mean of all the sigma points prediction
  x_ = weights_.transpose()*Xsig_pred_.transpose();
  x_(4) = tools_.phi_range(x_(4));

  // predict covariance matrix for process by finding mean variance of each of the 
  // sigma points from mean
  P_.fill(0.0);
  for(int i=0;i<(2*n_aug_+1);i++) {
      VectorXd diff = Xsig_pred_.col(i) - x_;
      diff(3) = tools_.phi_range(diff(3));
      P_ += weights_(i)*(diff*diff.transpose());
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  VectorXd z = meas_package.raw_measurements_.head(2);
  MatrixXd Zsig = MatrixXd(2, 2*n_aug_+1);
  for(int i=0; i<(2*n_aug_+1);i++) {
    Zsig.col(i) = Xsig_pred_.col(i).head(2);  
  }
  VectorXd z_pred = x_.head(2);

  // MatrixXd S = MatrixXd(2,2);
  // S.fill(0.0);
  // for(int i=0; i<(2*n_aug_+1); i++) {
  //   VectorXd diff = Zsig.col(i) - z_pred;
  //   S += weights_(i)*(diff*diff.transpose());
  // }
  // MatrixXd S = P_.topLeftCorner(2,2);

  MatrixXd S = MatrixXd(2,2);
  S.fill(0.0);
  // cross correlation between cartesian and polar measurements
  MatrixXd Tc = MatrixXd(n_x_, 2);
  Tc.fill(0.0);

  //calculate cross correlation matrix
  for(int i=0; i<(2*n_aug_ + 1); i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = tools_.phi_range(x_diff(3));

    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = tools_.phi_range(z_diff(1));

    S += weights_(i)*(z_diff*z_diff.transpose());
    Tc += weights_(i)*(x_diff*z_diff.transpose());
  }
  S += L_;

  MatrixXd K = Tc*S.inverse();

  VectorXd zdiff = z - z_pred;
  x_ = x_ + K*(zdiff);
  P_ = P_ + K*S*K.transpose();

  // calculate Laser NIS
  double nisl = zdiff.transpose()*S.inverse()*zdiff;
  nis_sample_laser_.push_back(nisl);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  VectorXd z = meas_package.raw_measurements_.head(3);

  MatrixXd Zsig = MatrixXd(3, 2*n_aug_+1);
  for(int i=0;i<(2*n_aug_+1);i++) {
    Zsig.col(i) = tools_.convert_to_polar(Xsig_pred_.col(i));
  }

  // find weighted mean of all polar converted sigma points
  VectorXd z_pred = weights_.transpose()*Zsig.transpose();
  z_pred(1) = tools_.phi_range(z_pred(1));
  VectorXd compare = tools_.convert_to_polar(x_);

  int n_z = 3;
  // calculate S matrix
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  // cross correlation between cartesian and polar measurements
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for(int i=0;i<(2*n_aug_+1);i++) {
      VectorXd z_diff = Zsig.col(i) - z_pred;
      z_diff(1) = tools_.phi_range(z_diff(1));
      S += weights_(i)*(z_diff*z_diff.transpose());

      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      x_diff(3) = tools_.phi_range(x_diff(3)); 
      Tc += weights_(i)*(x_diff*z_diff.transpose());
  }
  S += R_; // additive noise

  // find kalman gain
  MatrixXd K = Tc*S.inverse();

  // udpate state and covariances using the difference in predicted and observed
  VectorXd zdiff = z - z_pred;
  zdiff(1) = tools_.phi_range(zdiff(1));
  x_ = x_ + K*(zdiff);
  P_ = P_ - K*S*K.transpose();

  // calculate radar NIS
  double nisr = zdiff.transpose()*S.inverse()*zdiff;
  nis_sample_radar_.push_back(nisr);
}

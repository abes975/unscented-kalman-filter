#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <stdexcept>
#include <cfloat>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

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
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

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

   // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);;
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  weights_.tail(2*n_aug_).fill( 0.5 / (lambda_ + n_aug_) );

  // state vector
  x_ = VectorXd(n_x_);

  // state covariance matrix P
  // (Used Identity matrix as initialization)
  P_ = MatrixXd::Identity(n_x_, n_x_);

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // time when the state is true, in us
  time_us_ = 0;

  // the current NIS for radar
  NIS_radar_ = 0;

  // the current NIS for laser
  NIS_laser_ = 0;

  is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if (!is_initialized_) {
    // first measurement
    cout << "UKF: " << endl;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      x_ << Tools::PolarToCartesian(meas_package.raw_measurements_),
            0, 0;
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_,
            0, 0, 0;
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }


  // compute the time elapsed between the current and previous measurements
  //dt - expressed in seconds
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  if (use_radar_ &&
      meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR) {
        UpdateRadar(meas_package);
  } else if (use_laser_ &&
      meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER) {
        UpdateLidar(meas_package);
  } else {
    throw runtime_error("Unknown measurement type specified, "
        "nor lidar or radar not used");
  }

  cout << "x_ = " << endl << x_ << endl;
  cout << "P_ = " << endl << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
// Probably better separate this funcion in different function...
// Now go like this...but think about it.
void UKF::Prediction(double delta_t) {
  // The following code is extracted from the lessons...
  //augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  //augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  //sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //ugmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //set augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //calculate square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i=0; i<n_aug_; i++) {
    Xsig_aug.col(1+i) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(1+i+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  // From here predict sigma points
  for (int i=0; i<2 * n_aug_ + 1; i++) {
    double px = Xsig_aug(0,i);
    double py = Xsig_aug(1,i);
    double v  = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawd = Xsig_aug(6,i);
    double px_pred, py_pred;

    //avoid division by zero
    if (fabs(yawd) > FLT_EPSILON) {
      px_pred = px + v/yawd * ( sin (yaw + yawd * delta_t) - sin(yaw));
      py_pred = py + v/yawd * ( -cos(yaw + yawd * delta_t) + cos(yaw));
    } else {
      px_pred = px + v * delta_t * cos(yaw);
      py_pred = py + v * delta_t * sin(yaw);
    }

    // add noise
    double delta_square =  delta_t * delta_t;
    px_pred += 0.5 * nu_a * delta_square * cos(yaw);
    py_pred += 0.5 * nu_a * delta_square * sin(yaw);

    double v_pred = v + delta_t * nu_a;
    double yaw_pred = yaw + yawd * delta_t + 0.5 * delta_square * nu_yawd;
    double yawd_pred = yawd + delta_t * nu_yawd;

    // write predicted sigma points into right column
    Xsig_pred_(0,i) = px_pred;
    Xsig_pred_(1,i) = py_pred;
    Xsig_pred_(2,i) = v_pred;
    Xsig_pred_(3,i) = yaw_pred;
    Xsig_pred_(4,i) = yawd_pred;
  }


  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = Tools::normalize_angle(x_diff(3));
    P_ += weights_(i) * x_diff * x_diff.transpose() ;
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
  int n_z = 2;

  // transform sigma points into measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  for (int i=0; i<2 * n_aug_ + 1; i++) {
    Zsig(0,i) = Xsig_pred_(0,i);
    Zsig(1,i) = Xsig_pred_(1,i);
  }

  // calculate mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  // calculate measurement covariance matrix S and cross correlation matrix
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = Tools::normalize_angle(x_diff(3));

    VectorXd z_diff = Zsig.col(i) - z_pred;

    S += weights_(i) * z_diff * z_diff.transpose();
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  S(0,0) += std_laspx_ * std_laspx_;
  S(1,1) += std_laspy_ * std_laspy_;

  /* Update state with lidar measurement */
  // calculate Kalman gain K;
  MatrixXd S_inv = S.inverse();
  MatrixXd K = Tc * S_inv;

  // update state mean and covariance matrix
  VectorXd z_diff_mean =  meas_package.raw_measurements_ - z_pred;
  x_ += K * z_diff_mean;
  P_ -= K * S * K.transpose();

  // calculate the NIS
  NIS_laser_ = z_diff_mean.transpose() * S_inv * z_diff_mean;
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

  int n_z = 3;

  // transform sigma points into measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  for (int i=0; i<2 * n_aug_ + 1; i++) {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double yawd = Xsig_pred_(4,i);

    Zsig(0,i) = sqrt(px * px + py * py);
    Zsig(1,i) = atan2(py, px);
    Zsig(2,i) = v * (px * cos(yaw) + py * sin(yaw)) / Zsig(0,i);
  }

  // calculate mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  // calculate measurement covariance matrix S
  // and cross correlation matrix
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = Tools::normalize_angle(x_diff(3));

    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = Tools::normalize_angle(z_diff(1));

    S += weights_(i) * z_diff * z_diff.transpose();
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  S(0,0) += std_radr_ * std_radr_;
  S(1,1) += std_radphi_ * std_radphi_;
  S(2,2) += std_radrd_ * std_radrd_;

  MatrixXd S_inv = S.inverse();
  MatrixXd K = Tc * S_inv;
  VectorXd z_diff_mean =  meas_package.raw_measurements_ - z_pred;
  x_ += K * z_diff_mean;
  P_ -= K * S * K.transpose();

  NIS_radar_ = z_diff_mean.transpose() * S_inv * z_diff_mean;
}

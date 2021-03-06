#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  // state transition matrix F
  F_ = MatrixXd(4, 4);
  
  // measurement matrix H
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  
  //state covariance matrix P
  P_ = MatrixXd(4, 4);
  
  // process covariance matrix Q
  Q_ = MatrixXd(4,4);
  Q_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    
    // state covariance
    P_ << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1000, 0,
          0, 0, 0, 1000;
    ekf_.P_ = P_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      
      // Polar Conversion
      ekf_.x_ = tools.ConvertPolarToCartesian(measurement_pack.raw_measurements_);
      cout << "    EKF: Init w RADAR" << endl;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0],
                  measurement_pack.raw_measurements_[1],
                  0,0;
      cout << "    EKF: Init w LiDAR" << endl;
    }
    // Save initalization timestamp
    previous_timestamp_ = measurement_pack.timestamp_;
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   *****************************************************************************
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  
  long long timestamp = measurement_pack.timestamp_;
  float dt = (timestamp - previous_timestamp_)/ 1000000.0; // to seconds
  
  // Update F matrix
  F_ << 1, 0, dt, 0,
        0, 1, 0 ,dt,
        0, 0, 1, 0,
        0, 0, 0, 1;
  ekf_.F_ = F_;
  
  // Update process Covariance
  Q_  <<  noise_ax/4*powf(dt,4), 0, noise_ax/2*powf(dt,3), 0,
          0, noise_ay/4*powf(dt,4), 0, noise_ay/2*powf(dt,3),
          noise_ax/2*powf(dt,3), 0, noise_ax*powf(dt,2), 0,
          0, noise_ay/2*powf(dt,3), 0, noise_ay*powf(dt,2);
  ekf_.Q_ = Q_;
  
  // Predict new State + covariance (P)
  ekf_.Predict();
  
  // Update Time Stamp
  previous_timestamp_ = timestamp;
  /*****************************************************************************
   *  Update
   *****************************************************************************
   * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    cout << "    Measurement: RADAR" << endl;
    // Radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    cout << "    Measurement: LiDAR" << endl;
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

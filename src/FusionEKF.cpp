#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

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

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_<<1,0,0,0,
		    0,1,0,0;

  //set the acceleration noise components
  noise_ax = 9.0;
  noise_ay = 9.0;


  P_ = MatrixXd(4,4);
  P_<< 1,0,0,0,
	   0,1,0,0,
	   0,0,1000,0,
	   0,0,0,1000;

  F_ = MatrixXd(4,4);
  F_<< 1,0,1,0,
	   0,1,0,1,
	   0,0,1,0,
	   0,0,0,1;

  Q_ = MatrixXd(4,4);
  Q_ << 0,0,0,0,
		0,0,0,0,
		0,0,0,0,
		0,0,0,0;

  x_ = VectorXd(4);
  x_<< 1,1,1,1;

  Hj_<< 0,0,0,0,
		0,0,0,0,
		0,0,0,0;

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

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

    ekf_.Init(x_,P_,F_,Hj_,R_radar_,Q_);
    float rho,theta,rho_dot;
    rho = measurement_pack.raw_measurements_[0];
    theta = measurement_pack.raw_measurements_[1];
    rho_dot = measurement_pack.raw_measurements_[2];

    x_(0) = rho*cos(theta);
    x_(1) = rho*sin(theta);
    x_(2) = rho_dot*cos(theta);
    x_(3) = rho_dot*sin(theta);

    ekf_.x_ = x_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_, Hj_);
    Hj_ = ekf_.H_;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */

    ekf_.Init(x_,P_,F_,H_laser_,R_laser_,Q_);
    ekf_.x_<<measurement_pack.raw_measurements_[0],measurement_pack.raw_measurements_[1],
    		0,0;


    }

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }


  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
/*
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
*/
  float dt_2 = pow(dt,2);
  float dt_3 = pow(dt,3);
  float dt_4 = pow(dt,4);

  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  /*
  ekf_.F_ << 1, 0, dt, 0,
		     0, 1, 0, dt,
			 0, 0, 1, 0,
			 0, 0, 0, 1;
*/

  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
    		  0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
    		  dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
    		  0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	ekf_.H_ = tools.CalculateJacobian(ekf_.x_, Hj_);
	Hj_ = ekf_.H_;
	ekf_.R_ = R_radar_;

	ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
	ekf_.H_ = H_laser_;
	ekf_.R_ = R_laser_;
	ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

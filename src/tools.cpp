#include <iostream>
#include <cmath>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE here.
  */
  
  VectorXd rmse(4), e, e2;
  rmse << 0,0,0,0;
  
    // check the validity of the following inputs:
  if (estimations.size() == 0){
    return rmse;
  }
  
    //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size() ){
    return -rmse;
  }
  
    //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
    e = ground_truth[i] - estimations[i];
    e2 = e.array()*e.array();
    rmse = rmse + e2;
  }
  
  rmse = rmse / estimations.size();
  
    //calculate the squared root
  rmse = rmse.array().sqrt();
  
    //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian here.
  */
  
  MatrixXd Hj(3,4);
    //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  // Radial Distance
  float mag = sqrt(px*px + py*py);
  
  //check division by zero
  if (mag == 0){
    return Hj;
  }
  
  //compute the Jacobian matrix
  float h11 = px / mag;
  float h12 = py / mag;
  float h21 = -py / powf(mag,2);
  float h22 = px / powf(mag,2);
  
  // position
  Hj(0,0) = h11;
  Hj(0,1) = h12;
  Hj(0,2) = 0;
  Hj(0,3) = 0;
  
  Hj(1,0) = h21;
  Hj(1,1) = h22;
  Hj(1,2) = 0;
  Hj(1,3) = 0;
  
  // velocity
  float v1 = py*(vx*py - vy*px)/powf(mag,3);
  float v2 = -px*(vy*px - vx*py)/powf(mag,3);
  float v3 = px/mag;
  float v4 = py/mag;
  
  Hj(2,0) = v1;
  Hj(2,1) = v2;
  Hj(2,2) = v3;
  Hj(2,3) = v4;
  
  return Hj;
}

VectorXd Tools::ConvertCartesianToPolar(const VectorXd &x_state) {
  // Calcualte the Polar representation of a cartesian state.
  VectorXd polar(3);
  
  // Input size checking
  if (x_state.size() != 4){
    return polar;
  }
  
  // State Elements
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  float p0 = sqrt(powf(px,2) + powf(py,2));
  float p1 = atan2(py, px);
  float p2 = (px*vx + py*vy)/p0;
  
  polar << p0, p1, p2;
  return polar;
}

VectorXd Tools::ConvertPolarToCartesian(const VectorXd &x_polar){
  // Convert a polar measurement into a cartesian state.
  VectorXd x_state(4);
  float x_0, y_0, v_x, v_y;
  
  float p = x_polar(0);
  float theta = x_polar(1);
  float p_dot = x_polar(2);
  
  cout << x_polar[0] << x_polar[1] << endl;
  
  if (p != 0){
    x_0 = p*cos(theta);
    y_0 = p*sin(theta);
    v_x = p_dot*cos(theta);
    v_y = p_dot*sin(theta);
  }
  else {
    x_0 = 0;
    y_0 = 0;
    v_x = 0;
    v_y = 0;
  }
  std::cout << x_0 << y_0 << endl;
  x_state << x_0,y_0,v_x,v_y;
  return x_state;
}


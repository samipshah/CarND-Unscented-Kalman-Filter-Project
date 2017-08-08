#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  VectorXd convert_to_polar(const VectorXd& cartesian);

  VectorXd convert_to_cartesian(const VectorXd& polar);

  double phi_range(double phi);
  
  double percentile_check(const vector<double> &dist, double percentile_value);
};

#endif /* TOOLS_H_ */
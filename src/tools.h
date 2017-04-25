#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

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
  static Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

    /*
    * A helper method to convert polar coordinates to cartesian coordinates.
    * Used when the source of measurement is RADAR
    */
  static Eigen::VectorXd PolarToCartesian(const Eigen::VectorXd& polar);

    /*
    * A helper method to convert cartesian coordinates to polar.
    * Used when the source of measurement is  RADAR in UpdateEKF
    */
  static Eigen::VectorXd CartesianToPolar(const Eigen::VectorXd& cartesian);

    /* Used in oder to assure that the angle is in the range -PI -:-  PI */
  static double normalize_angle(double to_normalize);
};

#endif /* TOOLS_H_ */

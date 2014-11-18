/**
 * @file VectorAverage.h
 * @author Scott Rasmussen (scott.rasmussen@zaita.com)
 * @github https://github.com/Zaita
 * @date 28/10/2014
 * @section LICENSE
 *
 * Copyright NIWA Science �2014 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * vector average penalty = a penalty on the mean of a vector being
 * difference from some user defined value (i) square(mean(vector)-k).
 * Optionally in log space (ii) square(mean(log(vector))-k) or (iii)
 * square(log(mean(vector))/k). Has an arbitrary multiplier (default=1) as well.
 */
#ifndef PENALTIES_VECTORAVERAGE_H_
#define PENALTIES_VECTORAVERAGE_H_

// headers
#include "Penalties/Penalty.h"

// namespaces
namespace isam {
namespace penalties {

/**
 * Class definition
 */
class VectorAverage : public isam::Penalty {
public:
  // methods
  VectorAverage();
  virtual                     ~VectorAverage() = default;
  Double                      GetScore() override final;

protected:
  // methods
  void                        DoValidate() override final;
  void                        DoBuild()  override final;

private:
  // members
  string                      parameter_ = "";
  string                      method_ = "";
  map<unsigned, Double>*      estimable_map_ = 0;
  vector<Double>*             estimable_vector_ = 0;
  Double                      multiplier_ = 0.0;
  Double                      k_;
};

} /* namespace penalties */
} /* namespace isam */

#endif /* VECTORAVERAGE_H_ */
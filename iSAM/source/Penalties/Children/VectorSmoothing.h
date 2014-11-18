/**
 * @file VectorSmoothing.h
 * @author Scott Rasmussen (scott.rasmussen@zaita.com)
 * @github https://github.com/Zaita
 * @date 28/10/2014
 * @section LICENSE
 *
 * Copyright NIWA Science �2014 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * vector smoothing = a penalty on the between element differences of a
 * vector of values (or map). Optionally in log space,. Will also work
 * on allvalues ogives. Equal to sum of squares of ith differences,
 * where i is a user defined option that defaults to 3. Has an arbitrary
 * multiplier (default=1) as well.
 */
#ifndef PENALTIES_VECTORSMOOTHING_H_
#define PENALTIES_VECTORSMOOTHING_H_

// headers
#include "Estimates/Estimate.h"
#include "Penalties/Penalty.h"

// namespaces
namespace isam {
namespace penalties {

/**
 * Class definition
 */
class VectorSmoothing : public isam::Penalty {
public:
  // methods
  VectorSmoothing();
  virtual                     ~VectorSmoothing() = default;
  Double                      GetScore() override final;

protected:
  // methods
  void                        DoValidate() override final;
  void                        DoBuild() override final;

private:
  // members
  string                      parameter_ = "";
  map<unsigned, Double>*      estimable_map_ = 0;
  vector<Double>*             estimable_vector_ = 0;
  bool                        log_scale_ = false;
  Double                      multiplier_ = 0.0;
  unsigned                    step_size_ = 0;
};

} /* namespace penalties */
} /* namespace isam */

#endif /* VECTORSMOOTHING_H_ */
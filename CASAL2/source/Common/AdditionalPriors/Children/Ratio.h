/**
 * @file Ratio.h
 * @author C.Marsh
 * @github https://github.com/
 * @date 28/06/2017
 * @section LICENSE
 *
 * Copyright NIWA Science ©2017 - www.niwa.co.nz
 *
 * @section An additional prior that is analagous to ratio Q's in the old Casal. Penalises for when ratios are drastically different
 *
 */
#ifndef ADDITIONALPRIORS_RATIO_H_
#define ADDITIONALPRIORS_RATIO_H_

// headers
#include "Common/AdditionalPriors/AdditionalPrior.h"


// namespaces
namespace niwa {
namespace additionalpriors {

// classes
class Ratio : public AdditionalPrior {
public:
  // methods
  Ratio(Model* model);
  virtual                     ~Ratio() = default;
  void                        DoValidate() override final;
  void                        DoBuild() override final;
  Double                      GetScore() override final;

protected:
  // members
  Double* 										addressable_ = nullptr;

  Double                      mu_ = 0.0;
  Double                      sigma_ = 0.0;
  Double                      a_ = 0.0;
  Double                      b_ = 0.0;
  Double                      v_ = 0.0;
  Double                      t_ = 0.0;
  Double                      m_ = 0.0;
  Double                      n_ = 0.0;
};

} /* namespace additionalpriors */
} /* namespace niwa */

#endif /* ADDITIONALPRIORS_RATIO_H_ */

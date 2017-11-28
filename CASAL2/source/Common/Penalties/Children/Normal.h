/**
 * @file Process.h
 * @author Scott Rasmussen (scott.rasmussen@zaita.com)
 * @github https://github.com/Zaita
 * @date 28/10/2014
 * @section LICENSE
 *
 * Copyright NIWA Science ©2014 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * Process penalties (i.e., special flagged penalties like the catch limit penalty.
 * Has an arbitrary multiplier (default=1) as well.
 */
#ifndef PENALTIES_NORMAL_H_
#define PENALTIES_NORMAL_H_

// headers
#include "Common/Penalties/Penalty.h"

// namespaces
namespace niwa {
namespace penalties {

/**
 * Class definition
 */
class Normal : public niwa::Penalty {
public:
  // methods
  Normal(Model* model);
  virtual                     ~Normal() = default;
  void                        Trigger(const string& source_label, Double value, Double std_dev, Double bias_adjustment);
  Double                      GetScore() override final { return 0.0; }
protected:
  // methods
  void                        DoValidate() override final { };
  void                        DoBuild() override final { };

private:
  // members

};
} /* namespace penalties */
} /* namespace niwa */

#endif /* PENALTIES_NORMAL_H_ */

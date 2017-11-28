/**
 * @file Normal.cpp
 * @author C.Marsh
 * @github https://github.com/Zaita
 * @date 4/9/2017
 * @section LICENSE
 *
 * Copyright NIWA Science ©2017 - www.niwa.co.nz
 *
 * Applies a bias corrected standard normal contribution to the objective function.
 */

// headers
#include "Normal.h"

#include "Common/Model/Model.h"
#include "Common/Model/Managers.h"
#include "Common/Penalties/Manager.h"
#include "Common/Utilities/DoubleCompare.h"

// namespaces
namespace niwa {
namespace penalties {

/**
 * Default Constructor
 */
Normal::Normal(Model* model) : Penalty(model) {
  has_score_ = false;
}

/**
 * Trigger our penalty.
 * Basic value for the trigger will be: (value_1 - value_2)^2 * multiplier
 * logscale is: (log(value_1) - log(value_2))^2 * multiplier
 *
 * @param source_label The label for the source of the trigger
 * @param value_1 The first value to use in equation
 * @param value_2 The second valud to use in equatin
 */
void Normal::Trigger(const string& source_label, Double value, Double std_dev, Double bias_adjustment) {
  Double contribution;
  if (value == 0.0) {
    contribution = 0.0;
  } else
    contribution = 0.5 * (value * value / (std_dev * std_dev) + bias_adjustment * log(std_dev));
  string name  = label_ + "(" + source_label + ")";
  model_->managers().penalty()->FlagPenalty(name, contribution);
}

} /* namespace penalties */
} /* namespace niwa */

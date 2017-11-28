/**
 * @file Ratio.cpp
 * @author C.Marsh
 * @github https://github.com/
 * @date 28/06/2017
 * @section LICENSE
 *
 * Copyright NIWA Science ©2017 - www.niwa.co.nz
 *
 */

// headers
#include "Ratio.h"

#include "Common/Estimates/Manager.h"
#include "Common/Model/Model.h"
#include "Common/Model/Objects.h"

// namespaces
namespace niwa {
namespace additionalpriors {

/**
 * Default constructor
 *
 * Bind any parameters that are allowed to be loaded from the configuration files.
 * Set bounds on registered parameters
 * Register any parameters that can be an estimated or utilised in other run modes (e.g profiling, yields, projections etc)
 * Set some initial values
 *
 * Note: The constructor is parsed to generate Latex for the documentation.
 */
Ratio::Ratio(Model* model) : AdditionalPrior(model) {

}

/**
 * Populate any parameters,
 * Validate values are within expected ranges when we cannot use bind<>() overloads
 *
 * Note: all parameters are populated from configuration files
 */
void Ratio::DoValidate() {
  LOG_TRACE();

}

/**
 * Build
 */
void Ratio::DoBuild() {
  LOG_TRACE();

}

/**
 * Return the score for something
 */
Double Ratio::GetScore() {
  LOG_TRACE();

 return 0.0;
}

} /* namespace additionalpriors */
} /* namespace niwa */

/**
 * @file Maturation.h
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @date 11/09/2013
 * @section LICENSE
 *
 * Copyright NIWA Science �2013 - www.niwa.co.nz
 *
 * @section DESCRIPTION
 *
 * << Add Description >>
 */
#ifndef PROCESSES_MATURATION_H_
#define PROCESSES_MATURATION_H_

// headers
#include "Model/Model.h"
#include "Processes/Process.h"
#include "Partition/Accessors/Age/Categories.h"
#include "Selectivities/Selectivity.h"

// namespaces
namespace niwa {
namespace processes {

namespace accessor = niwa::partition::accessors;

/**
 *
 */
class Maturation : public niwa::Process {
public:
  // methods
  Maturation();
  virtual                     ~Maturation() = default;
  void                        DoValidate() override final;
  void                        DoBuild() override final;
  void                        DoReset() override final { };
  void                        DoExecute() override final;

protected:
  // Members
  ModelPtr                    model_;
  vector<string>              from_category_names_;
  vector<string>              to_category_names_;
  vector<Double>              rates_;
  vector<unsigned>            years_;
  vector<string>              selectivity_names_;
  accessor::Categories        from_partition_;
  accessor::Categories        to_partition_;
  vector<SelectivityPtr>      selectivities_;
  map<unsigned, Double>       rates_by_years_;
};

} /* namespace processes */
} /* namespace niwa */
#endif /* MATURATION_H_ */

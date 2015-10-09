/*
 * Derived.cpp
 *
 *  Created on: 2/09/2014
 *      Author: Admin
 */

#include "Derived.h"

#include <algorithm>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/algorithm/string/split.hpp>

#include "Categories/Categories.h"
#include "DerivedQuantities/Manager.h"
#include "Partition/Accessors/Categories.h"
#include "Processes/Manager.h"
#include "TimeSteps/Factory.h"
#include "Model/Model.h"

namespace niwa {
namespace initialisationphases {

namespace cached   = partition::accessors::cached;
namespace accessor = partition::accessors;

Derived::Derived() {
    parameters_.Bind<string>(PARAM_INSERT_PROCESSES, &insert_processes_, "The processes to insert in to target time steps", "", true);
    parameters_.Bind<string>(PARAM_EXCLUDE_PROCESSES, &exclude_processes_, "The processes to exclude from all time steps", "", true);
}

/*
 * Validate the format of insert_processes if any are given
 */
void Derived::DoValidate() {
  for (string insert : insert_processes_) {
    vector<string> pieces;
    boost::split(pieces, insert, boost::is_any_of("()="), boost::token_compress_on);
    if (pieces.size() != 2 && pieces.size() != 3)
      LOG_ERROR_P(PARAM_INSERT_PROCESSES) << " value " << insert << " does not match the format time_step(process)=new_process = " << pieces.size();
  }
}

/*
 * Build any runtime relationships needed for execution
 */
void Derived::DoBuild() {
  time_steps_ = timesteps::Manager::Instance().ordered_time_steps();

  // Set the default process labels for the time step for this phase
  for (TimeStepPtr time_step : time_steps_)
    time_step->SetInitialisationProcessLabels(label_, time_step->process_labels());

  // handle any new processes we want to insert
  for (string insert : insert_processes_) {
    vector<string> pieces;
    boost::split(pieces, insert, boost::is_any_of("()="), boost::token_compress_on);

    string target_process   = pieces.size() == 3 ? pieces[1] : "";
    string new_process      = pieces.size() == 3 ? pieces[2] : pieces[1];

    TimeStepPtr time_step = timesteps::Manager::Instance().GetTimeStep(pieces[0]);
    vector<string> process_labels = time_step->initialisation_process_labels(label_);

    if (target_process == "") {
      process_labels.insert(process_labels.begin(), new_process);
    } else {
      vector<string>::iterator iter = std::find(process_labels.begin(), process_labels.end(), target_process);
      if (iter == process_labels.end())
        LOG_ERROR_P(PARAM_INSERT_PROCESSES) << " process " << target_process << " does not exist in time step " << time_step->label();
      process_labels.insert(iter, new_process);
    }

    time_step->SetInitialisationProcessLabels(label_, process_labels);
  }

  // handle the excludes we've specified
  for (string exclude : exclude_processes_) {
    unsigned count = 0;
    for (TimeStepPtr time_step : time_steps_) {
      vector<string> process_labels = time_step->initialisation_process_labels(label_);
      unsigned size_before = process_labels.size();
      process_labels.erase(std::remove_if(process_labels.begin(), process_labels.end(), [exclude](string& ex) { return exclude == ex; }), process_labels.end());
      unsigned diff = size_before - process_labels.size();

      time_step->SetInitialisationProcessLabels(label_, process_labels);
      count += diff;
    }

    if (count == 0)
      LOG_ERROR_P(PARAM_EXCLUDE_PROCESSES) << " process " << exclude << " does not exist in any time steps to be excluded. Please check your spelling";
  }

  // Build our partition
  vector<string> categories = Categories::Instance()->category_names();
  partition_.Init(categories);
  cached_partition_.Init(categories);

  // Build a model pointer so that we can access min and max age that all category ages must be within


}

/*
 * Execute our Derived Initialisation phase
 */
void Derived::Execute() {

  ModelPtr model = Model::Instance();
  unsigned year_range = model->max_age() - model->min_age() + 1;
  // Maybe in dobuild minus the following condition off max_age ((recruitment_time < ageing_time) ? 1 : 0)
  vector<string> categories = Categories::Instance()->category_names();

  timesteps::Manager& time_step_manager = timesteps::Manager::Instance();
  time_step_manager.ExecuteInitialisation(label_, year_range);

  cached_partition_.BuildCache();
  // a shortcut to avoid running the model over more years to get the plus group right
  // calculate the annual change c for each element of the plus group
  if (model->age_plus()) {
    // Run the model for an extra year
    timesteps::Manager& time_step_manager = timesteps::Manager::Instance();
    time_step_manager.ExecuteInitialisation(label_, 1);
    Double c;
    auto cached_category = cached_partition_.begin();
    auto category = partition_.begin();

    for (; category != partition_.end(); ++cached_category, ++category) {

      // find the element in the data container (numbers at age) that contains the plus group for each category
      // We are assuming it is the last element I need to check with Scott with this.
      unsigned plus_index = (*category)->data_.size() - 1;

      LOG_FINEST() << "Did any numbers populate the plus group " << cached_category->data_[plus_index] << " " << (*category)->data_[plus_index];
      //if (the plus group has been populated)
      if (cached_category->data_[plus_index] > 0) {
        c = (*category)->data_[plus_index] / cached_category->data_[plus_index] - 1;
        LOG_FINEST() << "The value of c = " << c;
        if (c > 0.9) {
          c = 0.9;
        } else if (c < 0.0)
          c = 0.0;
        // reset the partition back to the original Cached partition
        (*category)->data_ = cached_category->data_;
        // now add the approximated plus group to the partition
        (*category)->data_[plus_index] *= 1 / (1 - c);
        LOG_FINEST() << "Adjustment based an approximation for the plus group" <<  (*category)->data_[plus_index];
      } else {
        // reset the partition back to the original Cached partition
        (*category)->data_ = cached_category->data_;
      }
    }
  }
  Double max_rel_diff = 1e18;
  vector<Double> plus_group(1,categories.size()), old_plus_group(1,categories.size());

  auto category = partition_.begin();
  unsigned iter = 0;
  for (; category != partition_.end(); ++category, ++iter)
    old_plus_group[iter] = (*category)->data_[(*category)->data_.size() - 1];
  LOG_FINEST() << " max diff " << max_rel_diff;
  while (max_rel_diff > 0.0005){
    timesteps::Manager& time_step_manager = timesteps::Manager::Instance();
    time_step_manager.ExecuteInitialisation(label_, 1);
    max_rel_diff = 0;
    auto category = partition_.begin();
    unsigned iter = 0;
    for (; category != partition_.end(); ++category, ++iter) {
      plus_group[iter] = (*category)->data_[(*category)->data_.size() - 1];
      LOG_FINEST() << " plus group " << plus_group[iter] << " Old plus group " << old_plus_group[iter];
      if (old_plus_group[iter] != 0){
        if (fabs((plus_group[iter]-old_plus_group[iter])/old_plus_group[iter]) > max_rel_diff)
          max_rel_diff = fabs((plus_group[iter] - old_plus_group[iter]) / old_plus_group[iter]);
        LOG_FINEST() << " max diff " << max_rel_diff;
      }
    }
    old_plus_group = plus_group;
  }
}


} /* namespace initialisationphases */
} /* namespace niwa */

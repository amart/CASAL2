/**
 * @file Categories.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 11/12/2012
 * @section LICENSE
 *
 * Copyright NIWA Science �2012 - www.niwa.co.nz
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */

// Headers
#include "Categories.h"

#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/algorithm/string/split.hpp>

#include "AgeLengths/Manager.h"
#include "Model/Model.h"
#include "Utilities/Logging/Logging.h"
#include "Utilities/String.h"
#include "Utilities/To.h"

// Namespaces
namespace niwa {

/**
 * Default constructor
 */
Categories::Categories() {
  LOG_TRACE();
  parameters_.Bind<string>(PARAM_FORMAT, &format_, "The format that the category names should adhere too", "");
  parameters_.Bind<string>(PARAM_NAMES, &names_, "The names of the categories to be used in the model", "");
  parameters_.Bind<string>(PARAM_YEARS, &years_, "The years that individual categories will be active for. This overrides the model values", "", true);
  parameters_.Bind<string>(PARAM_AGES, &ages_, "The ages that individual categories support. This overrides the model values", "", true);
  parameters_.Bind<string>(PARAM_AGE_SIZES, &age_size_labels_, R"(The labels of age\_size objects that are assigned to categories)", "", true);
}

/**
 * Our singleton accessor method
 *
 * @return singleton shared ptr
 */
shared_ptr<Categories> Categories::Instance() {
  static CategoriesPtr categories = CategoriesPtr(new Categories());
  return categories;
}

/**
 * Validate our Categories
 */
void Categories::Validate() {
  // Check that we actually had a categories block
  if (block_type_ == "")
    LOG_ERROR("The @categories block is missing from the configuration file. This block is required");

  parameters_.Populate();

  /**
   * Parameter: Names
   *
   * 1. Verify we have some names present
   * 2. Verify the names match the format string properly (right amount of segments)
   * 3. Create a new CategoryInfo object for each category and store it
   */

  // Parameter: Names
  ModelPtr model = Model::Instance();

  // build the default years
  vector<unsigned> default_years;
  for (unsigned i = model->start_year(); i <= model->final_year(); ++i)
    default_years.push_back(i);

  // get the age sizes
  if (age_size_labels_.size() > 0 && age_size_labels_.size() != names_.size())
    LOG_ERROR(parameters_.location(PARAM_AGE_SIZES) << " number defined (" << age_size_labels_.size() << ") must be the same as the number " <<
        " of categories defined (" << names_.size() << ")");

  // build our categories vector
  for (unsigned i = 0; i < names_.size(); ++i) {
    if (age_size_labels_.size() > i)
      category_age_size_labels_[names_[i]] = age_size_labels_[i];

    // TODO: Verify the name matches the format string properly
    // TODO: Expand the names

    // Create a new CategoryInfo object
    CategoryInfo new_category_info;
    new_category_info.name_     = names_[i];
    new_category_info.min_age_  = model->min_age();
    new_category_info.max_age_  = model->max_age();
    new_category_info.years_    = default_years;
    categories_[names_[i]] = new_category_info;

    category_names_.push_back(names_[i]);
  }
}

/**
 * Build runtime relationships between the categories
 * and other objects in the system
 */
void Categories::Build() {
  agelengths::Manager& age_sizes_manager = agelengths::Manager::Instance();

  auto iter = category_age_size_labels_.begin();
  for (; iter != category_age_size_labels_.end(); ++iter) {
    AgeLengthPtr age_size = age_sizes_manager.GetAgeLength(iter->second);
    if (!age_size)
      LOG_ERROR(parameters_.location(PARAM_AGE_SIZES) << "(" << iter->second << ") could not be found. Have you defined it?");

    categories_[iter->first].age_length_ = age_size;
  }
}

/**
 * This method will expand short-hand syntax in to a vector with all of the elements
 * separated.
 *
 * e.g male.immature sex=female == male.immature female.immature female.mature == 3 elements
 *
 * @param category_labels A vector of category definitions to iterate over and expand
 * @param source_parameter The parameter object which holds configuration file details for error reporting
 * @return a singular vector with each category as it's own element
 */
vector<string> Categories::ExpandLabels(const vector<string> &category_labels, const ParameterPtr source_parameter) {
  vector<string> result;

  vector<string> temp;
  for (const string& label : category_labels) {
    temp = GetCategoryLabelsV(label, source_parameter);
    result.insert(result.end(), temp.begin(), temp.end());
  }

  return result;
}

/**
 *
 */
vector<string> Categories::GetCategoryLabelsV(const string& lookup_string, const ParameterPtr source_parameter) {
  string temp = GetCategoryLabels(lookup_string, source_parameter);

  vector<string> result;
  boost::split(result, temp, boost::is_any_of(" "));
  return result;
}

/**
 * This method will take a lookup string and parse it looking for our short-hand
 * syntax. The source_parameter parameter is passed in to allow us to print
 * proper error messages without having to pass back any enum types or something
 * silly like that.
 *
 * @param lookup_string The category definition to parse, short-hand or not
 * @param source_parameter Source parameter object defined in configuration file
 * @return String containing the new lookup string once it's been parsed
 */
string Categories::GetCategoryLabels(const string& lookup_string, const ParameterPtr source_parameter) {
  /**
   * if we're asking for all categories then get a list of them all joined
   * by the + symbol
   */
  if (lookup_string == "*+") {
    string all = category_names_[0];
    for (unsigned i = 1; i < category_names_.size(); ++i)
      all += "+" + category_names_[i];

    return all;

  } else if (lookup_string == "*") {
    string all = category_names_[0];
    for (unsigned i = 1; i < category_names_.size(); ++i)
      all += " " + category_names_[i];

    return all;
  }

  /**
   * Check if we're just asking for a single category
   * or collection without any short-hand syntax
   */
  if (lookup_string.find("=") == string::npos)
    return lookup_string;

  /**
   * From here forward we're working with short-hand syntax.
   * Either
   * <format_chunk>=<lookup_string>[+]
   * or
   * format=<lookup_string_0>.<lookup_string_1>.etc[+]
   */
  string result = "";

  vector<string> pieces;
  boost::split(pieces, lookup_string, boost::is_any_of("="));

  if (pieces.size() != 2) {
    LOG_ERROR(source_parameter->location() << " short-hand category string (" << lookup_string
        << ") is not in the proper format (e.g <format_chunk>=<lookup_chunk>)");
  }

  boost::replace_all(pieces[0], " ", "");
  boost::replace_all(pieces[1], " ", "");

  string format = pieces[0];
  string lookup = pieces[1];

  bool use_plus_join = false;
  if (lookup.find("+") != string::npos) {
    use_plus_join = true;
    boost::replace_all(lookup, "+", "");
  }

  // organise categories split into segments to match against format
  map<string, vector<string>> category_pieces;
  for (string category : category_names_) {
    vector<string> temp;
    boost::split(temp, category, boost::is_any_of("."));
    category_pieces[category] = temp;
  }

  vector<string> matched_categories = category_names_;
  if (utilities::ToLowercase(format) == PARAM_FORMAT) {
    /**
     * this is a full length lookup string. It must have the
     * same number of pieces as our actual format string
     * e.g
     * format sex.stage.tag[+]
     * lookup_string male.*.*
     */
    boost::split(pieces, format_, boost::is_any_of("."));
    unsigned format_pieces = pieces.size();
    boost::split(pieces, lookup, boost::is_any_of("."));

    if (pieces.size() != format_pieces) {
      LOG_ERROR(source_parameter->location() << " short-hand category string ( " << lookup_string
          << ") does not have the correct number of sections. Expected " << format_pieces << " but got " << pieces.size() <<
          ". Pieces are chunks of the string separated with a '.' character");
    }

    for (unsigned i = 0; i < pieces.size(); ++i) {
      if (pieces[i] == "*")
        continue;

      vector<string> comma_separated_pieces;
      boost::split(comma_separated_pieces, pieces[i], boost::is_any_of(","));
      matched_categories.erase(
          std::remove_if(matched_categories.begin(), matched_categories.end(),
          [&i, &category_pieces, &comma_separated_pieces](string& category) {

            if (std::find(comma_separated_pieces.begin(), comma_separated_pieces.end(), category_pieces[category][i]) == comma_separated_pieces.end())
              return true;

            return false;
          }),
          matched_categories.end()
      );
    }

    LOG_INFO("Full format parse of categories returned " << matched_categories.size() << " results");
  } else {
    /**
     * Here we have the shorter form of syntax
     * category_chunk=filter
     * e.g
     * sex=male
     */
    if (lookup.find(".") != string::npos) {
      LOG_ERROR(source_parameter->location() << " short-hand category string (" << lookup_string
          << ") is not in the correct format. The lookup component (" << lookup
          << ") cannot contain any '.' characters");
    }

    // Verify we've actually got a good part of the format here.
    boost::split(pieces, format_, boost::is_any_of("."));
    unsigned format_offset = pieces.size();
    for (unsigned i = 0; i < pieces.size(); ++i) {
      if (pieces[i] == format) {
        format_offset = i;
        break;
      }
    }
    if (format_offset == pieces.size()) {
      LOG_ERROR(source_parameter->location() << " short-hand category syntax (" << lookup_string
          << ") is using an invalid format chunk (" << format << ") for it's lookup. "
          << "Valid format chunks must be taken from the format (" << format_ << ")");
    }

    matched_categories.erase(
        std::remove_if(matched_categories.begin(), matched_categories.end(),
        [&format_offset, &source_parameter, &lookup, &lookup_string](string& category) {
          vector<string> chunks;
          boost::split(chunks, category, boost::is_any_of("."));
          if (chunks.size() <= format_offset) {
            LOG_ERROR(source_parameter->location() << " short-hand category syntax (" << lookup_string
                << ") could not be compared to category (" << category << ") because category was malformed");
          }

          vector<string> comma_separated_pieces;
          boost::split(comma_separated_pieces, lookup, boost::is_any_of(","));

          if (std::find(comma_separated_pieces.begin(), comma_separated_pieces.end(), chunks[format_offset]) == comma_separated_pieces.end())
            return true;
          return false;
        }),
        matched_categories.end()
    );

    LOG_INFO("Short format parse of categories returned " << matched_categories.size() << " results");
  }


  if (matched_categories.size() == 0) {
    LOG_ERROR(source_parameter->location() << " short-hand format string (" << lookup_string <<
        ") did not match any of the categories. Please check your string to ensure it's accurate");
  }

  result = matched_categories[0];
  for (unsigned i = 1; i < matched_categories.size(); ++i) {
    if (use_plus_join)
      result += "+" + matched_categories[i];
    else
      result += " " + matched_categories[i];
  }

  return result;
}

/**
 * Check if the parameter category label is a valid
 * category or not.
 */
bool Categories::IsValid(const string& label) const {
  return std::find(category_names_.begin(), category_names_.end(), label) != category_names_.end();
}

/**
 *
 */
bool Categories::IsCombinedLabels(const string& label) const {
  return (label.find_first_of("+") != string::npos);
}

/**
 *
 */
unsigned Categories::GetNumberOfCategoriesDefined(const string& label) const {
  vector<string> category_labels;
  boost::split(category_labels, label, boost::is_any_of("+"));

  return category_labels.size();
}

/**
 *
 */
unsigned Categories::min_age(const string& category_name) {
  if (categories_.find(category_name) == categories_.end())
    LOG_CODE_ERROR("Could not find category_name: " << category_name << " in the list of loaded categories");

  return categories_[category_name].min_age_;
}

/**
 *
 */
unsigned Categories::max_age(const string& category_name) {
  if (categories_.find(category_name) == categories_.end())
    LOG_CODE_ERROR("Could not find category_name: " << category_name << " in the list of loaded categories");

  return categories_[category_name].max_age_;
}

/**
 *
 */
vector<unsigned> Categories::years(const string& category_name) {
  if (categories_.find(category_name) == categories_.end())
    LOG_CODE_ERROR("Could not find category_name: " << category_name << " in the list of loaded categories");

  return categories_[category_name].years_;
}

/**
 *
 */
AgeLengthPtr Categories::age_length(const string& category_name) {
  if (categories_.find(category_name) == categories_.end())
    LOG_CODE_ERROR("Could not find category_name: " << category_name << " in the list of loaded categories");
  if (!categories_[category_name].age_length_)
    LOG_CODE_ERROR("The age size pointer was null for category " << category_name);

  return categories_[category_name].age_length_;
}

/**
 * This method will remove all of the information from our categories
 */
void Categories::RemoveAllObjects() {
  parameters_.Clear();
  format_ = "";
  names_.clear();
  category_names_.clear();
  categories_.clear();
  age_size_labels_.clear();
  category_age_size_labels_.clear();
}


} /* namespace niwa */

















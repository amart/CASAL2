/*
 * LogisticNormal.cpp
 *
 *  Created on: Oct 26, 2016
 *      Author: Zaita
 */

#include <Likelihoods/Children/LogisticNormal.h>

#include "Utilities/DoubleCompare.h"
#include "Utilities/Math.h"
#include "Utilities/RandomNumberGenerator.h"

namespace niwa {
namespace likelihoods {

using std::set;

namespace dc = niwa::utilities::doublecompare;
namespace math = niwa::utilities::math;

LogisticNormal::LogisticNormal() {
  parameters_.Bind<string>(PARAM_LABEL, &label_, "Label for Logisitic Normal Likelihood", "");
  parameters_.Bind<string>(PARAM_TYPE, &type_, "Type of likelihood", "");
  covariance_table_ = new parameters::Table(PARAM_COVARIANCE_MATRIX);

  parameters_.Bind<Double>(PARAM_RHO, &rho_, "The auto-correlation parameter $\rho$", "");
  parameters_.BindTable(PARAM_COVARIANCE_MATRIX, covariance_table_, "User defined Covariance matrix", "",false,true);
  parameters_.Bind<Double>(PARAM_SIGMA, &sigma_, "Sigma parameter in the likelihood", "");
  parameters_.Bind<bool>(PARAM_ARMA, &arma_, "Defines if two rho parameters supplied then covar is assumed to have the correlation matrix of an ARMA(1,1) process", "");
  parameters_.Bind<unsigned>(PARAM_BIN_LABELS, &bins_, "If no covariance matrix parameter then list a vector of bin labels that the covariance matrix will be built for, can be ages or lengths.", "",false);
  parameters_.Bind<bool>(PARAM_SEXED, &sexed_, "Will the observation be split by sex?", "",false);
  parameters_.Bind<bool>(PARAM_ROBUST, &robust_, "Robustification term for zero observations", "",false);
  parameters_.Bind<bool>(PARAM_SEPERATE_BY_SEX, &sep_by_sex_, "If data is sexed, should the covariance matrix be seperated by sex?", "",false);
  parameters_.Bind<bool>(PARAM_SEX_LAG, &sex_lag_, "if T and data are sexed, then the AR(n) correlation structure ignores sex and sets lag = |i-j|+1, where i & j index the age or length classes in the data.  Ignored if data are not sexed.", "",false);

  RegisterAsEstimable(PARAM_SIGMA, &sigma_);
  RegisterAsEstimable(PARAM_RHO, &rho_);


}

LogisticNormal::~LogisticNormal() {
  delete covariance_table_;
}

void LogisticNormal::DoValidate() {
  LOG_TRACE();
  if(!arma_ & (rho_.size() > 1)) {
    if(rho_[1] <= -1 || rho_[1] >= (1- fabs(rho_[0])))
      LOG_ERROR_P(PARAM_RHO) << "Incorrect values for rho. Rho cannot be less than -1 or the first value cannot be greater than (1 - |first value|)";
  }

  if(rho_.size() > 2)
    LOG_ERROR_P(PARAM_RHO) << "Can only have one value or two";
  /*
   * Build our covariance matrix with user defined values
   */
  if (parameters_.GetTable(PARAM_COVARIANCE_MATRIX)->has_been_defined()) {
    LOG_FINEST() << "Converting user defined covariance matrix";
    vector<vector<string>>& covariance_values_data = covariance_table_->data();
    int rows = covariance_values_data.size();
    int cols = covariance_values_data[0].size();
    if(rows != cols) {
      LOG_ERROR_P(PARAM_COVARIANCE_MATRIX) << "Covariance matrix must be square (equal rows and columns). rows = " << rows << " cols = " << cols;
    }
    covariance_matrix_(rows,cols);
    // Covariance should be an age by age matrix or age x sex by age x sex matrix; if comp data given by sex.
    LOG_FINEST() << " has " << covariance_values_data.size() << " rows defined";
    unsigned col_iter = 0;
    for (vector<string>& covariance_values_data_line : covariance_values_data) {
      unsigned bin = 0;
      vector<Double> temp;
      if (!utilities::To<unsigned>(covariance_values_data_line[0], bin)) {
        LOG_ERROR_P(PARAM_COVARIANCE_MATRIX) << " value " << covariance_values_data_line[0] << " could not be converted in to an unsigned integer. It should be the year for this line";
      } else {
        for (unsigned i = 1; i < covariance_values_data_line.size(); ++i) {
          Double value = 0;
          if (!utilities::To<Double>(covariance_values_data_line[i], value)) {
            LOG_ERROR_P(PARAM_SCANNED) << " value (" << covariance_values_data_line[i] << ") could not be converted to a Double";
          }
          covariance_matrix_(i - 1, col_iter) = value;
        }
      }
      bins_.push_back(bin);
      ++col_iter;
    }
    n_bins_ = bins_.size();
    if (sexed_)
      unique_bins_ /= 2;
    else
      unique_bins_ = n_bins_;
  } else {
    LOG_FINEST() << "Calculating Covariance matrix with user specified parameters";
    unique_bins_ = bins_.size();
    if (sexed_)
      n_bins_ = unique_bins_ * 2;
    else
      n_bins_ = unique_bins_;

    //initialise covariance matrix
    LOG_FINEST() << "number of bins = " << n_bins_ << " number of unique bins = " << unique_bins_;
    // Create a covariance matrix with the user defined parameters
    calculate_covariance();
  }

  LOG_FINEST() << "Printing Covariance top left triangle matrix";
  for (unsigned k = 0; k < covariance_matrix_.size1(); ++k) {
    for (unsigned j = 0; j < covariance_matrix_.size2(); ++j) {
      LOG_FINEST() << "row = " << k << " col = " << j << " val = " << covariance_matrix_(j,k) << " ";
    }
  }
}
/**
 * Adjust the error value based on the process error
 *
 * @param process_error The observations process_error
 * @param error_value The observations error_value
 * @return An adjusted error value
 */
Double LogisticNormal::AdjustErrorValue(const Double process_error, const Double error_value) {
  if (process_error > 0.0 && error_value > 0.0)
    return (1.0/(1.0/error_value + 1.0/process_error));

  return error_value;
}

/**
 *
 */
void LogisticNormal::GetScores(map<unsigned, vector<observations::Comparison> >& comparisons) {
  LOG_TRACE();
  unsigned N_years = comparisons.size();
  // Convert from a comparsion structure to a matrix, this likelihood, initialise matricies
  ublas::matrix<Double> Observed(N_years,n_bins_,0.0);
  ublas::matrix<Double> Expected(N_years,n_bins_,0.0);
  ublas::matrix<Double> Kmat(n_bins_ - 1,n_bins_,0.0);
  ublas::matrix<Double> Vmat(n_bins_ - 1,n_bins_ - 1,0.0);
  ublas::matrix<Double> Vmat_inv(n_bins_ - 1,n_bins_ - 1,0.0);
  ublas::matrix<Double> ww(N_years,n_bins_ - 1,0.0);
  ublas::matrix<Double> log_obs(N_years,n_bins_,0.0);
  LOG_FINEST() << "Finished intialising temp matricies";
  Double Score = 0.0;
  vector<Double> N,weights;
  LOG_FINEST() << "Dim Observed rows = " << Observed.size1() << " cols = " << Observed.size2() << " n bins = " << n_bins_;


  unsigned year_index = 0;
  for (auto year_iterator = comparisons.begin(); year_iterator != comparisons.end(); ++year_iterator, ++year_index) {
    unsigned bin_iter = 0;
    bool first_comparison = true;
    for (observations::Comparison& comparison : year_iterator->second) {
      if (bin_iter >= n_bins_ && parameters_.Get(PARAM_BIN_LABELS)->has_been_defined()) {
        LOG_FATAL_P(PARAM_BIN_LABELS) << "Number of bins = " << n_bins_ << " there are more bins in the observation than supplied in the likleihood, please check that you have specified the correct number of bin labels";
      } else if (bin_iter >= n_bins_ && parameters_.GetTable(PARAM_COVARIANCE_MATRIX)->has_been_defined())
        LOG_FATAL_P(PARAM_COVARIANCE_MATRIX) <<"Number of bins = " << n_bins_ << " there are more columns in the user defined covariance matrix, please check that you have specified the correct dimensions for the covariance matrix";

      Double error_value = AdjustErrorValue(comparison.process_error_, comparison.error_value_) * error_value_multiplier_;
      if (first_comparison) {
        N.push_back(error_value);
        first_comparison = false;
      } else {
        // do a check that N is the same for all bins. if fail report an error
        if (N[year_index] != error_value) {
          LOG_FATAL() << "Logisitic normal requires the effective sample size to be the same for all bins. For category " << comparison.category_ <<
              " and age or length = " <<  comparison.age_ << " N was = "<< comparison.error_value_ << " which differs from the first N = " << N[year_index];
        }
      }
      Observed(year_index,bin_iter) = comparison.observed_;
      Expected(year_index,bin_iter) = comparison.expected_;
      ++bin_iter;
    }
  }
  LOG_FINEST() << "Finished filling matricies from comparisons.";

  Double ave_N = math::mean(N);
  for (auto this_N : N) {
    Double temp = sqrt(ave_N/this_N);
    LOG_FINEST() << "weight = " << temp;
    weights.push_back(temp);
  }

  for (unsigned i = 0; i < (n_bins_ - 1); ++i) {
    Kmat(i,i) = 1.0;
    Kmat(i,n_bins_ - 1) = -1.0;
  }
  LOG_FINEST() << "Are we here?";

  ublas::matrix<Double> Kmat_transposed = ublas::trans(Kmat);
  LOG_FINEST() << "Kmat_transposed  = " << Kmat_transposed;

  Kmat_transposed = ublas::prod(covariance_matrix_,Kmat_transposed);
  LOG_FINEST() << "Kmat_after covarmultiple transposed  = " << Kmat_transposed;

  Vmat = ublas::prod(Kmat,Kmat_transposed);
  LOG_FINEST() << "Vmat  = " << Vmat;

  bool inverted = InvertMatrix(Vmat,Vmat_inv);
  if (!inverted) {
    LOG_FATAL() << "could not invert convariance matrix matrix, if it is a user supplied covariance, check the matrix is invertable, else it is a code error";
  }
  LOG_FINEST() << "print Vmat_inv = " << Vmat_inv;

  // log and sum the logged observaton matrix
  Double Total_logged = 0.0;
  for (unsigned i = 0; i < Observed.size1(); ++i) {
    for (unsigned j = 0; j < Observed.size2(); ++j) {
      log_obs(i,j) = log(Observed(i,j));
      Total_logged += log_obs(i,j);
    }
  }
  LOG_FINEST() << "Are we here1?";

  Double detminant = det_fast(Vmat);

  // First value for the score
  Score = 0.5 * (Double)N.size() * (n_bins_ - 1) * log(2 * PI) + Total_logged + 0.5 * (Double)N.size() *  log(detminant);
  LOG_FINEST() << "total log = " << Total_logged << " log determinant = " << log(detminant);
  // check if any weights deviate form 1
  if (!math::all_ones(weights)) {
    vector<Double> log_weights;
    for (auto weight : weights)
      log_weights.push_back(log(weight));
    Score += (n_bins_ - 1) * math::Sum(log_weights);
  }

  // Sweep over the the obs and create this ww object.
  for (unsigned i = 0; i < ww.size1(); ++i) {
    for (unsigned j = 0; j < ww.size2(); ++j) {
      Double l_obs = Observed(i, j) / Observed(i, n_bins_ - 1);
      Double l_exp = Expected(i, j) / Expected(i, n_bins_ - 1);
      ww(i, j) = log(l_obs) - log(l_exp);
      //LOG_FINEST() << ww(i,j) << " i = "<< i << " j = " << j;
    }
  }

  // Now finish with a robustification
  LOG_FINEST() << "years = " << N.size() << " should equal rows of ww = " << ww.size1();
  bool first_year = true;
  if (robust_) {
    LOG_FINEST() << "Entering Robustification";
    // get the diaganol components of the inverted matrix
    vector<Double> diag_inv;

    for (unsigned i = 0; i < Vmat_inv.size1(); ++i)
      diag_inv.push_back((Double)Vmat_inv(i, i));

    unsigned year = 0;
    for (auto year_iterator = comparisons.begin(); year_iterator != comparisons.end(); ++year_iterator) {
      vector<Double> temp1, temp2;
      temp1.assign(n_bins_ - 1, 0.0);
      temp2.assign(1, 0.0);
      vector<Double> temp_holder, temp_holder1;
      Double temp3,robust_temp;
      // store in vectors for use in vector manipulation otherwise it would be shitty loops;
      for (unsigned row = 0; row < ww.size2(); ++row) {
        temp_holder.push_back(ww(year, row));
        temp_holder1.push_back(ww(year, row));
      }

      temp_holder = math::elem_prod(temp_holder, temp_holder);
      vector<Double> temp_holder2 = math::elem_prod(temp_holder, diag_inv);
      robust_temp = 0.0;
      for (unsigned i = 0; i < temp_holder2.size(); ++i) {
        Double exp_val = exp(-temp_holder2[i]) + 0.01;
        //LOG_FINEST() << log(exp_val);
        robust_temp += log(exp_val);
      }

      // Do some matrix multiplication
      for (unsigned i = 0; i < Vmat_inv.size1(); i++) {
        for (unsigned k = 0; k < Vmat_inv.size2(); k++)
          temp1[k] += temp_holder1[i] *  (Double)Vmat_inv(i, k);
      }
      for (unsigned i = 0; i < temp_holder1.size(); i++)
          temp2[0] += temp1[i] * temp_holder1[i];

      // Checking
      for(auto num : temp2)
        LOG_FINEST() << (Double)num;
      // Calculate the score
      temp3 = 0.5 / (weights[year] * weights[year]);
      Score += temp3 * (temp2[0] - math::Sum(temp_holder2) - robust_temp);
      LOG_FINEST() << "score after  = " << Score << " = " << (temp3 * temp2[0]) << " temp2 = " << math::Sum(temp_holder2) << " temp3 = " << AS_DOUBLE(robust_temp);
      for (observations::Comparison& comparison : year_iterator->second) {
        if (first_year) {
          comparison.score_ = Score;
          first_year = false;
        } else
          comparison.score_ = temp3 * (temp2[0] - math::Sum(temp_holder2) - robust_temp);
        break;
      }
      ++year;
    }
  } else {
    LOG_FINEST() << "Calculate score. number of years = " << N.size();
    unsigned year = 0;
    for (auto year_iterator = comparisons.begin(); year_iterator != comparisons.end(); ++year_iterator) {
      vector<Double> temp1, temp2;
      temp1.assign(n_bins_ - 1, 0.0);
      temp2.assign(1, 0.0);
      vector<Double> temp_holder, temp_holder1;
      Double temp3;
        for (unsigned row = 0; row < ww.size2(); ++row) {
        temp_holder.push_back(ww(year, row));
        temp_holder1.push_back(ww(year, row));
        //LOG_FINEST() << "equivalent t0 ww[i,] in Chris's code " << ww(year, row);
      }

      // Do some matrix multiplication
      for (unsigned i = 0; i < Vmat_inv.size1(); i++) {
        for (unsigned k = 0; k < Vmat_inv.size2(); k++)
          temp1[k] += temp_holder1[i] *  (Double)Vmat_inv(i, k);
      }

      LOG_FINEST() << "rows inverse mat = " << Vmat_inv.size1() << " cols = " << Vmat_inv.size2();

      //for(auto num: temp1)
      //  LOG_FINEST() << num;

      for (unsigned i = 0; i < temp_holder1.size(); i++)
          temp2[0] += temp1[i] * temp_holder1[i];

      // Calculate the score
      temp3 = 0.5 / (weights[year] * weights[year]);
      //LOG_FINEST() << "temp3 = " << AS_DOUBLE(temp3) << " temp2 " << AS_DOUBLE(temp2[0]) << " current NLL " << AS_DOUBLE(Score);
      Score += (temp3 * temp2[0]);
      for (observations::Comparison& comparison : year_iterator->second) {
        if (first_year) {
          comparison.score_ = Score;
          first_year = false;
        } else
          comparison.score_ = (temp3 * temp2[0]);
        break;
      }
      ++year;
    }
  }
}

/**
 * Grab the initial score for this likelihood
 *
 * @param comparisons A collection of comparisons passed by the observation
 * Returns a value for a likelihood that can calculated without the observation or expectation
 */

Double LogisticNormal::GetInitialScore(map<unsigned, vector<observations::Comparison> >& comparisons, unsigned year) {
  // This doesn't do anything in the logistic normal
  Double score = 0.0;

  return score;
}

/**
 * Simulate observed values
 *
 * @param comparisons A collection of comparisons passed by the observation
 */
void LogisticNormal::SimulateObserved(map<unsigned, vector<observations::Comparison> >& comparisons) {
  LOG_TRACE();
  // Generate a multivariate variabel X
  vector<Double> year_totals(comparisons.size(), 0.0);
  utilities::RandomNumberGenerator& rng = utilities::RandomNumberGenerator::Instance();
  unsigned year_storer = 0;

  if (!parameters_.GetTable(PARAM_COVARIANCE_MATRIX)->has_been_defined() && rho_[0] == 0.0) {
    LOG_FINEST() << "sigma_ = " << sigma_;
    // We just need to generate an independent multivariate normal distribtuion
    for (auto year_iterator = comparisons.begin(); year_iterator != comparisons.end(); ++year_iterator) {
      for (observations::Comparison& comparison : year_iterator->second) {
        comparison.observed_ = exp(rng.normal(0.0,1.0) * sigma_ + log(comparison.expected_));
        LOG_FINEST() << "random deviate " << rng.normal(0.0,1.0) << " age = " << comparison.age_ << " simuiulated val = " << comparison.observed_  << " expected = " << comparison.expected_ ;
        year_totals[year_storer] += comparison.observed_;
      }
      ++year_storer;
    }
  } else {
    if (!parameters_.GetTable(PARAM_COVARIANCE_MATRIX)->has_been_defined())
      calculate_covariance();

    LOG_FINEST() << "Printing Covariance top left triangle";
    for (unsigned k = 0; k < covariance_matrix_.size1(); ++k) {
      for (unsigned j = 0; j <= k; ++j) {
        LOG_FINEST() << "row = " << k << " col = " << j << " val = " << covariance_matrix_(j,k) << " ";
      }
    }
    if (!DoCholeskyDecmposition())
      LOG_FATAL()<< "Cholesky decomposition failed. Cannot continue Simulating from a logisitic-normal likelihood";
    // Calculate multivariate normal distribution
    year_storer = 0;
    vector<Double> normals(covariance_matrix_.size1(), 0.0);
    Double row_sum;
    for (auto year_iterator = comparisons.begin(); year_iterator != comparisons.end(); ++year_iterator) {
      unsigned age_category_iter = 0;
      unsigned nbins = year_iterator->second.size();
      for (unsigned i = 0; i < nbins; ++i) {
        normals[i] = rng.normal();
      }
      for (observations::Comparison& comparison : year_iterator->second) {
        row_sum = 0.0;
        for (unsigned j = 0; j < nbins; ++j) {
          row_sum += covariance_matrix_lt(j,age_category_iter) * normals[j];
        }
        comparison.observed_ = exp(row_sum + log(comparison.expected_));
        //LOG_FINEST() << " age = " << comparison.age_ << " simuiulated val = " << comparison.observed_  << " expected = " << comparison.expected_  << " multivariate offset = " << row_sum << " log expectations = " << log(comparison.expected_);
        year_totals[year_storer] += comparison.observed_;
        ++age_category_iter;
      }
      ++year_storer;
    }
  }
  // Do the logistic transformation to get our desired values.
  year_storer = 0;
  for (auto year_iterator = comparisons.begin(); year_iterator != comparisons.end(); ++year_iterator) {
    LOG_FINEST() << "year = " << year_iterator->first;
    for (observations::Comparison& comparison : year_iterator->second) {
      comparison.observed_ /= year_totals[year_storer];
      LOG_FINEST() << "Simulated val = " << comparison.observed_ << " expected = " << comparison.expected_;
    }
    ++year_storer;
  }
  LOG_FINEST() << "check out the totals";
  for(auto num :year_totals)
    LOG_FINEST() << num ;
}



void LogisticNormal::calculate_covariance() {
  LOG_TRACE();
  unsigned n_phi = rho_.size();
  vector<Double> rhovec;
  vector<vector<Double>> covar;
  covar.resize(n_bins_, vector<Double>(n_bins_, 0.0));

  LOG_FINEST() << "covariance rows = " << covar.size() << " cols = " << covar[0].size();
  // initialise covar as the identity matrix
  for (unsigned diag = 0; diag < covar.size(); ++ diag)
    covar[diag][diag] = 1.0 * sigma_ * sigma_;

  if (rho_[0] == 0.0) {
    // Do nothing all zeros in the off diagonal
  } else if (sexed_ && sep_by_sex_) {

      LOG_FINEST() << "Calculating covar for sexed patrition with sepbysex = true, covar dims = " << covar.size()  << " " << covar[0].size();
      for (unsigned i = 0; i <= 1; ++i) {
        rhovec = GetRho(rho_,unique_bins_,arma_);
          for (unsigned diag = 0; diag < (unique_bins_ - 1); ++diag) {
            for (int row = 0; row < (int)unique_bins_; ++row) {
              for (int col = 0; col < (int)unique_bins_; ++col) {
                if (fabs(row - col) == diag + 1)
                  covar[row + (i * unique_bins_)][col + (i * unique_bins_)] = rhovec[diag] * sigma_ * sigma_;
              }
            }
          }
        }

    } else if (sexed_ && sex_lag_) {
      vector<int> binlab;
      int ages = 1;
      for (unsigned i = 1; i <= n_bins_; ++i, ++ages) {
        if (ages > (int)unique_bins_)
          ages = 1;
        binlab.push_back(ages);
        cout << ages << " ";
      }
      if (arma_ && n_phi == 2)
        rhovec = GetRho(rho_,unique_bins_,arma_);
      else if (!arma_ && n_phi == 2)
        rhovec = GetRho(rho_,unique_bins_ + 1,arma_);
        vector<int> col_vec, row_vec, offset_vec;

        for (unsigned row = 0; row < n_bins_; ++row) {
          for (unsigned col = 0; col < n_bins_; ++col) {
            col_vec.push_back(binlab[col]);
            row_vec.push_back(binlab[row]);
          }
        }
        // now calculate an offset vector;
        for (unsigned index = 0; index < col_vec.size(); ++index) {
          offset_vec.push_back(fabs(row_vec[index] - col_vec[index]) + 1.0);
          cout << offset_vec[index] << " ";
        }
        if (covar.size() * covar[0].size() != offset_vec.size())
          LOG_CODE_ERROR() << "covar.size() * covar[0].size() != offset_vec.size(). get_mat_index will fail as index vector needs to be same size as rows * col of the mat.";
        for (unsigned index = 1; index <  unique_bins_; ++index) {
          for (unsigned index2 = 0; index2 < offset_vec.size(); ++index2) {
            if ((int)index == offset_vec[index2]) {
              vector<int> indexes = math::get_mat_index((int)covar.size(),(int)covar[0].size(), index2);
              covar[indexes[0]][indexes[1]] = rhovec[index - 1];
              //cout << "row = " << indexes[0] << " col = " << indexes[1] << " index = " << index2 << endl;
            }
          }
        }


      // Add the identity mat
      for (unsigned row = 0; row < n_bins_; ++row)
        covar[row][row] = 1.0;
      // Store in the Covariance matrix
      for (unsigned row = 0; row < n_bins_; ++row) {
        for (unsigned col = 0; col < n_bins_; ++col) {
          covar[row][col] *= sigma_ * sigma_;
        }
      }
    } else {
      // Unisex or sexed but treated like a single covariance matrix
      rhovec = GetRho(rho_,unique_bins_,arma_);
      for (int diag = 0; diag < (int)unique_bins_; ++ diag) {
        for (int row = 0; row <  (int)n_bins_; ++ row) {
          for (int col = 0; col < (int)n_bins_; ++ col) {
            if (fabs(row - col) == diag + 1) {
              covar[row][col] = rhovec[diag] * sigma_ * sigma_;
            }
          }
        }
      }
    }
    // I am struggling to intialise this matrix with ints. so do the way I know will work
    ublas::matrix<Double> temp_covar(covar.size(),covar.size(),0.0);
    for(unsigned i = 0; i < covar.size(); ++i){
      for(unsigned j = 0; j < covar[i].size(); ++j) {
        temp_covar(i,j) = covar[i][j] ;
      }
    }
    covariance_matrix_ = temp_covar;
    LOG_FINEST() << "Finished building covariance matrix";
}

vector<Double> LogisticNormal::GetRho(vector<Double>& Phi, unsigned nBin, bool ARMA) {
  LOG_TRACE();
  // declare all variables that will be used in this function
  vector<Double> rhovec(nBin, 0.0);
  if (Phi.size() == 1) {
    LOG_FINEST() <<  "Found single Rho parameter = " <<Phi[0] << " number of bins = " << nBin;
    //calculation of AR(1) acf for  LN2
    for(unsigned i = 1; i <= nBin - 1; i++)
      rhovec[i - 1] = pow(Phi[0],i);
    LOG_FINEST() << "Finished building rho vec";
  } else {
    vector<Double> ar, ma,final_acf,Cor;
    vector<vector<Double> > A, ind;
    vector<Double> psi, theta, Acf;
    // we are doing ARMAacf function
    unsigned p, q, r;
    if (ARMA) {
      q = 1;
      p = 1;
      ar.push_back(Phi[0]);
    } else {
      q = 0;
      p = 2;
      ar = Phi;
    }
    r = fmax(p, q + 1);
    if (p > 0) {
      if (r > 1) {
        if (r > p) {
          LOG_FINEST() << "calculating rho from an ARMA(1,1) process";
          for (unsigned i = 0; i < (r - p); ++i)
            ar.push_back(0.0);
          p = r;
        }
        LOG_FINEST() << "Structureing A";

        A.resize(p + 1, vector<Double>(2 * p + 1, 0.0));
        ind.resize(2 * p + 1, vector<Double>(p + 1, 0.0));
        for (int i = 0; i < (int)ind.size(); ++i) {
          for (int j = 0; j < (int)ind[i].size(); ++j) {
            ind[i][0] = i + 1;
            ind[i][1] = (i + 1) + (j + 1) - (i + 1);
          }
        }

        for (unsigned i = 0; i < A.size(); ++i) {
          A[i][i] = 1.0;
           A[i][i + 1] = -ar[0];
           A[i][i + 2] = -ar[1];
        }
        LOG_FINEST() << "Populate A. the second ar value" << ar[1];
        ublas::matrix<Double> A_eig1(3,3,0.0);
        ublas::matrix<Double> A_eig_inv(3,3,0.0);

        vector<Double> rhs(3,0.0);
        // initialise rhs, which will be used to solve the following problem, that is Ax = b where b = rhs, so x = A^-1 b
        rhs[0] = 1.0;
        rhs[1] = 0.0;
        rhs[2] = 0.0;
        if (q > 0) {
          LOG_FINEST() << "Calculate ma coef";
          // deal with the ma coeffecient
          psi.push_back(1.0);
          psi.push_back(Phi[0] + Phi[1]);
          theta.push_back(1.0);
          theta.push_back(Phi[1]);
          for (unsigned i = 0; i <= q; ++i)
            theta.push_back(0.0);
          // Calculate rhs
          for (unsigned i = 0; i <= q; ++i) {
            Double x1, x2;
            x1 = psi[0] * theta[i];
            x2 = psi[1] * theta[i + q];
            Double val = 0.0;
            if (!utilities::To<Double, Double>(math::Sum({ x1, x2 }), val))
              LOG_CODE_ERROR() << " val " << math::Sum({ x1, x2 }) << " could not be converted in to a Double";
            rhs[i] = val;
          }
          rhs[2] = 0.0;
        }
        LOG_FINEST() << "Calculate seq parameter";
        // Use the eigen library yo solve the inverse of for A with known vector B
        //vector<Double> Ind;
        vector<unsigned> seq;
        for (unsigned i = 0; i <= p; ++i) {
          seq.push_back(p - i);
        }
        for (unsigned i = 0; i <= p; ++i) {
          for (unsigned j = 0; j <= p; ++j) {
            //LOG_FINEST() << ": i = " << i << " j = " << j << " i index = " << seq[i] << " j index = " << seq[j] << " mat value = " << A[seq[i]][seq[j]];
            Double val = 0.0;
            if (j == 2) {
              if (!utilities::To<Double, Double>(A[i][j], val))
                LOG_CODE_ERROR() << "variable = " << val << " could not be converted in to a Double";
              A_eig1(i,j) = val;
            } else {
              if (!utilities::To<Double, Double>(A[i][j] + A[i][2 * p  - j], val))
                LOG_CODE_ERROR() << "variable = " << val << " could not be converted in to a Double";
              A_eig1(i,j) = val;
            }
          }
        }
        for (unsigned i = 0; i <= p; ++i) {
          for (unsigned j = 0; j <= p ; ++j) {
            A_eig1(i,j)=  A_eig1(seq[i],seq[j]);
          }
        }
        // the bodge
        A_eig1(1,2) = 0.0;


        LOG_FINEST() << "Check A mat that we are inverting\n" << A_eig1;

        // Now try mine
        bool inverted = InvertMatrix(A_eig1,A_eig_inv);
        if (!inverted) {
          LOG_FATAL() << "could not invert convariance matrix matrix, if it is a user supplied covariance, check the matrix is invertable, else it is a code error";
        }
        // Matrix multiplication to solve for x
        vector<Double> result(A_eig1.size1(),0.0);

        for (unsigned i = 0; i < A_eig_inv.size1(); i++) {
         for (unsigned k = 0; k < A_eig_inv.size2(); k++) {
            result[i] += rhs[k] * A_eig_inv(i,k);
         }
        }
        LOG_FINEST() << "solution = ";
        for(auto num : result)
          LOG_FINEST() << num;

        for (unsigned i = 1; i <= 2; ++i) {
          final_acf.push_back(result[i] / result[0]);
        }
        LOG_FINEST() << "Final Acf";
        for (auto num : final_acf)
          LOG_FINEST() << num << " ";

        Cor = RecursiveFilter(ar, nBin, final_acf);

        // Add the initial coeffiecients to the beginning of the series.
        Cor.insert(Cor.begin(), final_acf[1]);
        Cor.insert(Cor.begin(), final_acf[0]);
        // Print results to screen
        vector<Double>::const_iterator first = Cor.begin();
        vector<Double>::const_iterator last = Cor.begin() + nBin;
        vector<Double> corvec(first, last);
        rhovec = corvec;
        LOG_FINEST() << " Print rhovec";
        for (auto num : rhovec)
          LOG_FINEST()  << num << " ";

      }

    }
  }
return rhovec;
}

vector<Double> LogisticNormal::RecursiveFilter(vector<Double>& ar_coef, unsigned nBins, vector<Double>& initial_vals) {
  LOG_TRACE();
  vector<Double> store_vec(nBins + 1,0.0);

  if (ar_coef.size() > 2) {
    LOG_CODE_ERROR() <<  "RecursiveFilter(): has not been coded for more than 2 AR coeffiecients, ar_coef.size() > 2" << endl;
  }
  store_vec[0] = initial_vals[1];
  store_vec[1] = initial_vals[0];
  for (unsigned i = 1; i < nBins + 1; ++i) {
    if (i == 1) {
      store_vec[i] =   store_vec[i - 1] *ar_coef[0]  + store_vec[i] *  ar_coef[1];
    } else {
      store_vec[i] = store_vec[i - 1] *  ar_coef[0] + store_vec[i - 2] * ar_coef[1];
    }
    LOG_FINEST() << "value = " << store_vec[i];
  }
  // remove the first value
  store_vec.erase(store_vec.begin());
  return store_vec;
}

/**
 * Perform cholesky decomposition on our covariance
 * matrix before it's used in the simulations aspect
 * when simualting a multivariate normal distribution
 *
 * @return true on success, false on failure
 */

bool LogisticNormal::DoCholeskyDecmposition() {
  if (covariance_matrix_.size1() != covariance_matrix_.size2() )
      LOG_ERROR() << "Invalid covariance matrix (size1!=size2)";
    unsigned matrix_size1 = covariance_matrix_.size1();
    covariance_matrix_lt = covariance_matrix_;

    for (unsigned i = 0; i < matrix_size1; ++i) {
      for (unsigned j = 0; j < matrix_size1; ++j) {
        covariance_matrix_lt(i,j) = 0.0;
      }
    }

    for (unsigned i = 0; i < matrix_size1; ++i) {
      covariance_matrix_lt(i,i) = 1.0;
    }

    if (covariance_matrix_(0,0) < 0)
      return false;
    Double sum = 0.0;

    covariance_matrix_lt(0,0) = sqrt(covariance_matrix_(0,0));

    for (unsigned i = 1; i < matrix_size1; ++i)
      covariance_matrix_lt(i,0) = covariance_matrix_(i,0)/covariance_matrix_lt(0,0);

    for (unsigned i = 1; i < matrix_size1; ++i) {
      sum = 0.0;
      for (unsigned j = 0; j < i; ++j)
        sum += covariance_matrix_lt(i,j) * covariance_matrix_lt(i,j);

      if (covariance_matrix_(i,i) <= sum)
        return false;
      covariance_matrix_lt(i,i) = sqrt(covariance_matrix_(i,i) - sum);
      for (unsigned j = i+1; j < matrix_size1; ++j) {
        sum = 0.0;
        for (unsigned k = 0; k < i; ++k)
          sum += covariance_matrix_lt(j,k) * covariance_matrix_lt(i,k);
        covariance_matrix_lt(j,i) = (covariance_matrix_(j,i) - sum) / covariance_matrix_lt(i,i);
      }
    }
    sum = 0.0;
    for (unsigned i = 0; i < (matrix_size1 - 1); ++i)
      sum += covariance_matrix_lt(matrix_size1 - 1,i) * covariance_matrix_lt(matrix_size1-1,i);
    if (covariance_matrix_(matrix_size1 - 1, matrix_size1 - 1) <= sum)
      return false;
    covariance_matrix_lt(matrix_size1 - 1, matrix_size1 - 1) = sqrt(covariance_matrix_(matrix_size1 - 1, matrix_size1 - 1) - sum);

   return true;
}

/*
 * Invert a square boost matrix
*/
bool LogisticNormal::InvertMatrix(const ublas::matrix<Double>& input, ublas::matrix<Double>& inverse) {
  typedef ublas::permutation_matrix<std::size_t> pmatrix;
  // create a working copy of the input
  ublas::matrix<Double> A(input);

  // create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());

  // perform LU-factorization
  int res = ublas::lu_factorize(A, pm);
  if (res != 0)
    return false;

  // create identity matrix of "inverse"
  inverse.assign(ublas::identity_matrix<Double> (A.size1()));

  // backsubstitute to get the inverse
  ublas::lu_substitute(A, pm, inverse);
  return true;
}

/*
 * Calcuate the determinant of a square boost matrix
*/
Double LogisticNormal::det_fast(const ublas::matrix<Double>& matrix) {
  // create a working copy of the input
 ublas::matrix<Double> mLu(matrix);
 ublas::permutation_matrix<size_t> pivots(matrix.size1());

 auto isSingular = ublas::lu_factorize(mLu, pivots);

   if (isSingular)
       return static_cast<Double>(0);

   Double det = static_cast<Double>(1);
   for (size_t i = 0; i < pivots.size(); ++i)
   {
       if (pivots(i) != i)
           det *= static_cast<Double>(-1);

       det *= mLu(i, i);
   }
   return det;
}


} /* namespace likelihoods */
} /* namespace niwa */

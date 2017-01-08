/*
 * LogisticNormal.h
 *
 *  Created on: Oct 26, 2016
 *      Author: Zaita
 */

#ifndef SOURCE_LIKELIHOODS_CHILDREN_LOGISTICNORMAL_H_
#define SOURCE_LIKELIHOODS_CHILDREN_LOGISTICNORMAL_H_

#include "Likelihoods/Likelihood.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>

namespace niwa {
namespace likelihoods {
namespace ublas = boost::numeric::ublas;


class LogisticNormal : public niwa::Likelihood {
public:
  // Methods
  LogisticNormal();
  virtual                     ~LogisticNormal();
  void                        DoValidate() override final;
  Double                      AdjustErrorValue(const Double process_error, const Double error_value) override final;
  void                        SimulateObserved(map<unsigned, vector<observations::Comparison> >& comparisons) override final;
  Double                      GetInitialScore(map<unsigned, vector<observations::Comparison> >& comparisons, unsigned year) override final;
  void                        GetScores(map<unsigned, vector<observations::Comparison> >& comparisons) override final;

protected:
  // Estimable parameters
  Double                      sigma_;
  vector<unsigned>            bins_;
  vector<Double>              rho_;
  bool                        arma_;
  bool                        robust_;
  unsigned                    n_bins_;
  bool                        sep_by_sex_;
  bool                        sex_lag_;
  bool                        sexed_;
  unsigned                    unique_bins_;

  // Covariance containers
  ublas::matrix<Double>       covariance_matrix_;
  ublas::matrix<Double>       covariance_matrix_lt;
  parameters::Table*          covariance_table_ = nullptr;
  // Methods
  void                        calculate_covariance();
  vector<Double>              GetRho(vector<Double>& Phi, unsigned nBin, bool ARMA);
  vector<Double>              RecursiveFilter(vector<Double>& ar_coef, unsigned nBins, vector<Double>& initial_vals);
  bool                        DoCholeskyDecmposition();
  bool                        InvertMatrix(const ublas::matrix<Double>& input, ublas::matrix<Double>& inverse);
  Double                      det_fast(const ublas::matrix<Double>& matrix);


};

} /* namespace likelihoods */
} /* namespace niwa */
#endif /* SOURCE_LIKELIHOODS_CHILDREN_LOGISTICNORMAL_H_ */

// round_regularized_params_preload.cpp
// -----------------------------------------------------------------------------
// Same fast logic as before, but now C++ pulls data out of `results` and
// `hrf_params`—no manual prep in R.
// -----------------------------------------------------------------------------
//
// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
#include <limits>

using namespace Rcpp;

// ── helper: nearest grid combo ───────────────────────────────────────────────
inline int nearest_idx(double a1_i, double b1_i, double c_i,
                       const NumericVector& a_ok,
                       const NumericVector& b_ok,
                       const NumericVector& c_ok,
                       double a_rng, double b_rng)
{
  double best = std::numeric_limits<double>::infinity();
  int    best_j = -1;
  for (int j = 0; j < c_ok.size(); ++j) {
    if (c_ok[j] != c_i) continue;
    double da = (a1_i - a_ok[j]) / a_rng;
    double db = (b1_i - b_ok[j]) / b_rng;
    double d2 = da*da + db*db;
    if (d2 < best) { best = d2; best_j = j; }
  }
  return best_j;
}

// [[Rcpp::export]]
DataFrame round_regularized_params(
    List        results,                 // model_outputs
    DataFrame   hrf_params,              // grid dataframe (required)
    std::string model        = "A",      // defaults stay the same
    bool        inside_grid  = true,
    int         max_voxels   = -1)
{
  // ── choose column name ─────────────────────────────────────────────────────
  std::string colname;
  if      (model == "A") colname = "param_pred_A";
  else if (model == "B") colname = "param_pred_B";
  else   Rcpp::stop("model must be \"A\" or \"B\"");

  // ── pull sub-data-frames from results list ────────────────────────────────
  DataFrame df_a1 = as<DataFrame>( as<List>(results["a1"])["results_WLS"] );
  DataFrame df_b1 = as<DataFrame>( as<List>(results["b1"])["results_WLS"] );
  DataFrame df_c  = as<DataFrame>( as<List>(results["c"])["results_WLS"] );
  // parameters & metadata
  NumericVector   a1_vals = df_a1[colname];
  NumericVector   b1_vals = df_b1[colname];
  NumericVector   c_vals  = df_c [colname];

  IntegerVector   voxel   = df_a1["voxel"];
  CharacterVector subject = df_a1["subject"];
  LogicalVector   mask    = df_a1["mask"];

  const int n_total = a1_vals.size();
  const int n       = (max_voxels > 0 && max_voxels < n_total) ? max_voxels : n_total;

  // ── initial rounding ──────────────────────────────────────────────────────
  NumericVector new_a1 = clone(a1_vals);
  NumericVector new_b1(n_total);
  NumericVector new_c (n_total);
  const double sixth = 1.0 / 6.0;

  for (int i = 0; i < n_total; ++i) {
    new_a1[i] = std::round(a1_vals[i]);                  // nearest int
    new_b1[i] = std::round(b1_vals[i] * 4.0) / 4.0;      // nearest 0.25
    new_c [i] = (std::abs(c_vals[i]) <= std::abs(c_vals[i] - sixth))
      ? 0.0 : sixth;                         // 0 or 1/6
  }

  // ── optional grid snapping ────────────────────────────────────────────────
  if (inside_grid) {
    NumericVector ok_a1 = hrf_params["a1"];
    NumericVector ok_b1 = hrf_params["b1"];
    NumericVector ok_c  = hrf_params["c"];

    const double a_rng = max(ok_a1) - min(ok_a1);
    const double b_rng = max(ok_b1) - min(ok_b1);

    for (int i = 0; i < n; ++i) {
      // skip if already valid
      bool exact = false;
      for (int j = 0; j < ok_c.size(); ++j) {
        if (ok_c[j] == new_c[i] &&
            ok_a1[j] == new_a1[i] &&
            ok_b1[j] == new_b1[i]) { exact = true; break; }
      }
      if (exact) continue;

      int best = nearest_idx(new_a1[i], new_b1[i], new_c[i],
                             ok_a1, ok_b1, ok_c, a_rng, b_rng);
      if (best != -1) {
        new_a1[i] = ok_a1[best];
        new_b1[i] = ok_b1[best];
      }
    }
  }

  // ── return ────────────────────────────────────────────────────────────────
  return DataFrame::create(
    _["a1"]      = new_a1,
    _["b1"]      = new_b1,
    _["c"]       = new_c,
    _["voxel"]   = voxel,
    _["subject"] = subject,
    _["mask"]    = mask);
}

// [[Rcpp::export]]
void say_hello() {
  Rcpp::Rcout << "Hello from the pre-load C++ HRF helper!\n";
}

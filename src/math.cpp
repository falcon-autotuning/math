#include "falcon-math/math.hpp"
#include <algorithm>
#include <armadillo>
#include <cmath>
#include <ensmallen.hpp>
#include <numeric>

namespace falcon {
namespace math {

// Internal helper for numeric differentiation
template <typename Function>
void numeric_gradient(Function &f, const arma::mat &parameters,
                      arma::mat &gradient, double eps = 1e-6) {
  gradient.set_size(parameters.n_rows, parameters.n_cols);
  arma::mat params = parameters;
  double base_val = f.Evaluate(params);
  for (size_t i = 0; i < parameters.n_elem; ++i) {
    double old_val = params(i);
    params(i) += eps;
    double new_val = f.Evaluate(params);
    gradient(i) = (new_val - base_val) / eps;
    params(i) = old_val;
  }
}

// Objective function for 1D curve fitting
class CurveFit1DObjective {
  func1D model;
  std::vector<double> x;
  std::vector<double> y;

public:
  CurveFit1DObjective(func1D m, const std::vector<double> &xv,
                      const std::vector<double> &yv)
      : model(m), x(xv), y(yv) {}

  double Evaluate(const arma::mat &parameters) {
    std::vector<double> params(parameters.n_elem);
    for (size_t i = 0; i < parameters.n_elem; ++i)
      params[i] = parameters(i);

    std::vector<double> y_fit = model(x, params);
    double rss = 0;
    for (size_t i = 0; i < y.size(); ++i) {
      double diff = y[i] - y_fit[i];
      rss += diff * diff;
    }
    return rss;
  }

  void Gradient(const arma::mat &parameters, arma::mat &gradient) {
    numeric_gradient(*this, parameters, gradient);
  }
};

// Objective function for 2D curve fitting
class CurveFit2DObjective {
  func2D model;
  std::vector<std::vector<double>> x;
  std::vector<std::vector<double>> y;
  std::vector<std::vector<double>> z;

  std::optional<std::vector<std::pair<double, double>>> bounds;

public:
  CurveFit2DObjective(
      func2D m, const std::vector<std::vector<double>> &xv,
      const std::vector<std::vector<double>> &yv,
      const std::vector<std::vector<double>> &zv,
      std::optional<std::vector<std::pair<double, double>>> b = std::nullopt)
      : model(m), x(xv), y(yv), z(zv), bounds(b) {}

  double Evaluate(const arma::mat &parameters) {
    std::vector<double> params(parameters.n_elem);
    double penalty = 0.0;

    for (size_t i = 0; i < parameters.n_elem; ++i) {
      double val = parameters(i);
      if (bounds.has_value() && i < bounds->size()) {
        double low = (*bounds)[i].first;
        double high = (*bounds)[i].second;
        if (val < low) {
          penalty += 1e8 * (low - val) * (low - val);
          val = low; // Clamp for evaluation to prevent model instability
        } else if (val > high) {
          penalty += 1e8 * (val - high) * (val - high);
          val = high; // Clamp for evaluation
        }
      }
      params[i] = val;
    }

    std::vector<std::vector<double>> p_wrapped = {params};
    auto z_fit_res = model(x, y, p_wrapped);

    double rss = 0;
    for (size_t i = 0; i < z.size(); ++i) {
      for (size_t j = 0; j < z[i].size(); ++j) {
        double diff = z[i][j] - z_fit_res[i][j];
        if (std::isnan(diff) || std::isinf(diff)) {
          return 1e18; // Large value instead of NaN
        }
        rss += diff * diff;
      }
    }

    return rss + penalty;
  }

  void Gradient(const arma::mat &parameters, arma::mat &gradient) {
    numeric_gradient(*this, parameters, gradient);
  }
};

FALCON_MATH_API curvefit_result curvefit1D(func1D model,
                                          const std::vector<double> &x,
                                          const std::vector<double> &y,
                                          fitting_parameters params) {
  size_t num_params = params.initial_guess.has_value()
                          ? params.initial_guess->size()
                          : 4; // Default guess size?
  arma::mat p_arma(num_params, 1);
  if (params.initial_guess) {
    for (size_t i = 0; i < num_params; ++i)
      p_arma(i) = (*params.initial_guess)[i];
  } else {
    p_arma.fill(1.0);
  }

  CurveFit1DObjective obj(model, x, y);
  ens::L_BFGS optimizer;

  try {
    optimizer.Optimize(obj, p_arma);
  } catch (const std::exception &e) {
    return {{}, e.what(), false, 0.0};
  }

  curvefit_result res;
  res.coefficients.resize(num_params);
  for (size_t i = 0; i < num_params; ++i)
    res.coefficients[i] = p_arma(i);
  res.success = true;

  // Calculate R-squared
  double mean_y = std::accumulate(y.begin(), y.end(), 0.0) / y.size();
  double ss_tot = 0, ss_res = 0;
  std::vector<double> y_fit = model(x, res.coefficients);
  for (size_t i = 0; i < y.size(); ++i) {
    ss_tot += (y[i] - mean_y) * (y[i] - mean_y);
    ss_res += (y[i] - y_fit[i]) * (y[i] - y_fit[i]);
  }
  res.r_squared = 1.0 - (ss_res / ss_tot);

  return res;
}

FALCON_MATH_API curvefit_result curvefit2D(
    func2D model, const std::vector<std::vector<double>> &x,
    const std::vector<std::vector<double>> &y,
    const std::vector<std::vector<double>> &z, fitting_parameters params) {
  size_t num_params =
      params.initial_guess.has_value() ? params.initial_guess->size() : 12;
  arma::mat p_arma(num_params, 1);
  if (params.initial_guess) {
    for (size_t i = 0; i < num_params; ++i)
      p_arma(i) = (*params.initial_guess)[i];
  } else {
    p_arma.fill(1.0);
  }

  CurveFit2DObjective obj(model, x, y, z, params.bounds);

  // User requested DE. Using robust parameters.
  ens::DE optimizer(200, 3000, 0.6, 0.8, -1.0);

  try {
    optimizer.Optimize(obj, p_arma);
  } catch (const std::exception &e) {
    return {{}, e.what(), false, 0.0};
  }

  curvefit_result res;
  res.coefficients.resize(num_params);
  for (size_t i = 0; i < num_params; ++i)
    res.coefficients[i] = p_arma(i);
  res.success = true;

  // R-squared
  double sum_z = 0;
  size_t count = 0;
  for (const auto &row : z) {
    for (double val : row) {
      sum_z += val;
      count++;
    }
  }
  double mean_z = sum_z / count;
  double ss_tot = 0, ss_res = 0;
  auto z_fit = model(x, y, {res.coefficients});
  for (size_t i = 0; i < z.size(); ++i) {
    for (size_t j = 0; j < z[i].size(); ++j) {
      ss_tot += (z[i][j] - mean_z) * (z[i][j] - mean_z);
      ss_res += (z[i][j] - z_fit[i][j]) * (z[i][j] - z_fit[i][j]);
    }
  }
  res.r_squared = 1.0 - (ss_res / ss_tot);

  return res;
}

// Predefined models
FALCON_MATH_API double sigmoid(double x, double A, double x0, double k,
                               double b) {
  return A / (1.0 + std::exp(-k * (x - x0))) + b;
}

FALCON_MATH_API double piecewise_linear(double x, double x0, double x1,
                                        double m1, double m2, double y0) {
  if (x < x0)
    return m1 * (x - x0) + y0;
  if (x < x1)
    return m2 * (x - x0) + y0;
  return m2 * (x1 - x0) + y0;
}

FALCON_MATH_API double
channel_accumulation_2d(double x, double y, double cx, double cy, double cm,
                        double cr, double m1, double m2, double m3, double bx,
                        double dx, double dy, double dm, double dr) {
  auto sig = [](double val, double mm, double kk, double bb) {
    double arg = -kk * (val - bb);
    if (arg > 50.0)
      return 0.0; // Avoid overflow, return small value
    if (arg < -50.0)
      return mm; // Avoid underflow
    return mm / (1.0 + std::exp(arg));
  };

  double b1 = dx - m1 * dy;
  double b2 = dx - m2 * dy;
  double by = (bx - b2) / m2;
  double b3 = bx - m3 * by;

  double diag = 0;
  double val2 = m2 * y + b2;
  double val3 = m3 * y + b3;
  double val1 = m1 * y + b1;

  if (val2 <= dx && val2 >= bx) {
    diag = std::max(diag, sig(x, dm, dr, val2));
  }
  if (val3 <= bx) {
    diag = std::max(diag, sig(x, dm, dr, val3));
  }
  if (val1 >= dx) {
    diag = std::max(diag, sig(x, dm, dr, val1));
  }

  double main_sig = sig(y, cm, cr, cy) * sig(x, cm, cr, cx);
  return std::max(main_sig, diag);
}

} // namespace math
} // namespace falcon

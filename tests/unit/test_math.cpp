#include "falcon-math/math.hpp"
#include <armadillo>
#include <cmath>
#include <gtest/gtest.h>
#include <random>
#include <vector>

using namespace falcon::math;

TEST(MathTest, SigmoidEvaluation) {
  double val = sigmoid(0.0, 1.0, 0.0, 1.0, 0.0);
  EXPECT_NEAR(val, 0.5, 1e-6);
}

TEST(MathTest, PiecewiseLinearEvaluation) {
  // x0=0, x1=10, m1=0, m2=1, y0=0
  EXPECT_NEAR(piecewise_linear(-5, 0, 10, 0, 1, 0), 0, 1e-6);
  EXPECT_NEAR(piecewise_linear(5, 0, 10, 0, 1, 0), 5, 1e-6);
  EXPECT_NEAR(piecewise_linear(15, 0, 10, 0, 1, 0), 10, 1e-6);
}

TEST(MathTest, CurveFit1DSimple) {
  // Linear model: y = m*x + b
  auto model = [](std::vector<double> x, std::vector<double> p) {
    std::vector<double> y(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
      y[i] = p[0] * x[i] + p[1];
    }
    return y;
  };

  std::vector<double> x = {0, 1, 2, 3, 4, 5};
  std::vector<double> y = {1, 3, 5, 7, 9, 11}; // m=2, b=1

  fitting_parameters params;
  params.initial_guess = {1.0, 0.0};

  auto res = curvefit1D(model, x, y, params);

  ASSERT_TRUE(res.success);
  EXPECT_NEAR(res.coefficients[0], 2.0, 1e-3);
  EXPECT_NEAR(res.coefficients[1], 1.0, 1e-3);
  EXPECT_GT(res.r_squared, 0.99);
}

TEST(MathTest, CurveFitSigmoid) {
  auto model = [](std::vector<double> x, std::vector<double> p) {
    std::vector<double> y(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
      y[i] = sigmoid(x[i], p[0], p[1], p[2], p[3]);
    }
    return y;
  };

  std::vector<double> x;
  std::vector<double> y;
  double A = 1.0, x0 = 5.0, k = 2.0, b = 0.1;
  for (double val = 0; val <= 10; val += 0.5) {
    x.push_back(val);
    y.push_back(sigmoid(val, A, x0, k, b));
  }

  fitting_parameters params;
  params.initial_guess = {0.8, 4.5, 1.5, 0.0};

  auto res = curvefit1D(model, x, y, params);

  ASSERT_TRUE(res.success);
  EXPECT_NEAR(res.coefficients[0], A, 1e-2);
  EXPECT_NEAR(res.coefficients[1], x0, 1e-2);
  EXPECT_NEAR(res.coefficients[2], k, 1e-2);
  EXPECT_NEAR(res.coefficients[3], b, 1e-2);
}

TEST(MathTest, CurveFit2DChannelAccumulation) {
  auto model = [](std::vector<std::vector<double>> x,
                  std::vector<std::vector<double>> y,
                  std::vector<std::vector<double>> p) {
    std::vector<double> params = p[0];
    std::vector<std::vector<double>> z(x.size(),
                                       std::vector<double>(x[0].size()));
    for (size_t i = 0; i < x.size(); ++i) {
      for (size_t j = 0; j < x[i].size(); ++j) {
        z[i][j] = channel_accumulation_2d(
            x[i][j], y[i][j], params[0], params[1], params[2], params[3],
            params[4], params[5], params[6], params[7], params[8], params[9],
            params[10], params[11]);
      }
    }
    return z;
  };

  // Generate grid
  std::vector<std::vector<double>> x(10, std::vector<double>(10));
  std::vector<std::vector<double>> y(10, std::vector<double>(10));
  for (int i = 0; i < 10; ++i) {
    for (int j = 0; j < 10; ++j) {
      x[i][j] = i * 0.1;
      y[i][j] = j * 0.1;
    }
  }

  // Parameters for channel_accumulation_2d (12 total)
  // double cx, double cy, double cm, double cr,
  // double m1, double m2, double m3,
  // double bx, double dx, double dy,
  // double dm, double dr
  std::vector<double> true_p = {
      0.5, 0.5,  1.0, 10.0, // cx, cy, cm, cr
      1.0, -1.0, 0.5,       // m1, m2, m3
      0.2, 0.8,  0.5,       // bx, dx, dy
      0.5, 5.0              // dm, dr
  };

  std::vector<std::vector<double>> z = model(x, y, {true_p});

  fitting_parameters params;
  params.initial_guess = true_p;
  // Perturb initial guess slightly to verify convergence
  for (double &val : *params.initial_guess)
    val *= 1.05;

  // Provide reasonable bounds to help DE converge
  std::vector<std::pair<double, double>> bounds;
  for (double val : true_p) {
    bounds.push_back({val - 2.0, val + 2.0});
  }
  params.bounds = bounds;

  // Set seed for reproducibility
  arma::arma_rng::set_seed(42);

  auto res = curvefit2D(model, x, y, z, params);

  ASSERT_TRUE(res.success);
  for (size_t i = 0; i < true_p.size(); ++i) {
    EXPECT_NEAR(res.coefficients[i], true_p[i], 2.5);
  }
  EXPECT_GT(res.r_squared, 0.90);
}

TEST(MathTest, CurveFit2DGaussian) {
  auto model = [](std::vector<std::vector<double>> x,
                  std::vector<std::vector<double>> y,
                  std::vector<std::vector<double>> p) {
    double A = p[0][0];
    double x0 = p[0][1];
    double y0 = p[0][2];
    double sigma = p[0][3];
    std::vector<std::vector<double>> z(x.size(),
                                       std::vector<double>(x[0].size()));
    for (size_t i = 0; i < x.size(); ++i) {
      for (size_t j = 0; j < x[i].size(); ++j) {
        double dx = x[i][j] - x0;
        double dy = y[i][j] - y0;
        z[i][j] = A * std::exp(-(dx * dx + dy * dy) / (2 * sigma * sigma));
      }
    }
    return z;
  };

  std::vector<std::vector<double>> x(10, std::vector<double>(10));
  std::vector<std::vector<double>> y(10, std::vector<double>(10));
  for (int i = 0; i < 10; ++i) {
    for (int j = 0; j < 10; ++j) {
      x[i][j] = i * 1.0;
      y[i][j] = j * 1.0;
    }
  }

  std::vector<double> true_p = {10.0, 5.0, 5.0, 2.0};
  std::vector<std::vector<double>> z = model(x, y, {true_p});

  fitting_parameters params;
  params.initial_guess = {9.0, 4.5, 4.5, 1.5};
  params.bounds = {{0.0, 20.0}, {0.0, 10.0}, {0.0, 10.0}, {0.1, 5.0}};

  // Set seed for reproducibility
  arma::arma_rng::set_seed(42);

  auto res = curvefit2D(model, x, y, z, params);

  ASSERT_TRUE(res.success);
  EXPECT_NEAR(res.coefficients[0], true_p[0], 0.5);
  EXPECT_NEAR(res.coefficients[1], true_p[1], 0.5);
  EXPECT_NEAR(res.coefficients[2], true_p[2], 0.5);
  EXPECT_NEAR(res.coefficients[3], true_p[3], 0.5);
  EXPECT_GT(res.r_squared, 0.95);
}

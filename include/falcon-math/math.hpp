#pragma once
#include "falcon-math/export.h"
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace falcon {
namespace math {

// The first vector is the domain (x values) and the second vector is the
// parameters of the model
using func1D = std::function<std::vector<double>(std::vector<double>,
                                                 std::vector<double>)>;
using func2D = std::function<std::vector<std::vector<double>>(
    std::vector<std::vector<double>>, std::vector<std::vector<double>>,
    std::vector<std::vector<double>>)>;

struct FALCON_MATH_API fitting_parameters {
  std::optional<std::vector<double>> initial_guess;
  std::optional<std::vector<std::pair<double, double>>> bounds;
  std::optional<std::string> method;
};

struct FALCON_MATH_API curvefit_result {
  std::vector<double> coefficients;
  std::string error_message;
  bool success;
  double r_squared;
};

enum class AnalysisType {
  TURN_ON,
  PINCH_OFF,
  CHANNEL_ACCUMULATION_2D
};

struct FALCON_MATH_API AnalysisResult {
  AnalysisType type;
  curvefit_result fit;
  std::map<std::string, double> extracted_parameters;
};

/**
 * @brief Performs a curve fit on the given data.
 */
FALCON_MATH_API curvefit_result curvefit1D(func1D model,
                                           const std::vector<double>& x,
                                           const std::vector<double>& y,
                                           fitting_parameters params = {});

/**
 * @brief Performs a 2D curve fit on the given data.
 */
FALCON_MATH_API curvefit_result curvefit2D(func2D model,
                                           const std::vector<std::vector<double>>& x,
                                           const std::vector<std::vector<double>>& y,
                                           const std::vector<std::vector<double>>& z,
                                           fitting_parameters params = {});

// Predefined models
FALCON_MATH_API double sigmoid(double x, double A, double x0, double k, double b);

FALCON_MATH_API double piecewise_linear(double x, double x0, double x1,
                                        double m1, double m2, double y0);

FALCON_MATH_API double channel_accumulation_2d(double x, double y,
                                               double cx, double cy, double cm, double cr,
                                               double m1, double m2, double m3,
                                               double bx, double dx, double dy,
                                               double dm, double dr);

} // namespace math
// Backward compatibility alias: allow falcon::routine::...
namespace routine {
  using func1D = ::falcon::math::func1D;
  using func2D = ::falcon::math::func2D;
  using fitting_parameters = ::falcon::math::fitting_parameters;
  using curvefit_result = ::falcon::math::curvefit_result;
  using AnalysisType = ::falcon::math::AnalysisType;
  using AnalysisResult = ::falcon::math::AnalysisResult;
  using ::falcon::math::curvefit1D;
  using ::falcon::math::curvefit2D;
  using ::falcon::math::sigmoid;
  using ::falcon::math::piecewise_linear;
  using ::falcon::math::channel_accumulation_2d;
} // namespace routine
} // namespace falcon

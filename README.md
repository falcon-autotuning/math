# math

Mathematical models, nonlinear curve fitting routines, and FFI bindings for Falcon autotuning and quantum transport analysis.

## Overview

`math` provides high-performance C++ mathematical models and robust nonlinear curve fitting routines (powered by `ensmallen` and `armadillo`), exposed directly to the Falcon DSL:

- **Sigmoid Model**: Direct evaluation and 1D fitting with initial parameter estimation and bounds.
- **Piecewise Linear Model**: Direct evaluation and 1D fitting for transition/regime identification.
- **2D Channel Accumulation Model**: Direct evaluation and Differential Evolution (DE) optimization for 2D transport maps.

## Falcon DSL Interface (`math.fal`)

```fal
// Analytical evaluations
routine Sigmoid (float x, float A, float x0, float k, float b) -> (float y)
routine PiecewiseLinear (float x, float x0, float x1, float m1, float m2, float y0) -> (float y)
routine ChannelAccumulation2D (...) -> (float z)

// Math fitting routines
struct MathUtils {
    routine FitSigmoid (array::Array<float> x, array::Array<float> y) -> (float A, float x0, float k, float b, float r_squared, bool success)
    routine FitSigmoidWithGuess (array::Array<float> x, array::Array<float> y, float guess_A, float guess_x0, float guess_k, float guess_b) -> (float A, float x0, float k, float b, float r_squared, bool success)
    routine FitPiecewiseLinear (array::Array<float> x, array::Array<float> y) -> (float x0, float x1, float m1, float m2, float y0, float r_squared, bool success)
    routine FitPiecewiseLinearWithGuess (array::Array<float> x, array::Array<float> y, float guess_x0, float guess_x1, float guess_m1, float guess_m2, float guess_y0) -> (float x0, float x1, float m1, float m2, float y0, float r_squared, bool success)
    routine FitChannelAccumulation2D (array::Array<array::Array<float>> x, array::Array<array::Array<float>> y, array::Array<array::Array<float>> z) -> (array::Array<float> params, float r_squared, bool success)
}
```

## Building & Testing

```bash
# Build C++ library and FFI wrapper
make build

# Run C++ unit tests
make test-cpp

# Run Falcon DSL tests
make test-fal

# Run all tests
make test

# Package release tarball
make archive
```

## License

MPL-2.0

# Pure C++ Ellipse Detector

This is a pure C++ implementation of the LU method for ellipse detection, converted from the original MATLAB implementation.

## Overview

The ellipse detector uses the LU method (Lu & Xie, 2020) to detect ellipses in images. This implementation is a complete C++ port that removes all MATLAB/MEX dependencies.

## Features

- Pure C++ implementation (no MATLAB required)
- Uses OpenCV for image I/O and visualization
- Efficient ellipse detection using arc-support line segments
- Multiple refinement methods (least-squares and Ahn's orthogonal distance fitting)

## Dependencies

- OpenCV (>= 3.0)
- Eigen3 (>= 3.0)
- CMake (>= 3.10)
- C++11 compatible compiler

## Building

### Ubuntu/Linux

```bash
# Install dependencies
sudo apt-get install libopencv-dev libeigen3-dev cmake

# Build
mkdir build
cd build
cmake ..
make
```

### macOS

```bash
# Install dependencies using Homebrew
brew install opencv eigen cmake

# Build
mkdir build
cd build
cmake ..
make
```

### Windows

1. Install OpenCV and Eigen3
2. Update CMakeLists.txt with correct paths if needed
3. Use CMake GUI or command line:

```bash
mkdir build
cd build
cmake ..
cmake --build . --config Release
```

## Usage

```bash
# Basic usage
./ellipse_detector <image_path>

# With parameters
./ellipse_detector <image_path> [Tac] [Tr] [polarity]

# Example
./ellipse_detector test_image.jpg 165 0.6 0
```

### Parameters

- `Tac` (default: 165): Threshold of elliptic angular coverage (0-360 degrees). Higher values require more complete ellipses.
- `Tr` (default: 0.6): Ratio of support inliers to ellipse (0-1). Higher values require more sufficient support.
- `polarity` (default: 0): 
  - 1: Detect ellipses with positive polarity only
  - -1: Detect ellipses with negative polarity only
  - 0: Detect all ellipses

## Converted MATLAB Functions

The following MATLAB functions have been converted to C++:

- `LCS_ellipse.m` → `src/main.cpp`
- `ellipseDetectionLU.m` → `src/ellipse_detection_lu.cpp`
- `ellipseDetection.m` → `src/ellipse_detection_main.cpp`
- `subEllipseDetection.m` → `src/sub_ellipse_detection.cpp`
- `fitEllipse.m` → `src/fit_ellipse.cpp`
- `fitAhn.m` → `src/fit_ahn.cpp`
- `computePointAngle.m` → `src/ellipse_utils.cpp`
- `dRosin_square.m` → `src/ellipse_utils.cpp`
- `takeInliers.m` → `src/ellipse_utils.cpp`
- `calcuCompleteness.m` → `src/ellipse_utils.cpp`
- `drawEllipses.m` → `src/draw_ellipses.cpp`
- `generateEllipseCandidates` (MEX) → `src/generate_candidates_standalone.cpp`

## Output

The program will:
1. Read the input image
2. Generate ellipse candidates
3. Detect ellipses using the LU method
4. Refine ellipses using Ahn's method
5. Display the detected ellipses on the image
6. Print ellipse parameters to console

Ellipse parameters are output in the format: (x, y, a, b, φ)
- (x, y): Center coordinates
- a: Semi-major axis
- b: Semi-minor axis  
- φ: Orientation angle in degrees

## Testing

A comprehensive test suite is available to compare results between MATLAB and C++ implementations:

```bash
# Run automated test suite
./tests/run_tests.sh

# Or run individual tests
./build/create_test_image
./build/test_cpp_detection
python3 tests/compare_results.py
```

See [tests/README.md](tests/README.md) for detailed testing documentation.

## Documentation

- [论文与实验重排方案（中文）](PAPER_EXPERIMENT_REORGANIZATION.md)

## References

Lu, C., Xia, S., Shao, M., & Fu, Y. (2019). Arc-support line segments revisited: An efficient high-quality ellipse detection. IEEE Transactions on Image Processing, 29, 768-781.

## License

This code is based on the original implementation by the authors. Please refer to their repository for licensing information.

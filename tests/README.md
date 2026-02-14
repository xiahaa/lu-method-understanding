# Ellipse Detection Test Suite

This directory contains test cases to compare ellipse detection results between the MATLAB and C++ implementations of the LU method.

## Overview

The test suite provides:
1. **Test image generation** - Creates synthetic images with known ellipses
2. **MATLAB test script** - Runs detection using the MATLAB implementation
3. **C++ test program** - Runs detection using the C++ implementation
4. **Comparison tool** - Compares and analyzes the results from both implementations

## Directory Structure

```
tests/
├── README.md                     # This file
├── run_tests.sh                  # Main test runner script
├── create_test_image.cpp         # C++ program to create test images
├── test_matlab_detection.m       # MATLAB test script
├── test_cpp_detection.cpp        # C++ test program
├── compare_results.py            # Python script to compare results
└── test_data/                    # Test data directory (created during test)
    ├── test_ellipses.png         # Test image with known ellipses
    ├── matlab_results.csv        # MATLAB detection results
    ├── cpp_results.csv           # C++ detection results
    ├── matlab_detection_result.png  # Visualization of MATLAB results
    └── cpp_detection_result.png     # Visualization of C++ results
```

## Prerequisites

### For C++ tests:
- CMake >= 3.10
- OpenCV >= 3.0
- Eigen3 >= 3.0
- LAPACK/LAPACKE
- C++11 compatible compiler

### For MATLAB tests:
- MATLAB (any recent version)
- The MATLAB implementation files in the parent directory

### For comparison:
- Python 3.x

## Building the Tests

The test programs are built as part of the main project:

```bash
cd /path/to/lu-method-understanding
mkdir -p build
cd build
cmake ..
make
```

This will create:
- `build/create_test_image` - Test image generator
- `build/test_cpp_detection` - C++ test program

## Running the Tests

### Automated Test Suite

The easiest way to run all tests is using the automated script:

```bash
./tests/run_tests.sh
```

This script will:
1. Create the test image
2. Run C++ detection
3. Run MATLAB detection (if MATLAB is available)
4. Compare the results

### Manual Testing

#### Step 1: Create Test Image

```bash
./build/create_test_image
```

This creates `tests/test_data/test_ellipses.png` with three known ellipses:
- Ellipse 1: center=(150, 150), a=60, b=40, angle=30°
- Ellipse 2: center=(350, 150), a=50, b=30, angle=0°
- Ellipse 3: center=(250, 350), a=70, b=50, angle=45°

#### Step 2: Run C++ Test

```bash
./build/test_cpp_detection
```

This will:
- Detect ellipses in the test image
- Save results to `tests/test_data/cpp_results.csv`
- Save visualization to `tests/test_data/cpp_detection_result.png`

#### Step 3: Run MATLAB Test

In MATLAB:
```matlab
cd /path/to/lu-method-understanding/tests
test_matlab_detection
```

This will:
- Detect ellipses in the test image
- Save results to `tests/test_data/matlab_results.csv`
- Save visualization to `tests/test_data/matlab_detection_result.png`

#### Step 4: Compare Results

```bash
python3 tests/compare_results.py
```

This will:
- Load both CSV result files
- Match corresponding ellipses
- Report differences in parameters
- Determine if the implementations are consistent

## Test Output Format

Both MATLAB and C++ tests save results in CSV format with the following columns:

```
x,y,a,b,phi_degrees
```

Where:
- `x, y` - Ellipse center coordinates
- `a` - Semi-major axis length
- `b` - Semi-minor axis length
- `phi_degrees` - Orientation angle in degrees

## Comparison Criteria

The comparison script matches ellipses based on:
- **Position similarity** - Center coordinates should be close
- **Axes similarity** - Semi-major and semi-minor axes should match
- **Angle similarity** - Orientation should be consistent (accounting for 180° symmetry)

A successful test means:
- Same number of ellipses detected
- All ellipses match between implementations
- Parameter differences are within acceptable tolerances

## Testing with Custom Images

You can test with your own images:

### C++ Test:
```bash
./build/test_cpp_detection /path/to/your/image.jpg /path/to/output.csv
```

### MATLAB Test:
Modify `test_matlab_detection.m` to change the `test_image` variable.

### Comparison:
```bash
python3 tests/compare_results.py /path/to/matlab.csv /path/to/cpp.csv
```

## Expected Results

For the default test image with three ellipses:
- Both implementations should detect all three ellipses
- Parameters should match within a few pixels/degrees
- Typical differences:
  - Position: < 1 pixel
  - Axes: < 2 pixels
  - Angle: < 5 degrees

## Troubleshooting

### "Test image not found" error
- Make sure you've run `create_test_image` first
- Check that `tests/test_data/` directory exists

### "MATLAB test failed" error
- Ensure all MATLAB `.m` files are in the parent directory
- Check that the image path in the test script is correct
- Make sure MATLAB can find all required functions

### "Results don't match" warning
- Small differences are normal due to floating-point precision
- Check the visualizations to see if detections look reasonable
- Large differences may indicate a bug in one implementation

### Build errors
- Ensure all dependencies are installed
- Check that OpenCV and Eigen3 are properly configured
- See main README.md for detailed build instructions

## Adding More Tests

To add additional test cases:

1. Create more test images in `test_data/`
2. Modify test scripts to process multiple images
3. Update comparison script to handle multiple test cases
4. Add assertions for expected results

Example test images to add:
- Partial ellipses (incomplete arcs)
- Multiple overlapping ellipses
- Ellipses with different sizes and orientations
- Real-world images from datasets

## Continuous Integration

For CI/CD pipelines, you can run the tests non-interactively:

```bash
./tests/run_tests.sh
echo $?  # 0 for success, 1 for failure
```

The comparison script returns:
- Exit code 0: Tests passed
- Exit code 1: Tests failed

## References

- Original paper: Lu, C., et al. (2019). "Arc-support line segments revisited: An efficient high-quality ellipse detection." IEEE TIP.
- MATLAB implementation: This repository's `.m` files
- C++ implementation: This repository's `src/` and `include/` files

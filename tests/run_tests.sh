#!/bin/bash
# Run all tests to compare MATLAB and C++ ellipse detection

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

echo "====================================================================="
echo "Ellipse Detection Test Suite"
echo "Comparing MATLAB and C++ implementations"
echo "====================================================================="
echo ""

# Step 1: Create test image
echo "Step 1: Creating test image..."
if [ ! -f "build/create_test_image" ]; then
    echo "Error: build/create_test_image not found. Please build the project first."
    echo "Run: cd build && cmake .. && make"
    exit 1
fi

./build/create_test_image
echo ""

# Step 2: Run C++ test
echo "Step 2: Running C++ ellipse detection test..."
if [ ! -f "build/test_cpp_detection" ]; then
    echo "Error: build/test_cpp_detection not found. Please build the project first."
    exit 1
fi

./build/test_cpp_detection
echo ""

# Step 3: Run MATLAB test (if MATLAB is available)
echo "Step 3: Running MATLAB ellipse detection test..."
if command -v matlab &> /dev/null; then
    cd tests
    matlab -nodisplay -nosplash -nodesktop -r "try; test_matlab_detection; catch ME; disp('Error occurred:'); disp(getReport(ME)); exit(1); end; exit(0);"
    cd ..
    echo ""
    
    # Step 4: Compare results
    echo "Step 4: Comparing results..."
    python3 tests/compare_results.py
    echo ""
else
    echo "MATLAB not found. Skipping MATLAB test."
    echo "Please run the MATLAB test manually:"
    echo "  1. Open MATLAB"
    echo "  2. cd to the repository directory"
    echo "  3. Run: cd tests; test_matlab_detection"
    echo "  4. Run comparison: python3 tests/compare_results.py"
    echo ""
fi

echo "====================================================================="
echo "Test suite completed!"
echo "====================================================================="

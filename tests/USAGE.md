# Example: Testing with Custom Images

This example demonstrates how to run ellipse detection tests with your own images.

## Quick Start

### 1. Test with the Default Synthetic Image

```bash
# Build the project
cd /path/to/lu-method-understanding
mkdir -p build && cd build
cmake ..
make

# Run the test
cd ..
./tests/run_tests.sh
```

### 2. Test with Your Own Image

#### Using C++ Test Program:

```bash
# Detect ellipses in your image
./build/test_cpp_detection /path/to/your/image.jpg tests/test_data/your_cpp_results.csv

# Results will be saved to: tests/test_data/your_cpp_results.csv
```

#### Using MATLAB Test Script:

In MATLAB:
```matlab
% Edit the test script to use your image
cd /path/to/lu-method-understanding/tests
% Modify test_matlab_detection.m to change test_image variable
test_matlab_detection
```

#### Comparing Results:

```bash
# Compare MATLAB and C++ results
python3 tests/compare_results.py \
    tests/test_data/matlab_results.csv \
    tests/test_data/cpp_results.csv
```

## Testing Parameters

Both MATLAB and C++ implementations use the same default parameters:
- `Tac = 165`: Elliptic angular coverage threshold (0-360 degrees)
- `Tr = 0.6`: Ratio of support inliers (0-1)
- `polarity = 0`: Detect all ellipses (1=positive, -1=negative, 0=all)

### Adjusting Parameters

For images with different characteristics, you may need to adjust parameters:

**For partial ellipses:**
```bash
# Lower Tac to detect incomplete ellipses
./build/ellipse_detector your_image.jpg 120 0.6 0
```

**For noisy images:**
```bash
# Increase Tr to require more support
./build/ellipse_detector your_image.jpg 165 0.8 0
```

**For small ellipses:**
```bash
# Lower both thresholds
./build/ellipse_detector your_image.jpg 100 0.4 0
```

## Image Requirements

For best results, your test images should:

1. **Have clear edges**: The algorithm detects ellipse boundaries, not filled regions
2. **Sufficient contrast**: Dark ellipses on light background (or vice versa)
3. **Reasonable size**: Ellipses should be at least 20-30 pixels in diameter
4. **Complete or nearly complete**: Very partial ellipses may not be detected

## Example Test Images

### Create Simple Test Images

You can create your own test images using:

**Python (with OpenCV):**
```python
import cv2
import numpy as np

img = np.ones((500, 500, 3), dtype=np.uint8) * 255
cv2.ellipse(img, (250, 250), (100, 60), 30, 0, 360, (0, 0, 0), 5)
cv2.imwrite('test_ellipse.png', img)
```

**C++:**
```cpp
Mat img = Mat::ones(500, 500, CV_8UC3) * 255;
ellipse(img, Point(250, 250), Size(100, 60), 30, 0, 360, Scalar(0, 0, 0), 5);
imwrite("test_ellipse.png", img);
```

### Download Real Test Images

You can also use real-world datasets:
- [AAMED Ellipse Datasets](https://github.com/AlanLuSun/AAMED)
- [Prasad Dataset](https://www.dropbox.com/sh/y7z5l30kp7yqmzo/AABZi6TKfBKlUmZDj9c9c9bna)

## Interpreting Results

### Comparison Output

The comparison script reports:
- **Number of ellipses detected** by each implementation
- **Matched ellipses** and their parameter differences
- **Unmatched ellipses** (detected by only one implementation)

### Expected Differences

Small differences are normal:
- **Position**: < 1-2 pixels
- **Axes**: < 2-3 pixels  
- **Angle**: < 5 degrees

Larger differences may indicate:
- Different parameter settings
- Numerical precision differences
- Implementation bugs

### Success Criteria

Tests pass when:
- ✓ Same number of ellipses detected
- ✓ All ellipses match between implementations
- ✓ Parameter differences are within tolerance

## Troubleshooting

### No Ellipses Detected

If neither implementation detects ellipses:
1. Check image has clear edges (use edge detection to verify)
2. Try lowering Tac parameter
3. Try lowering Tr parameter
4. Ensure ellipses are large enough (>20 pixels diameter)
5. Check image format and loading

### Different Number of Ellipses

If implementations detect different numbers:
1. Check parameter settings are identical
2. Look at visualizations to see which is correct
3. May indicate edge cases or bugs

### Large Parameter Differences

If matched ellipses have large differences:
1. Verify both are using same refinement (Ahn's method)
2. Check for numerical stability issues
3. May indicate implementation differences

## Advanced Testing

### Batch Testing

Test multiple images:

```bash
#!/bin/bash
for img in /path/to/images/*.jpg; do
    echo "Testing $img"
    ./build/test_cpp_detection "$img" "results/$(basename $img).csv"
done
```

### Automated Validation

Create a validation script:

```python
import sys
from compare_results import compare_ellipses

test_cases = [
    ('image1.jpg', 'matlab1.csv', 'cpp1.csv'),
    ('image2.jpg', 'matlab2.csv', 'cpp2.csv'),
]

all_passed = True
for name, matlab_csv, cpp_csv in test_cases:
    print(f"\nTesting {name}...")
    if not compare_ellipses(matlab_csv, cpp_csv, verbose=False):
        all_passed = False
        print(f"FAILED: {name}")

sys.exit(0 if all_passed else 1)
```

### Ground Truth Validation

For synthetic images with known ellipses, validate detection accuracy:

```python
def validate_against_ground_truth(detected_csv, ground_truth):
    """
    Compare detected ellipses against known ground truth.
    
    Args:
        detected_csv: Path to CSV with detected ellipses
        ground_truth: List of known ellipse parameters
                      [(x, y, a, b, phi), ...]
    """
    detected = read_ellipses(detected_csv)
    
    for gt in ground_truth:
        # Find closest match
        # Check parameters within tolerance
        # Report accuracy
    
    return accuracy_metrics
```

## Integration with CI/CD

Add to your CI pipeline:

```yaml
# .github/workflows/test.yml
- name: Build tests
  run: |
    mkdir build && cd build
    cmake ..
    make

- name: Run C++ tests
  run: |
    ./build/create_test_image
    ./build/test_cpp_detection

- name: Compare results (if MATLAB available)
  run: |
    python3 tests/compare_results.py
```

## References

- Main README: [../README.md](../README.md)
- Test Suite README: [tests/README.md](tests/README.md)
- Algorithm: Lu, C., et al. (2019). "Arc-support line segments revisited: An efficient high-quality ellipse detection." IEEE TIP.

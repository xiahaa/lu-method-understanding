#!/usr/bin/env python3
"""
Validation test for the comparison script.

This test verifies that the comparison functionality works correctly
by testing with known mock data.
"""

import sys
import os
import tempfile

# Add the tests directory to path
sys.path.insert(0, os.path.dirname(__file__))

from compare_results import compare_ellipses, read_ellipses, ellipse_distance, match_ellipses

def test_read_ellipses():
    """Test reading ellipses from CSV."""
    print("Testing CSV reading...")
    
    # Create a temporary CSV file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("x,y,a,b,phi_degrees\n")
        f.write("100.5,200.3,50.2,30.1,45.0\n")
        f.write("300.0,400.0,60.0,40.0,90.0\n")
        temp_file = f.name
    
    try:
        ellipses = read_ellipses(temp_file)
        assert ellipses is not None, "Failed to read ellipses"
        assert len(ellipses) == 2, f"Expected 2 ellipses, got {len(ellipses)}"
        assert abs(ellipses[0]['x'] - 100.5) < 0.01, "X coordinate mismatch"
        assert abs(ellipses[0]['a'] - 50.2) < 0.01, "Axis 'a' mismatch"
        print("✓ CSV reading test passed")
        return True
    finally:
        os.unlink(temp_file)

def test_ellipse_distance():
    """Test ellipse distance calculation."""
    print("Testing ellipse distance calculation...")
    
    e1 = {'x': 100, 'y': 200, 'a': 50, 'b': 30, 'phi': 45}
    e2 = {'x': 101, 'y': 201, 'a': 51, 'b': 31, 'phi': 46}
    
    dist = ellipse_distance(e1, e2)
    
    # Check that distances are reasonable
    assert 0 < dist['position'] < 2, f"Unexpected position distance: {dist['position']}"
    assert 0 < dist['axes'] < 2, f"Unexpected axes distance: {dist['axes']}"
    assert 0 < dist['angle'] < 2, f"Unexpected angle distance: {dist['angle']}"
    
    print("✓ Ellipse distance calculation test passed")
    return True

def test_angle_symmetry():
    """Test that angle comparison handles 180-degree symmetry."""
    print("Testing angle symmetry...")
    
    from compare_results import angle_difference
    
    # Test basic difference
    assert abs(angle_difference(45, 50) - 5) < 0.1, "Basic angle difference failed"
    
    # Test 180-degree symmetry (ellipses are symmetric)
    assert abs(angle_difference(45, 225) - 0) < 0.1, "180-degree symmetry failed"
    assert abs(angle_difference(0, 180) - 0) < 0.1, "0/180 symmetry failed"
    
    # Test wrap-around
    assert abs(angle_difference(5, 355) - 10) < 0.1, "Wrap-around failed"
    
    print("✓ Angle symmetry test passed")
    return True

def test_matching():
    """Test ellipse matching algorithm."""
    print("Testing ellipse matching...")
    
    matlab_ellipses = [
        {'x': 100, 'y': 200, 'a': 50, 'b': 30, 'phi': 45},
        {'x': 300, 'y': 400, 'a': 60, 'b': 40, 'phi': 90},
    ]
    
    cpp_ellipses = [
        {'x': 100.5, 'y': 200.5, 'a': 50.2, 'b': 30.1, 'phi': 45.2},
        {'x': 300.3, 'y': 400.2, 'a': 60.1, 'b': 40.0, 'phi': 90.1},
    ]
    
    matches, unmatched_m, unmatched_c = match_ellipses(matlab_ellipses, cpp_ellipses)
    
    assert len(matches) == 2, f"Expected 2 matches, got {len(matches)}"
    assert len(unmatched_m) == 0, f"Expected 0 unmatched MATLAB, got {len(unmatched_m)}"
    assert len(unmatched_c) == 0, f"Expected 0 unmatched C++, got {len(unmatched_c)}"
    
    print("✓ Ellipse matching test passed")
    return True

def test_comparison_integration():
    """Test full comparison with temporary files."""
    print("Testing full comparison integration...")
    
    # Create temporary CSV files
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("x,y,a,b,phi_degrees\n")
        f.write("100.0,200.0,50.0,30.0,45.0\n")
        matlab_file = f.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("x,y,a,b,phi_degrees\n")
        f.write("100.5,200.5,50.2,30.1,45.2\n")
        cpp_file = f.name
    
    try:
        # Should pass with small differences
        result = compare_ellipses(matlab_file, cpp_file, verbose=False)
        assert result == True, "Comparison should pass with small differences"
        print("✓ Full comparison integration test passed")
        return True
    finally:
        os.unlink(matlab_file)
        os.unlink(cpp_file)

def main():
    """Run all validation tests."""
    print("="*70)
    print("Validation Tests for Ellipse Detection Comparison")
    print("="*70)
    print()
    
    tests = [
        test_read_ellipses,
        test_ellipse_distance,
        test_angle_symmetry,
        test_matching,
        test_comparison_integration,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
            else:
                failed += 1
                print(f"✗ {test.__name__} failed")
        except Exception as e:
            failed += 1
            print(f"✗ {test.__name__} failed with exception: {e}")
            import traceback
            traceback.print_exc()
    
    print()
    print("="*70)
    print(f"Results: {passed}/{len(tests)} tests passed")
    print("="*70)
    
    if failed == 0:
        print("✓ All validation tests passed!")
        return 0
    else:
        print(f"✗ {failed} test(s) failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())

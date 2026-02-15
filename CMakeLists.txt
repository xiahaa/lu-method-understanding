cmake_minimum_required(VERSION 3.10)
project(EllipseDetector)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(OpenCV REQUIRED)
find_package(Eigen3 REQUIRED)

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${OpenCV_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
)

# Source files
set(ELLIPSE_SOURCES
    src/canny.cpp
    src/clustering.cpp
    src/ellipse_geometry.cpp
    src/elsd.cpp
    src/gauss.cpp
    src/gen_init_set.cpp
    src/gradient.cpp
    src/group_forming.cpp
    src/image.cpp
    src/nfa.cpp
    src/rect.cpp
    src/tuple.cpp
    src/utils.cpp
    src/ellipse_utils.cpp
    src/fit_ellipse.cpp
    src/fit_ahn.cpp
    src/draw_ellipses.cpp
    src/sub_ellipse_detection.cpp
    src/ellipse_detection_main.cpp
    src/ellipse_detection_lu.cpp
    src/generate_candidates_standalone.cpp
    src/logging_stats.cpp
)

# Executable
add_executable(ellipse_detector
    src/main.cpp
    ${ELLIPSE_SOURCES}
)

# Link libraries
target_link_libraries(ellipse_detector
    ${OpenCV_LIBS}
    Eigen3::Eigen
    lapack
    blas
)

# Install
install(TARGETS ellipse_detector DESTINATION bin)

# Test executables
add_executable(create_test_image
    tests/create_test_image.cpp
)

add_executable(test_cpp_detection
    tests/test_cpp_detection.cpp
    ${ELLIPSE_SOURCES}
)

# Link libraries for test executables
target_link_libraries(create_test_image
    ${OpenCV_LIBS}
)

target_link_libraries(test_cpp_detection
    ${OpenCV_LIBS}
    Eigen3::Eigen
    lapack
    blas
)

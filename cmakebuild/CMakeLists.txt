cmake_minimum_required(VERSION 3.12)
project(RDKitDescriptors)

# Set C++ Standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set policies for compatibility
cmake_policy(SET CMP0167 NEW)
cmake_policy(SET CMP0144 NEW)

# Find Python and pybind11
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED)

# Find RDKit (users must set RDKit_DIR manually or use CMake's cache)
find_package(RDKit REQUIRED RDKitGraphMol RDKitDescriptors RDKitFingerprints RDKitRDGeneral RDKitDataStructs)

# Include RDKit directories
include_directories(${RDKit_INCLUDE_DIRS})
link_directories("/Users/guillaume-osmo/miniconda3/pkgs/rdkit-2023.09.3-py311h64bc748_1/lib")
link_directories(${RDKit_LIBRARY_DIRS})

# Find Eigen3 (users should install Eigen3 in a discoverable path)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

# Find Boost (ensure Python support)
find_package(Boost REQUIRED COMPONENTS python)

# LAPACK configuration
find_package(LAPACK REQUIRED)

if (LAPACK_FOUND)
    include_directories("/opt/homebrew/Cellar/lapack/3.12.0/include")
    set(LAPACK_LIBRARIES_LIBS
        "/opt/homebrew/Cellar/lapack/3.12.0/lib/liblapacke.dylib"
        "/opt/homebrew/Cellar/lapack/3.12.0/lib/liblapack.dylib"
        "/opt/homebrew/Cellar/lapack/3.12.0/lib/libblas.dylib"
    )
    message(STATUS "LAPACK found: ${LAPACK_LIBRARIES}")
else()
    message(FATAL_ERROR "LAPACK not found. Install LAPACK using your package manager.")
endif()

# Create the pybind11 module
pybind11_add_module(cppmordred cppmordred/cppmordred.cpp)

# Link libraries
target_link_libraries(cppmordred PRIVATE
    RDKitGraphMol
    RDKitDescriptors
    RDKitFingerprints
    RDKitRDGeneral
    RDKitDataStructs
    ${LAPACK_LIBRARIES}
    ${LAPACK_LIBRARIES_LIBS}
    ${Boost_LIBRARIES}
    ${Boost_PYTHON_LIBRARY} 
)

# Include directories
target_include_directories(cppmordred PRIVATE
    ${Eigen3_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
)
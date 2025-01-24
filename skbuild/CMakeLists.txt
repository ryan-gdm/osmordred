cmake_minimum_required(VERSION 3.21)

project(cppmordred)

# Set policies
cmake_policy(SET CMP0167 NEW)
cmake_policy(SET CMP0144 NEW)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
#set(CMAKE_MACOSX_RPATH 1)

# define the C++ library "_cppmordred"
add_library(_cppmordred MODULE cppmordred/cppmordred.cpp)

find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
#execute_process(
#    COMMAND ${Python3_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
#    OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
#    OUTPUT_STRIP_TRAILING_WHITESPACE
#)


include_directories(${PYTHON_INCLUDE_DIR})



find_package(RDKit REQUIRED)
include_directories(${RDKit_INCLUDE_DIRS})
# RDKit isn't reporting its libraries, so set a path and use it to point to them directly below.
# This relies on the fact that we're installing it via Conda.
set(RDKit_HACK_LIBRARY_DIR "${Python3_ROOT_DIR}/lib")

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "${Python3_ROOT_DIR}/lib;/opt/homebrew/opt/lapack/lib")

# set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
# set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Find and link Eigen
find_path(EIGEN3_INCLUDE_DIR Eigen/Dense
    PATHS /opt/homebrew/include/eigen3 /usr/local/include/eigen3 /usr/include/eigen3
    REQUIRED
)
include_directories(${EIGEN3_INCLUDE_DIR})

# Find and link Boost
find_package(Boost REQUIRED python)
include_directories(${Boost_INCLUDE_DIRS})


# Find and link LAPACK
find_package(LAPACK REQUIRED)
if (LAPACK_FOUND)
    message(STATUS "LAPACK found: ${LAPACK_LIBRARIES}")
    set(LAPACK_INCLUDEDIR "${LAPACK_ROOT}/include")
    set(LAPACK_LIBRARYDIR "${LAPACK_ROOT}/lib")
    include_directories(${LAPACK_INCLUDE_DIRS})
    # Explicit LAPACK library paths
    set(LAPACK_LIBRARIES_LIBS
    "/opt/homebrew/opt/lapack/lib/liblapacke.dylib"
    "/opt/homebrew/opt/lapack/lib/liblapack.dylib"
    "/opt/homebrew/opt/lapack/lib/libblas.dylib"
    )
else()
    message(FATAL_ERROR "LAPACK not found. Install LAPACK via Homebrew or appropriate package manager.")
endif()


# Find and link RDKit libraries
target_link_libraries(_cppmordred PRIVATE
    ${RDKit_HACK_LIBRARY_DIR}/libRDKitDataStructs.dylib
    ${RDKit_HACK_LIBRARY_DIR}/libRDKitDescriptors.dylib
    ${RDKit_HACK_LIBRARY_DIR}/libRDKitFingerprints.dylib
    ${RDKit_HACK_LIBRARY_DIR}/libRDKitGraphMol.dylib
    ${RDKit_HACK_LIBRARY_DIR}/libRDKitPartialCharges.dylib
    ${RDKit_HACK_LIBRARY_DIR}/libRDKitRDGeneral.dylib
    ${RDKit_HACK_LIBRARY_DIR}/libRDKitSmilesParse.dylib
    ${RDKit_HACK_LIBRARY_DIR}/libRDKitSubgraphs.dylib
    ${RDKit_HACK_LIBRARY_DIR}/libRDKitSubstructMatch.dylib
    ${LAPACK_LIBRARIES_LIBS}
    ${Boost_LIBRARIES}
    ${Boost_PYTHON_LIBRARY}  # Ensure Boost.Python is linked
    ${Python3_LIBRARIES}
    )

# set(CMAKE_INSTALL_RPATH "${RDKit_HACK_LIBRARY_DIR}")
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set_target_properties(_cppmordred PROPERTIES
    PREFIX ""
    SUFFIX ".so"
)

install(TARGETS _cppmordred LIBRARY DESTINATION cppmordred)
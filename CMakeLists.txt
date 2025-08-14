cmake_minimum_required(VERSION 3.28)
project(osmordred LANGUAGES CXX)

# Prefer conda prefix for discovery
if(DEFINED ENV{CONDA_PREFIX})
  list(PREPEND CMAKE_PREFIX_PATH "$ENV{CONDA_PREFIX}")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(_osmordred MODULE osmordred/osmordred.cpp)
set_target_properties(_osmordred PROPERTIES PREFIX "" SUFFIX ".so" POSITION_INDEPENDENT_CODE ON)

# ---------- Python ----------
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
target_link_libraries(_osmordred PRIVATE Python3::Module)
if(TARGET Python3::Python)
  target_link_libraries(_osmordred PRIVATE Python3::Python)
endif()
if(DEFINED Python3_INCLUDE_DIRS)
  target_include_directories(_osmordred PRIVATE ${Python3_INCLUDE_DIRS})
endif()
string(REGEX REPLACE "^([0-9]+)\\.([0-9]+).*" "\\1\\2" PYXY "${Python3_VERSION}")

# ---------- Eigen ----------
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(_osmordred PRIVATE Eigen3::Eigen)

# ---------- Boost (non-Python) ----------
find_package(Boost 1.86 REQUIRED COMPONENTS headers system iostreams serialization)
target_link_libraries(_osmordred PRIVATE Boost::system Boost::iostreams Boost::serialization)

# ---------- Boost.Python / NumPy (robust: find actual .so names) ----------
set(_bpy_candidates "boost_python${PYXY}" "boost_numpy${PYXY}" "boost_python3" "boost_numpy3" "boost_python" "boost_numpy")
foreach(kind IN ITEMS python numpy)
  set(_names "")
  foreach(n IN LISTS _bpy_candidates)
    if(n MATCHES "boost_${kind}.*")
      list(APPEND _names "${n}")
    endif()
  endforeach()
  if(DEFINED ENV{CONDA_PREFIX})
    find_library(BOOST_${kind}_LIB NAMES ${_names} PATHS "$ENV{CONDA_PREFIX}/lib" NO_DEFAULT_PATH)
  endif()
  if(NOT BOOST_${kind}_LIB)
    find_library(BOOST_${kind}_LIB NAMES ${_names})
  endif()
  if(NOT BOOST_${kind}_LIB)
    message(FATAL_ERROR "Could not find libboost_${kind} for Python ${Python3_VERSION}. Try: conda install -c conda-forge boost")
  endif()
endforeach()
target_link_libraries(_osmordred PRIVATE ${BOOST_python_LIB} ${BOOST_numpy_LIB})

# ---------- RDKit: try config/module; else manual headers+libs ----------
find_package(RDKit QUIET)

if(TARGET RDKit::RDKit)
  target_link_libraries(_osmordred PRIVATE RDKit::RDKit)
else()
  # If module didn’t give variables, try to locate manually
  if(NOT RDKit_INCLUDE_DIRS)
    # Look for a directory that contains 'rdkit/GraphMol/ROMol.h'
    find_path(RDKIT_INCLUDE_BASE
      NAMES rdkit/GraphMol/ROMol.h
      PATHS "$ENV{CONDA_PREFIX}/include" "/usr/local/include" "/usr/include"
    )
    if(NOT RDKIT_INCLUDE_BASE)
      message(FATAL_ERROR "RDKit headers not found. Install RDKit in this env, e.g.: conda install -c conda-forge rdkit")
    endif()
    set(RDKit_INCLUDE_DIRS "${RDKIT_INCLUDE_BASE}")
  endif()
  target_include_directories(_osmordred PRIVATE ${RDKit_INCLUDE_DIRS})

  # Core libs typically needed by descriptor/fingerprint code
  set(_rdkit_core_libs
    RDKitGraphMol
    RDKitRDGeneral
    RDKitDataStructs
    RDKitSmilesParse
    RDKitSubstructMatch
    RDKitSubgraphs
    RDKitFingerprints
    RDKitDescriptors
    RDKitPartialCharges
  )

  set(RDKIT_MANUAL_LIBS "")
  foreach(lib ${_rdkit_core_libs})
    if(DEFINED ENV{CONDA_PREFIX})
      find_library(${lib}_PATH NAMES ${lib} PATHS "$ENV{CONDA_PREFIX}/lib" NO_DEFAULT_PATH)
    endif()
    if(NOT ${lib}_PATH)
      find_library(${lib}_PATH NAMES ${lib})
    endif()
    if(${lib}_PATH)
      list(APPEND RDKIT_MANUAL_LIBS ${${lib}_PATH})
    endif()
  endforeach()

  # Require at least the absolute minimum set:
  set(_min_required RDKitGraphMol RDKitRDGeneral RDKitDataStructs RDKitSmilesParse RDKitSubstructMatch)
  foreach(req IN LISTS _min_required)
    if(NOT DEFINED ${req}_PATH)
      message(FATAL_ERROR "Missing required RDKit library: ${req}. Check RDKit installation in this env (conda install -c conda-forge rdkit).")
    endif()
  endforeach()

  target_link_libraries(_osmordred PRIVATE ${RDKIT_MANUAL_LIBS})
endif()

# ---------- LAPACKE / LAPACK (remove if unused) ----------
find_path(LAPACKE_INCLUDE_DIR lapacke.h PATHS "$ENV{CONDA_PREFIX}/include" NO_DEFAULT_PATH)
find_library(LAPACKE_LIBRARY lapacke PATHS "$ENV{CONDA_PREFIX}/lib" NO_DEFAULT_PATH)
if(LAPACKE_INCLUDE_DIR AND LAPACKE_LIBRARY)
  target_include_directories(_osmordred PRIVATE ${LAPACKE_INCLUDE_DIR})
  target_link_libraries(_osmordred PRIVATE ${LAPACKE_LIBRARY})
else()
  find_package(LAPACK QUIET)
  if(LAPACK_FOUND)
    target_link_libraries(_osmordred PRIVATE LAPACK::LAPACK)
  endif()
endif()

# ---------- RPATH ----------
set_target_properties(_osmordred PROPERTIES
  BUILD_RPATH "$ORIGIN;$ORIGIN/..;$ENV{CONDA_PREFIX}/lib"
  INSTALL_RPATH "$ORIGIN;$ORIGIN/..;$ENV{CONDA_PREFIX}/lib"
)

install(TARGETS _osmordred LIBRARY DESTINATION osmordred)

rdkit_library(Osmordred Osmordred.cpp LINK_LIBRARIES DataStructs Descriptors Fingerprints GraphMol PartialCharges RDGeneral SmilesParse Subgraphs SubstructMatch)
target_compile_definitions(Osmordred PRIVATE RDKIT_OSMORDRED_BUILD)

rdkit_headers(Osmordred.h DEST Osmordred)

find_path(EIGEN3_INCLUDE_DIR Eigen/Dense
    PATH_SUFFIXES /eigen
    REQUIRED
)
include_directories(Osmordred ${EIGEN3_INCLUDE_DIR})

find_package(LAPACK REQUIRED)
include_directories(Osmordred ${LAPACK_INCLUDE_DIRS})
if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    # ${LAPACK_INCLUDE_DIRS} isn't working yet on linux so hack it:
    include_directories(Osmordred /usr/include/lapacke)
endif ()

link_libraries(Osmordred ${LAPACK_LIBRARIES})

# rdkit_test(testOsmordred testOsmordred.cpp LINK_LIBRARIES Osmordred )

if(RDK_BUILD_PYTHON_WRAPPERS)
add_subdirectory(Wrap)
endif()
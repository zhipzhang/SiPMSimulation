cmake_minimum_required(VERSION 3.10)
project(SiPMSimulation)

set(CMAKE_CXX_STANDARD 20)

# Add the source files

add_library(sipm_simulation SHARED
    ${CMAKE_CURRENT_SOURCE_DIR}/src/LACTSIPM.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/TELPIXEL.cpp
)
target_include_directories(sipm_simulation PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(sipm_simulation PRIVATE TMath TSpline TH1D TSpectrum TMultiGraph)
if(NOT TARGET Eigen3::Eigen)
    find_package(Eigen3 REQUIRED)
    if(NOT Eigen3_FOUND)
        message(FATAL_ERROR "Eigen3 not found")
    endif()
    target_link_libraries(sipm_simulation PRIVATE Eigen3::Eigen)
endif()
add_subdirectory(test)


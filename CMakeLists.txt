cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

set(Run "4.0" CACHE STRING "Gm2 data period")

project(ROOTEXAMPLE VERSION ${Run})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# ROOT
find_package(ROOT) # REQUIRED COMPONENTS Net)
list(APPEND EXTRA_LIBS ROOT::Net ROOT::Hist)

# Fitter
add_subdirectory(Fitter)
list(APPEND EXTRA_LIBS FITTER)

add_executable(MAIN test.cxx)
target_link_libraries(MAIN PUBLIC ${EXTRA_LIBS})
cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

set(VERSION_STR "4.0" CACHE STRING "Gm2 data period")
option(USE_JSON "Use json for initial values" ON)

project(ROOTEXAMPLE VERSION ${VERSION_STR})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# ROOT
find_package(ROOT) # REQUIRED COMPONENTS Net)
list(APPEND EXTRA_LIBS ROOT::Net ROOT::Hist)

# json
if(USE_JSON)
    list(APPEND EXTRA_INCLUDE_DIRS /usr/local/include)
    list(APPEND EXTRA_LINK_DIRS /usr/local/lib)
    list(APPEND EXTRA_LIBS jsoncpp)
endif()

# Fitter
add_subdirectory(Fitter)
list(APPEND EXTRA_LIBS FITTER)

add_executable(MAIN test.cxx)

target_include_directories(MAIN PUBLIC ${EXTRA_INCLUDE_DIRS})
target_link_directories(MAIN PUBLIC ${EXTRA_LINK_DIRS})
target_link_libraries(MAIN PUBLIC ${EXTRA_LIBS})
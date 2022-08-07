add_subdirectory(Blinders)

add_library(FITTER FitTools.cpp)

set(data_version_major ${PROJECT_VERSION_MAJOR})
set(data_version_minor ${PROJECT_VERSION_MINOR})

configure_file(configure_version.h.in configure_version.h)

target_link_libraries(FITTER PUBLIC BLINDERS ROOT::Net ROOT::Hist)
target_include_directories(FITTER PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(FITTER INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
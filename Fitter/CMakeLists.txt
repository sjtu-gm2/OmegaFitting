add_subdirectory(Blinders)

add_library(FITTER FitTools.cpp)

option(USE_JSON "Use json for initial values" OFF)
set(data_version_major ${PROJECT_VERSION_MAJOR})
set(data_version_minor ${PROJECT_VERSION_MINOR})
configure_file(configure.h.in configure.h)

list(APPEND EXTRA_LIBS BLINDERS ROOT::Net ROOT::Hist ROOT::MathCore ROOT::MathMore)
list(APPEND EXTRA_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR})

# json
if(USE_JSON)
    list(APPEND EXTRA_INCLUDE_DIRS /usr/local/include)
    list(APPEND EXTRA_LINK_DIRS /usr/local/lib)
    list(APPEND EXTRA_LIBS jsoncpp)
endif()

target_link_directories(FITTER PUBLIC ${EXTRA_LINK_DIRS})
target_link_libraries(FITTER PUBLIC ${EXTRA_LIBS})
target_include_directories(FITTER PUBLIC ${EXTRA_INCLUDE_DIRS})

target_include_directories(FITTER INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
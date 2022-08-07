# random lib
add_library(RLIB src/Random.cpp)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(RLIB INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Blinders
add_library(BLINDERS src/Blinders.cc)

target_include_directories(BLINDERS PUBLIC /usr/local/opt/openssl/include)
target_link_directories(BLINDERS PUBLIC /usr/local/opt/openssl/lib)
target_link_libraries(BLINDERS PUBLIC RLIB ssl crypto)
target_include_directories(BLINDERS INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
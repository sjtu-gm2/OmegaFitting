set(CMAKE_CXX_STANDARD 11)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(RLIB src/Random.cpp)

# blinder static
add_library(BLINDERS src/Blinders.cc)
target_include_directories(BLINDERS PUBLIC /usr/local/opt/openssl/include)
target_link_directories(BLINDERS PUBLIC /usr/local/opt/openssl/lib)
target_link_libraries(BLINDERS PUBLIC RLIB crypto)

# blinder shared
add_library(BLINDERS2 SHARED src/Blinders.cc)
target_include_directories(BLINDERS2 PUBLIC /usr/local/opt/openssl/include)
target_link_directories(BLINDERS2 PUBLIC /usr/local/opt/openssl/lib)
target_link_libraries(BLINDERS2 PUBLIC RLIB crypto)


# install blinder
# RLIB headers
file (GLOB HEADERSRLIB include/RandomLib/*.h*)

# Blinder headers
file (GLOB HEADERSBLINDER include/*.h*)

# install headers
install (FILES ${HEADERSRLIB} DESTINATION include/RandomLib)
install (FILES ${HEADERSBLINDER} DESTINATION include/)

# install libraries
install(TARGETS BLINDERS BLINDERS2 RLIB)
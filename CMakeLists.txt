project(PADVANTG)
cmake_minimum_required(VERSION 2.8)

# Location of the cnpy library
if (NOT CNPY_DIR)
  set(CNPY_DIR $ENV{HOME}/PADVANTG/build/cnpy)
endif ()
include_directories(${CNPY_DIR}/include)
link_directories(${CNPY_DIR}/lib)

# Build type -- Release or Debug
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif ()

# Automatically include linked directories in the rpath
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

# Suppress warning
list(APPEND CMAKE_CXX_FLAGS -Wno-narrowing)

# Print some variables to the screen
message("")
message(STATUS "CNPY_DIR:             ${CNPY_DIR}"            )
message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "CMAKE_BUILD_TYPE:     ${CMAKE_BUILD_TYPE}"    )
message(STATUS "CMAKE_CXX_FLAGS:      ${CMAKE_CXX_FLAGS}"     )
message("")

# Install executables
add_executable(calculate calculate.cpp)
target_link_libraries(calculate cnpy)
install(TARGETS calculate DESTINATION ".")

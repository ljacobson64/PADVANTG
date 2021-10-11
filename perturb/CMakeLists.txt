project(PADVANTG)
cmake_minimum_required(VERSION 2.8)

# Use C++17
set(CMAKE_CXX_STANDARD 17)

# Find OpenMP
find_package(OpenMP REQUIRED)

# Find Boost
find_package(Boost 1.65 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# Find HDF5
find_package(HDF5 COMPONENTS CXX REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})

# Default to Release build
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif ()

# Automatically include linked directories in the rpath
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

# Print some variables to the screen
message("")
message(STATUS "Boost_LIBRARIES:      ${Boost_LIBRARIES}"     )
message(STATUS "Boost_INCLUDE_DIR:    ${Boost_INCLUDE_DIR}"   )
message(STATUS "HDF5_LIBRARIES:       ${HDF5_LIBRARIES}"      )
message(STATUS "HDF5_INCLUDE_DIRS:    ${HDF5_INCLUDE_DIRS}"   )
message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "CMAKE_BUILD_TYPE:     ${CMAKE_BUILD_TYPE}"    )
message("")

# Install executables
add_executable(calculate main.cpp CalculateDR.cpp)
target_link_libraries(calculate OpenMP::OpenMP_CXX ${Boost_LIBRARIES} ${HDF5_LIBRARIES})
install(TARGETS calculate DESTINATION ".")

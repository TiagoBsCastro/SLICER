# Find required packages
find_package(GSL REQUIRED)
find_package(CFITSIO REQUIRED)
find_package(CCFITS REQUIRED)

# Include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${GSL_INCLUDE_DIRS}
    ${CFITSIO_INCLUDE_DIR}
    ${CCFITS_INCLUDE_DIR}
)

# Source files
set(SOURCES
    slicer-v2.cpp
    utilities.cpp
    w0waCDM.cpp
    densitymaps.cpp
    gadget2io.cpp
    writeplc.cpp
    data.cpp
)

# Add executable target
add_executable(SLICER ${SOURCES})

# Link libraries
target_link_libraries(SLICER
    ${GSL_LIBRARIES}
    ${CFITSIO_LIBRARIES}
    ${CCFITS_LIBRARY}
    -Wl,-rpath
    -lstdc++
    -lgsl
    -lgslcblas
    -lm
)

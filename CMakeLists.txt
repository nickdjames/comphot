# comphot top level CMakeLists
cmake_minimum_required (VERSION 2.8.11)
project (COMPHOT)

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

# Find required modules
include (gd)
include (cfitsio)
include (boost)

set (SRC_FILES
    src/app.cpp
    src/comphot.c
    src/proclib.c
)

# Output binary comphot
add_executable (comphot ${SRC_FILES})
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

# Set compiler requirements
target_compile_features(comphot PUBLIC
    cxx_nullptr
    cxx_nonstatic_member_init
)

# Set include directories
set (INCLUDE_PATHS ${INCLUDE_PATHS} ${CFITSIO_INCLUDE_PATH})
set (INCLUDE_PATHS ${INCLUDE_PATHS} ${GD_INCLUDE_PATH})

# Set library dependencies from included modules
set (LIBRARY_PATHS ${LIBRARY_PATHS} ${GD_LIBRARY_PATH})
set (LIBRARY_PATHS ${LIBRARY_PATHS} ${CFITSIO_LIBRARY_PATH})

# Handle platform specific configuration
if (WIN32)
    add_definitions (-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS)
else (WIN32)
    set_target_properties (comphot PROPERTIES COMPILE_OPTIONS "-Wall;-Wextra;-Wfloat-conversion")
    find_library (M_LIBRARY_PATH m)
    set(LIBRARY_PATHS ${LIBRARY_PATHS} ${M_LIBRARY_PATH})
endif (WIN32)

target_include_directories (comphot PUBLIC ${INCLUDE_PATHS})
target_link_libraries (comphot ${LIBRARY_PATHS})

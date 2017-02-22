# comphot top level CMakeLists
cmake_minimum_required (VERSION 2.8.11)
project (COMPHOT)

set (COMPHOT_VERSION_MAJOR 1)
set (COMPHOT_VERSION_MINOR 2)
set (COMPHOT_VERSION_PATCH 2)
set (COMPHOT_VERSION ${COMPHOT_VERSION_MAJOR}.${COMPHOT_VERSION_MINOR}.${COMPHOT_VERSION_PATCH})

set (CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/cmake/modules
)

set (COMPHOT_DESTINATION Comphot)
set (CMAKE_INSTALL_SYSTEM_RUNTIME_DESTINATION ${COMPHOT_DESTINATION})

# Find required modules
include (InstallRequiredSystemLibraries)
include (gd)
include (cfitsio)
include (boost)
include (m)

set (SRC_FILES
    ${CMAKE_SOURCE_DIR}/src/app.cpp
    ${CMAKE_SOURCE_DIR}/src/comphot.c
    ${CMAKE_SOURCE_DIR}/src/proclib.c
)

# Output binary comphot
add_executable (comphot ${SRC_FILES})
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

# Set compiler requirements
target_compile_features (comphot PUBLIC
    cxx_nullptr
    cxx_nonstatic_member_init
)

# Set common definitions
add_definitions (-DVERSION="${COMPHOT_VERSION}")

# Handle platform specific configuration
if (WIN32)
    # Set Windows only definitions
    add_definitions (-D_CRT_SECURE_NO_WARNINGS)
else (WIN32)
    # Set Unix compile options
    set_target_properties (comphot PROPERTIES COMPILE_OPTIONS -Wall -Wextra -Wfloat-conversion)
endif (WIN32)

target_include_directories (comphot PUBLIC ${INCLUDE_PATHS})
target_link_libraries (comphot ${LIBRARY_PATHS})

# Set installation locations
install (FILES ${RUNTIME_LIBRARY_PATHS} DESTINATION ${COMPHOT_DESTINATION})
install (TARGETS comphot RUNTIME DESTINATION ${COMPHOT_DESTINATION})

# CPack configuration for packaging comphot
set (CPACK_RESOURCE_FILE_LICENSE ${CMAKE_SOURCE_DIR}/LICENSE.txt)
set (CPACK_PACKAGE_NAME "Comphot")
set (CPACK_PACKAGE_VENDOR "BAA")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Comet photometry developed for the BAA comet section")
set (CPACK_PACKAGE_VERSION ${COMPHOT_VERSION})
set (CPACK_PACKAGE_VERSION_MAJOR ${COMPHOT_VERSION_MAJOR})
set (CPACK_PACKAGE_VERSION_MINOR ${COMPHOT_VERSION_MINOR})
set (CPACK_PACKAGE_VERSION_PATCH ${COMPHOT_VERSION_PATCH})
set (CPACK_PACKAGE_INSTALL_DIRECTORY ${COMPHOT_DESTINATION})
# Windows WiX specific configuration
set (CPACK_WIX_UPGRADE_GUID "C71583B4-CF14-487F-B27A-BAA86F9038DE")
set (CPACK_WIX_PRODUCT_ICON ${CMAKE_SOURCE_DIR}/images/baa.ico)
set (CPACK_WIX_UI_BANNER ${CMAKE_SOURCE_DIR}/images/baa_banner.bmp)
set (CPACK_WIX_UI_DIALOG ${CMAKE_SOURCE_DIR}/images/baa_welcome.bmp)
set (CPACK_WIX_PATCH_FILE ${CMAKE_SOURCE_DIR}/cmake/cpack_wix_patch.xml)

include (CPack)


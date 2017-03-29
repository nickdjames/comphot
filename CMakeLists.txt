# comphot top level CMakeLists
cmake_minimum_required (VERSION 3.1.3)
project (COMPHOT)

set (PROJECT_VERSION_MAJOR 1)
set (PROJECT_VERSION_MINOR 2)
set (PROJECT_VERSION_PATCH 5)
set (COMPHOT_VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH})
set (COMPHOT_DESTINATION Comphot)

# Find required modules
set (CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
include (gd)
include (cfitsio)
include (boost)
include (m)

# Handle platform specific configuration
if (WIN32)
    # Set Windows only definitions
    add_definitions (-D_CRT_SECURE_NO_WARNINGS)
else (WIN32)
    # Set Unix compile options
    add_compile_options (-Wall -Wextra -Wfloat-conversion)
endif (WIN32)

add_subdirectory (lib)
add_subdirectory (app)

# CPack configuration for packaging comphot
set (CPACK_RESOURCE_FILE_LICENSE ${PROJECT_SOURCE_DIR}/LICENSE.txt)
set (CPACK_PACKAGE_NAME "Comphot")
set (CPACK_PACKAGE_VENDOR "BAA")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Comet photometry developed for the BAA comet section")
set (CPACK_PACKAGE_VERSION ${COMPHOT_VERSION})
set (CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set (CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set (CPACK_PACKAGE_VERSION_PATCH ${PROJECT_VERSION_PATCH})
set (CPACK_PACKAGE_INSTALL_DIRECTORY ${COMPHOT_DESTINATION})
# Windows WiX specific configuration
set (CPACK_WIX_UPGRADE_GUID "C71583B4-CF14-487F-B27A-BAA86F9038DE")
set (CPACK_WIX_PRODUCT_ICON ${PROJECT_SOURCE_DIR}/images/baa.ico)
set (CPACK_WIX_UI_BANNER ${PROJECT_SOURCE_DIR}/images/baa_banner.bmp)
set (CPACK_WIX_UI_DIALOG ${PROJECT_SOURCE_DIR}/images/baa_welcome.bmp)
set (CPACK_WIX_PATCH_FILE ${PROJECT_SOURCE_DIR}/cmake/cpack_wix_patch.xml)

include (CPack)


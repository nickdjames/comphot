# app CMakeLists
set (CMAKE_INSTALL_SYSTEM_RUNTIME_DESTINATION ${COMPHOT_DESTINATION})
include (InstallRequiredSystemLibraries)

add_executable (app app.cpp)

set_target_properties (app PROPERTIES
    OUTPUT_NAME comphot
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
)

target_compile_features (app PUBLIC
    cxx_nullptr
    cxx_nonstatic_member_init
)

target_include_directories (app PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${INCLUDE_PATHS})
target_link_libraries (app LINK_PUBLIC comphot ${LIBRARY_PATHS})


# Set installation locations
install (FILES ${RUNTIME_LIBRARY_PATHS} DESTINATION ${COMPHOT_DESTINATION})
install (TARGETS app RUNTIME DESTINATION ${COMPHOT_DESTINATION})


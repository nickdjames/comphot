# comphot library CMakeLists

add_library (comphot SHARED comphot.c proclib.c)

set_target_properties (comphot PROPERTIES
    COMPILE_DEFINITIONS VERSION="${COMPHOT_VERSION}"
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
)

target_include_directories (comphot PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${INCLUDE_PATHS})
target_link_libraries (comphot LINK_PUBLIC ${LIBRARY_PATHS})

# Set installation locations
install (TARGETS comphot RUNTIME DESTINATION ${COMPHOT_DESTINATION})


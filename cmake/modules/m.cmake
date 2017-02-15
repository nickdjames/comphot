# Find additional libraries required on Unix

if (NOT WIN32)
    find_library (M_LIBRARY_PATH m)
    set (LIBRARY_PATHS ${LIBRARY_PATHS} ${M_LIBRARY_PATH})
else (NOT WIN32)
    add_definitions (-D_USE_MATH_DEFINES)
endif (NOT WIN32)


# Find gd include and lib paths

if (WIN32)
    set (GD_NAMES ${GD_NAMES} gd bgd libgd)
else (WIN32)
    set (GD_NAMES ${GD_NAMES} gd)
endif (WIN32)

find_path (GD_INCLUDE_PATH gd.h)
find_library (GD_LIBRARY_PATH NAMES ${GD_NAMES})

if (NOT GD_INCLUDE_PATH OR NOT GD_LIBRARY_PATH)
    message (FATAL_ERROR "Failed to find library gd")
endif (NOT GD_INCLUDE_PATH OR NOT GD_LIBRARY_PATH)

if (WIN32)
    get_filename_component (GD_NAME ${GD_LIBRARY_PATH} NAME_WE)
    get_filename_component (GD_PATH ${GD_LIBRARY_PATH} DIRECTORY)
    set (RUNTIME_LIBRARY_PATHS ${RUNTIME_LIBRARY_PATHS} ${GD_PATH}/${GD_NAME}.dll)
endif (WIN32)

set (INCLUDE_PATHS ${INCLUDE_PATHS} ${GD_INCLUDE_PATH})
set (LIBRARY_PATHS ${LIBRARY_PATHS} ${GD_LIBRARY_PATH})


# Find gd include and lib paths

find_path (GD_INCLUDE_PATH
    gd.h
    PATHS /usr/local/include /usr/include
)

if (WIN32)
    SET(GD_NAMES ${GD_NAMES} bgd)
else (WIN32)
    SET(GD_NAMES ${GD_NAMES} gd)
endif (WIN32)

find_library (GD_LIBRARY_PATH
  NAMES ${GD_NAMES}
  PATHS /usr/local/lib /usr/lib
)

if (NOT GD_INCLUDE_PATH OR NOT GD_LIBRARY_PATH)
    message(FATAL_ERROR "Failed to find library gd")
endif (NOT GD_INCLUDE_PATH OR NOT GD_LIBRARY_PATH)

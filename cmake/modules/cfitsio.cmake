# Find cfitsio include and lib paths

find_path (CFITSIO_INCLUDE_PATH
    fitsio.h
    PATHS /usr/local/include /usr/include
)

find_library (CFITSIO_LIBRARY_PATH
  cfitsio
  PATHS /usr/local/lib /usr/lib
)

if (NOT CFITSIO_INCLUDE_PATH OR NOT CFITSIO_LIBRARY_PATH)
    message(FATAL_ERROR "Failed to find library cfitsio")
endif (NOT CFITSIO_INCLUDE_PATH OR NOT CFITSIO_LIBRARY_PATH)

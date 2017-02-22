# Find Boost include and lib paths

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED COMPONENTS program_options)

set (INCLUDE_PATHS ${INCLUDE_PATHS} ${Boost_INCLUDE_DIRS})
set (LIBRARY_PATHS ${LIBRARY_PATHS} ${Boost_LIBRARIES})

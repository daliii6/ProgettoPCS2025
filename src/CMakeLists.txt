list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/Polyhedron.hpp)

# sources = .cpp
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/Utils.cpp)

# eventuali altri header
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/Utils.hpp)

# info sul path locale
list(APPEND src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(src_sources ${src_sources} PARENT_SCOPE)
set(src_headers ${src_headers} PARENT_SCOPE)
set(src_includes ${src_includes} PARENT_SCOPE)

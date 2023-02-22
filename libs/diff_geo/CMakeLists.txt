add_library(diff_geo
  diff_geo.h  diff_geo.cpp
)

target_include_directories(diff_geo PRIVATE .)
target_include_directories(diff_geo PUBLIC ${CMAKE_SOURCE_DIR}/libs)

set_target_properties(diff_geo PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

target_link_libraries(diff_geo yocto)

if(UNIX AND NOT APPLE)
  find_package(Threads REQUIRED)
  target_link_libraries(diff_geo Threads::Threads)
endif(UNIX AND NOT APPLE)

# warning flags
if(APPLE)
  # target_compile_options(src PRIVATE -Wall)
endif(APPLE)
if(MSVC)
# target_compile_options(src PRIVATE /D_CRT_SECURE_NO_WARNINGS /wd4018 /wd4244 /wd4305 /wd4800 /wd4267)
# target_compile_options(src PUBLIC /EHsc)
endif(MSVC)

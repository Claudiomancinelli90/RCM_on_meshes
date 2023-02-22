
add_executable(RCM_GUI main.cpp app.h app.cpp)
add_definitions(-DSHADERS_PATH="${PROJECT_SOURCE_DIR}/shaders/")
add_definitions(-DPROJECT_PATH="${PROJECT_SOURCE_DIR}/")
set_target_properties(RCM_GUI     PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
target_include_directories(RCM_GUI    PUBLIC ${CMAKE_SOURCE_DIR}/libs)
target_link_libraries(RCM_GUI yocto diff_geo yocto_gui realtime utils)

add_subdirectory(yocto)
add_subdirectory(diff_geo)
include_directories(nanosvg/src)
add_subdirectory(realtime)
add_subdirectory(utils)



if(YOCTO_OPENGL)
add_subdirectory(yocto_gui)
endif(YOCTO_OPENGL)

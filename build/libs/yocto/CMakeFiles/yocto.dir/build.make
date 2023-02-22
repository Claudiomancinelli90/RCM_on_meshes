# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.24

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.24.2/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.24.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build

# Include any dependencies generated for this target.
include libs/yocto/CMakeFiles/yocto.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include libs/yocto/CMakeFiles/yocto.dir/compiler_depend.make

# Include the progress variables for this target.
include libs/yocto/CMakeFiles/yocto.dir/progress.make

# Include the compile flags for this target's objects.
include libs/yocto/CMakeFiles/yocto.dir/flags.make

libs/yocto/CMakeFiles/yocto.dir/yocto_modelio.cpp.o: libs/yocto/CMakeFiles/yocto.dir/flags.make
libs/yocto/CMakeFiles/yocto.dir/yocto_modelio.cpp.o: /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_modelio.cpp
libs/yocto/CMakeFiles/yocto.dir/yocto_modelio.cpp.o: libs/yocto/CMakeFiles/yocto.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object libs/yocto/CMakeFiles/yocto.dir/yocto_modelio.cpp.o"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT libs/yocto/CMakeFiles/yocto.dir/yocto_modelio.cpp.o -MF CMakeFiles/yocto.dir/yocto_modelio.cpp.o.d -o CMakeFiles/yocto.dir/yocto_modelio.cpp.o -c /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_modelio.cpp

libs/yocto/CMakeFiles/yocto.dir/yocto_modelio.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/yocto.dir/yocto_modelio.cpp.i"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_modelio.cpp > CMakeFiles/yocto.dir/yocto_modelio.cpp.i

libs/yocto/CMakeFiles/yocto.dir/yocto_modelio.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/yocto.dir/yocto_modelio.cpp.s"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_modelio.cpp -o CMakeFiles/yocto.dir/yocto_modelio.cpp.s

libs/yocto/CMakeFiles/yocto.dir/yocto_bvh.cpp.o: libs/yocto/CMakeFiles/yocto.dir/flags.make
libs/yocto/CMakeFiles/yocto.dir/yocto_bvh.cpp.o: /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_bvh.cpp
libs/yocto/CMakeFiles/yocto.dir/yocto_bvh.cpp.o: libs/yocto/CMakeFiles/yocto.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object libs/yocto/CMakeFiles/yocto.dir/yocto_bvh.cpp.o"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT libs/yocto/CMakeFiles/yocto.dir/yocto_bvh.cpp.o -MF CMakeFiles/yocto.dir/yocto_bvh.cpp.o.d -o CMakeFiles/yocto.dir/yocto_bvh.cpp.o -c /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_bvh.cpp

libs/yocto/CMakeFiles/yocto.dir/yocto_bvh.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/yocto.dir/yocto_bvh.cpp.i"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_bvh.cpp > CMakeFiles/yocto.dir/yocto_bvh.cpp.i

libs/yocto/CMakeFiles/yocto.dir/yocto_bvh.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/yocto.dir/yocto_bvh.cpp.s"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_bvh.cpp -o CMakeFiles/yocto.dir/yocto_bvh.cpp.s

libs/yocto/CMakeFiles/yocto.dir/yocto_shape.cpp.o: libs/yocto/CMakeFiles/yocto.dir/flags.make
libs/yocto/CMakeFiles/yocto.dir/yocto_shape.cpp.o: /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_shape.cpp
libs/yocto/CMakeFiles/yocto.dir/yocto_shape.cpp.o: libs/yocto/CMakeFiles/yocto.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object libs/yocto/CMakeFiles/yocto.dir/yocto_shape.cpp.o"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT libs/yocto/CMakeFiles/yocto.dir/yocto_shape.cpp.o -MF CMakeFiles/yocto.dir/yocto_shape.cpp.o.d -o CMakeFiles/yocto.dir/yocto_shape.cpp.o -c /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_shape.cpp

libs/yocto/CMakeFiles/yocto.dir/yocto_shape.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/yocto.dir/yocto_shape.cpp.i"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_shape.cpp > CMakeFiles/yocto.dir/yocto_shape.cpp.i

libs/yocto/CMakeFiles/yocto.dir/yocto_shape.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/yocto.dir/yocto_shape.cpp.s"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_shape.cpp -o CMakeFiles/yocto.dir/yocto_shape.cpp.s

libs/yocto/CMakeFiles/yocto.dir/yocto_mesh.cpp.o: libs/yocto/CMakeFiles/yocto.dir/flags.make
libs/yocto/CMakeFiles/yocto.dir/yocto_mesh.cpp.o: /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_mesh.cpp
libs/yocto/CMakeFiles/yocto.dir/yocto_mesh.cpp.o: libs/yocto/CMakeFiles/yocto.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object libs/yocto/CMakeFiles/yocto.dir/yocto_mesh.cpp.o"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT libs/yocto/CMakeFiles/yocto.dir/yocto_mesh.cpp.o -MF CMakeFiles/yocto.dir/yocto_mesh.cpp.o.d -o CMakeFiles/yocto.dir/yocto_mesh.cpp.o -c /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_mesh.cpp

libs/yocto/CMakeFiles/yocto.dir/yocto_mesh.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/yocto.dir/yocto_mesh.cpp.i"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_mesh.cpp > CMakeFiles/yocto.dir/yocto_mesh.cpp.i

libs/yocto/CMakeFiles/yocto.dir/yocto_mesh.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/yocto.dir/yocto_mesh.cpp.s"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_mesh.cpp -o CMakeFiles/yocto.dir/yocto_mesh.cpp.s

libs/yocto/CMakeFiles/yocto.dir/yocto_image.cpp.o: libs/yocto/CMakeFiles/yocto.dir/flags.make
libs/yocto/CMakeFiles/yocto.dir/yocto_image.cpp.o: /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_image.cpp
libs/yocto/CMakeFiles/yocto.dir/yocto_image.cpp.o: libs/yocto/CMakeFiles/yocto.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object libs/yocto/CMakeFiles/yocto.dir/yocto_image.cpp.o"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT libs/yocto/CMakeFiles/yocto.dir/yocto_image.cpp.o -MF CMakeFiles/yocto.dir/yocto_image.cpp.o.d -o CMakeFiles/yocto.dir/yocto_image.cpp.o -c /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_image.cpp

libs/yocto/CMakeFiles/yocto.dir/yocto_image.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/yocto.dir/yocto_image.cpp.i"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_image.cpp > CMakeFiles/yocto.dir/yocto_image.cpp.i

libs/yocto/CMakeFiles/yocto.dir/yocto_image.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/yocto.dir/yocto_image.cpp.s"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_image.cpp -o CMakeFiles/yocto.dir/yocto_image.cpp.s

libs/yocto/CMakeFiles/yocto.dir/yocto_trace.cpp.o: libs/yocto/CMakeFiles/yocto.dir/flags.make
libs/yocto/CMakeFiles/yocto.dir/yocto_trace.cpp.o: /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_trace.cpp
libs/yocto/CMakeFiles/yocto.dir/yocto_trace.cpp.o: libs/yocto/CMakeFiles/yocto.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object libs/yocto/CMakeFiles/yocto.dir/yocto_trace.cpp.o"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT libs/yocto/CMakeFiles/yocto.dir/yocto_trace.cpp.o -MF CMakeFiles/yocto.dir/yocto_trace.cpp.o.d -o CMakeFiles/yocto.dir/yocto_trace.cpp.o -c /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_trace.cpp

libs/yocto/CMakeFiles/yocto.dir/yocto_trace.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/yocto.dir/yocto_trace.cpp.i"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_trace.cpp > CMakeFiles/yocto.dir/yocto_trace.cpp.i

libs/yocto/CMakeFiles/yocto.dir/yocto_trace.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/yocto.dir/yocto_trace.cpp.s"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_trace.cpp -o CMakeFiles/yocto.dir/yocto_trace.cpp.s

libs/yocto/CMakeFiles/yocto.dir/yocto_sceneio.cpp.o: libs/yocto/CMakeFiles/yocto.dir/flags.make
libs/yocto/CMakeFiles/yocto.dir/yocto_sceneio.cpp.o: /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_sceneio.cpp
libs/yocto/CMakeFiles/yocto.dir/yocto_sceneio.cpp.o: libs/yocto/CMakeFiles/yocto.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object libs/yocto/CMakeFiles/yocto.dir/yocto_sceneio.cpp.o"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT libs/yocto/CMakeFiles/yocto.dir/yocto_sceneio.cpp.o -MF CMakeFiles/yocto.dir/yocto_sceneio.cpp.o.d -o CMakeFiles/yocto.dir/yocto_sceneio.cpp.o -c /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_sceneio.cpp

libs/yocto/CMakeFiles/yocto.dir/yocto_sceneio.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/yocto.dir/yocto_sceneio.cpp.i"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_sceneio.cpp > CMakeFiles/yocto.dir/yocto_sceneio.cpp.i

libs/yocto/CMakeFiles/yocto.dir/yocto_sceneio.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/yocto.dir/yocto_sceneio.cpp.s"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_sceneio.cpp -o CMakeFiles/yocto.dir/yocto_sceneio.cpp.s

libs/yocto/CMakeFiles/yocto.dir/yocto_commonio.cpp.o: libs/yocto/CMakeFiles/yocto.dir/flags.make
libs/yocto/CMakeFiles/yocto.dir/yocto_commonio.cpp.o: /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_commonio.cpp
libs/yocto/CMakeFiles/yocto.dir/yocto_commonio.cpp.o: libs/yocto/CMakeFiles/yocto.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object libs/yocto/CMakeFiles/yocto.dir/yocto_commonio.cpp.o"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT libs/yocto/CMakeFiles/yocto.dir/yocto_commonio.cpp.o -MF CMakeFiles/yocto.dir/yocto_commonio.cpp.o.d -o CMakeFiles/yocto.dir/yocto_commonio.cpp.o -c /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_commonio.cpp

libs/yocto/CMakeFiles/yocto.dir/yocto_commonio.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/yocto.dir/yocto_commonio.cpp.i"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_commonio.cpp > CMakeFiles/yocto.dir/yocto_commonio.cpp.i

libs/yocto/CMakeFiles/yocto.dir/yocto_commonio.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/yocto.dir/yocto_commonio.cpp.s"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/yocto_commonio.cpp -o CMakeFiles/yocto.dir/yocto_commonio.cpp.s

libs/yocto/CMakeFiles/yocto.dir/ext/stb_image.cpp.o: libs/yocto/CMakeFiles/yocto.dir/flags.make
libs/yocto/CMakeFiles/yocto.dir/ext/stb_image.cpp.o: /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/ext/stb_image.cpp
libs/yocto/CMakeFiles/yocto.dir/ext/stb_image.cpp.o: libs/yocto/CMakeFiles/yocto.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object libs/yocto/CMakeFiles/yocto.dir/ext/stb_image.cpp.o"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT libs/yocto/CMakeFiles/yocto.dir/ext/stb_image.cpp.o -MF CMakeFiles/yocto.dir/ext/stb_image.cpp.o.d -o CMakeFiles/yocto.dir/ext/stb_image.cpp.o -c /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/ext/stb_image.cpp

libs/yocto/CMakeFiles/yocto.dir/ext/stb_image.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/yocto.dir/ext/stb_image.cpp.i"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/ext/stb_image.cpp > CMakeFiles/yocto.dir/ext/stb_image.cpp.i

libs/yocto/CMakeFiles/yocto.dir/ext/stb_image.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/yocto.dir/ext/stb_image.cpp.s"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/ext/stb_image.cpp -o CMakeFiles/yocto.dir/ext/stb_image.cpp.s

libs/yocto/CMakeFiles/yocto.dir/ext/cgltf.cpp.o: libs/yocto/CMakeFiles/yocto.dir/flags.make
libs/yocto/CMakeFiles/yocto.dir/ext/cgltf.cpp.o: /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/ext/cgltf.cpp
libs/yocto/CMakeFiles/yocto.dir/ext/cgltf.cpp.o: libs/yocto/CMakeFiles/yocto.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object libs/yocto/CMakeFiles/yocto.dir/ext/cgltf.cpp.o"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT libs/yocto/CMakeFiles/yocto.dir/ext/cgltf.cpp.o -MF CMakeFiles/yocto.dir/ext/cgltf.cpp.o.d -o CMakeFiles/yocto.dir/ext/cgltf.cpp.o -c /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/ext/cgltf.cpp

libs/yocto/CMakeFiles/yocto.dir/ext/cgltf.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/yocto.dir/ext/cgltf.cpp.i"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/ext/cgltf.cpp > CMakeFiles/yocto.dir/ext/cgltf.cpp.i

libs/yocto/CMakeFiles/yocto.dir/ext/cgltf.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/yocto.dir/ext/cgltf.cpp.s"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/ext/cgltf.cpp -o CMakeFiles/yocto.dir/ext/cgltf.cpp.s

libs/yocto/CMakeFiles/yocto.dir/ext/tinyexr.cpp.o: libs/yocto/CMakeFiles/yocto.dir/flags.make
libs/yocto/CMakeFiles/yocto.dir/ext/tinyexr.cpp.o: /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/ext/tinyexr.cpp
libs/yocto/CMakeFiles/yocto.dir/ext/tinyexr.cpp.o: libs/yocto/CMakeFiles/yocto.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object libs/yocto/CMakeFiles/yocto.dir/ext/tinyexr.cpp.o"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT libs/yocto/CMakeFiles/yocto.dir/ext/tinyexr.cpp.o -MF CMakeFiles/yocto.dir/ext/tinyexr.cpp.o.d -o CMakeFiles/yocto.dir/ext/tinyexr.cpp.o -c /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/ext/tinyexr.cpp

libs/yocto/CMakeFiles/yocto.dir/ext/tinyexr.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/yocto.dir/ext/tinyexr.cpp.i"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/ext/tinyexr.cpp > CMakeFiles/yocto.dir/ext/tinyexr.cpp.i

libs/yocto/CMakeFiles/yocto.dir/ext/tinyexr.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/yocto.dir/ext/tinyexr.cpp.s"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto/ext/tinyexr.cpp -o CMakeFiles/yocto.dir/ext/tinyexr.cpp.s

# Object files for target yocto
yocto_OBJECTS = \
"CMakeFiles/yocto.dir/yocto_modelio.cpp.o" \
"CMakeFiles/yocto.dir/yocto_bvh.cpp.o" \
"CMakeFiles/yocto.dir/yocto_shape.cpp.o" \
"CMakeFiles/yocto.dir/yocto_mesh.cpp.o" \
"CMakeFiles/yocto.dir/yocto_image.cpp.o" \
"CMakeFiles/yocto.dir/yocto_trace.cpp.o" \
"CMakeFiles/yocto.dir/yocto_sceneio.cpp.o" \
"CMakeFiles/yocto.dir/yocto_commonio.cpp.o" \
"CMakeFiles/yocto.dir/ext/stb_image.cpp.o" \
"CMakeFiles/yocto.dir/ext/cgltf.cpp.o" \
"CMakeFiles/yocto.dir/ext/tinyexr.cpp.o"

# External object files for target yocto
yocto_EXTERNAL_OBJECTS =

/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/bin/libyocto.a: libs/yocto/CMakeFiles/yocto.dir/yocto_modelio.cpp.o
/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/bin/libyocto.a: libs/yocto/CMakeFiles/yocto.dir/yocto_bvh.cpp.o
/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/bin/libyocto.a: libs/yocto/CMakeFiles/yocto.dir/yocto_shape.cpp.o
/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/bin/libyocto.a: libs/yocto/CMakeFiles/yocto.dir/yocto_mesh.cpp.o
/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/bin/libyocto.a: libs/yocto/CMakeFiles/yocto.dir/yocto_image.cpp.o
/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/bin/libyocto.a: libs/yocto/CMakeFiles/yocto.dir/yocto_trace.cpp.o
/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/bin/libyocto.a: libs/yocto/CMakeFiles/yocto.dir/yocto_sceneio.cpp.o
/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/bin/libyocto.a: libs/yocto/CMakeFiles/yocto.dir/yocto_commonio.cpp.o
/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/bin/libyocto.a: libs/yocto/CMakeFiles/yocto.dir/ext/stb_image.cpp.o
/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/bin/libyocto.a: libs/yocto/CMakeFiles/yocto.dir/ext/cgltf.cpp.o
/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/bin/libyocto.a: libs/yocto/CMakeFiles/yocto.dir/ext/tinyexr.cpp.o
/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/bin/libyocto.a: libs/yocto/CMakeFiles/yocto.dir/build.make
/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/bin/libyocto.a: libs/yocto/CMakeFiles/yocto.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Linking CXX static library /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/bin/libyocto.a"
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && $(CMAKE_COMMAND) -P CMakeFiles/yocto.dir/cmake_clean_target.cmake
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/yocto.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
libs/yocto/CMakeFiles/yocto.dir/build: /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/bin/libyocto.a
.PHONY : libs/yocto/CMakeFiles/yocto.dir/build

libs/yocto/CMakeFiles/yocto.dir/clean:
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto && $(CMAKE_COMMAND) -P CMakeFiles/yocto.dir/cmake_clean.cmake
.PHONY : libs/yocto/CMakeFiles/yocto.dir/clean

libs/yocto/CMakeFiles/yocto.dir/depend:
	cd /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/libs/yocto /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto /Users/claudiomancinelli/Documents/GitHub/RCM_on_meshes/build/libs/yocto/CMakeFiles/yocto.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : libs/yocto/CMakeFiles/yocto.dir/depend

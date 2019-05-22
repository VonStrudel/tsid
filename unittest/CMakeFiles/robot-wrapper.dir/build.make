# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Produce verbose output by default.
VERBOSE = 1

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /opt/openrobots/src/tsid

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /opt/openrobots/src/tsid

# Include any dependencies generated for this target.
include unittest/CMakeFiles/robot-wrapper.dir/depend.make

# Include the progress variables for this target.
include unittest/CMakeFiles/robot-wrapper.dir/progress.make

# Include the compile flags for this target's objects.
include unittest/CMakeFiles/robot-wrapper.dir/flags.make

unittest/CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.o: unittest/CMakeFiles/robot-wrapper.dir/flags.make
unittest/CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.o: unittest/robot-wrapper.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/opt/openrobots/src/tsid/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object unittest/CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.o"
	cd /opt/openrobots/src/tsid/unittest && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.o -c /opt/openrobots/src/tsid/unittest/robot-wrapper.cpp

unittest/CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.i"
	cd /opt/openrobots/src/tsid/unittest && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /opt/openrobots/src/tsid/unittest/robot-wrapper.cpp > CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.i

unittest/CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.s"
	cd /opt/openrobots/src/tsid/unittest && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /opt/openrobots/src/tsid/unittest/robot-wrapper.cpp -o CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.s

unittest/CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.o.requires:

.PHONY : unittest/CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.o.requires

unittest/CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.o.provides: unittest/CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.o.requires
	$(MAKE) -f unittest/CMakeFiles/robot-wrapper.dir/build.make unittest/CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.o.provides.build
.PHONY : unittest/CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.o.provides

unittest/CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.o.provides.build: unittest/CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.o


# Object files for target robot-wrapper
robot__wrapper_OBJECTS = \
"CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.o"

# External object files for target robot-wrapper
robot__wrapper_EXTERNAL_OBJECTS =

unittest/robot-wrapper: unittest/CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.o
unittest/robot-wrapper: unittest/CMakeFiles/robot-wrapper.dir/build.make
unittest/robot-wrapper: src/libtsid.so
unittest/robot-wrapper: /usr/lib/x86_64-linux-gnu/libboost_unit_test_framework.so
unittest/robot-wrapper: unittest/CMakeFiles/robot-wrapper.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/opt/openrobots/src/tsid/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable robot-wrapper"
	cd /opt/openrobots/src/tsid/unittest && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/robot-wrapper.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
unittest/CMakeFiles/robot-wrapper.dir/build: unittest/robot-wrapper

.PHONY : unittest/CMakeFiles/robot-wrapper.dir/build

unittest/CMakeFiles/robot-wrapper.dir/requires: unittest/CMakeFiles/robot-wrapper.dir/robot-wrapper.cpp.o.requires

.PHONY : unittest/CMakeFiles/robot-wrapper.dir/requires

unittest/CMakeFiles/robot-wrapper.dir/clean:
	cd /opt/openrobots/src/tsid/unittest && $(CMAKE_COMMAND) -P CMakeFiles/robot-wrapper.dir/cmake_clean.cmake
.PHONY : unittest/CMakeFiles/robot-wrapper.dir/clean

unittest/CMakeFiles/robot-wrapper.dir/depend:
	cd /opt/openrobots/src/tsid && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /opt/openrobots/src/tsid /opt/openrobots/src/tsid/unittest /opt/openrobots/src/tsid /opt/openrobots/src/tsid/unittest /opt/openrobots/src/tsid/unittest/CMakeFiles/robot-wrapper.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : unittest/CMakeFiles/robot-wrapper.dir/depend

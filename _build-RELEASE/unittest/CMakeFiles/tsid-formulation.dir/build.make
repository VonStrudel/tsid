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
CMAKE_BINARY_DIR = /opt/openrobots/src/tsid/_build-RELEASE

# Include any dependencies generated for this target.
include unittest/CMakeFiles/tsid-formulation.dir/depend.make

# Include the progress variables for this target.
include unittest/CMakeFiles/tsid-formulation.dir/progress.make

# Include the compile flags for this target's objects.
include unittest/CMakeFiles/tsid-formulation.dir/flags.make

unittest/CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.o: unittest/CMakeFiles/tsid-formulation.dir/flags.make
unittest/CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.o: ../unittest/tsid-formulation.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/opt/openrobots/src/tsid/_build-RELEASE/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object unittest/CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.o"
	cd /opt/openrobots/src/tsid/_build-RELEASE/unittest && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.o -c /opt/openrobots/src/tsid/unittest/tsid-formulation.cpp

unittest/CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.i"
	cd /opt/openrobots/src/tsid/_build-RELEASE/unittest && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /opt/openrobots/src/tsid/unittest/tsid-formulation.cpp > CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.i

unittest/CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.s"
	cd /opt/openrobots/src/tsid/_build-RELEASE/unittest && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /opt/openrobots/src/tsid/unittest/tsid-formulation.cpp -o CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.s

unittest/CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.o.requires:

.PHONY : unittest/CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.o.requires

unittest/CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.o.provides: unittest/CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.o.requires
	$(MAKE) -f unittest/CMakeFiles/tsid-formulation.dir/build.make unittest/CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.o.provides.build
.PHONY : unittest/CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.o.provides

unittest/CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.o.provides.build: unittest/CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.o


# Object files for target tsid-formulation
tsid__formulation_OBJECTS = \
"CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.o"

# External object files for target tsid-formulation
tsid__formulation_EXTERNAL_OBJECTS =

unittest/tsid-formulation: unittest/CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.o
unittest/tsid-formulation: unittest/CMakeFiles/tsid-formulation.dir/build.make
unittest/tsid-formulation: src/libtsid.so
unittest/tsid-formulation: /usr/lib/x86_64-linux-gnu/libboost_unit_test_framework.so
unittest/tsid-formulation: unittest/CMakeFiles/tsid-formulation.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/opt/openrobots/src/tsid/_build-RELEASE/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable tsid-formulation"
	cd /opt/openrobots/src/tsid/_build-RELEASE/unittest && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tsid-formulation.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
unittest/CMakeFiles/tsid-formulation.dir/build: unittest/tsid-formulation

.PHONY : unittest/CMakeFiles/tsid-formulation.dir/build

unittest/CMakeFiles/tsid-formulation.dir/requires: unittest/CMakeFiles/tsid-formulation.dir/tsid-formulation.cpp.o.requires

.PHONY : unittest/CMakeFiles/tsid-formulation.dir/requires

unittest/CMakeFiles/tsid-formulation.dir/clean:
	cd /opt/openrobots/src/tsid/_build-RELEASE/unittest && $(CMAKE_COMMAND) -P CMakeFiles/tsid-formulation.dir/cmake_clean.cmake
.PHONY : unittest/CMakeFiles/tsid-formulation.dir/clean

unittest/CMakeFiles/tsid-formulation.dir/depend:
	cd /opt/openrobots/src/tsid/_build-RELEASE && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /opt/openrobots/src/tsid /opt/openrobots/src/tsid/unittest /opt/openrobots/src/tsid/_build-RELEASE /opt/openrobots/src/tsid/_build-RELEASE/unittest /opt/openrobots/src/tsid/_build-RELEASE/unittest/CMakeFiles/tsid-formulation.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : unittest/CMakeFiles/tsid-formulation.dir/depend


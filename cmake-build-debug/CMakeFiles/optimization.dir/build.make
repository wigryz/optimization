# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/macbook/Pliki/AGH/optimization

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/macbook/Pliki/AGH/optimization/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/optimization.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/optimization.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/optimization.dir/flags.make

CMakeFiles/optimization.dir/src/main.cpp.o: CMakeFiles/optimization.dir/flags.make
CMakeFiles/optimization.dir/src/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/macbook/Pliki/AGH/optimization/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/optimization.dir/src/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/optimization.dir/src/main.cpp.o -c /Users/macbook/Pliki/AGH/optimization/src/main.cpp

CMakeFiles/optimization.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/optimization.dir/src/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/macbook/Pliki/AGH/optimization/src/main.cpp > CMakeFiles/optimization.dir/src/main.cpp.i

CMakeFiles/optimization.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/optimization.dir/src/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/macbook/Pliki/AGH/optimization/src/main.cpp -o CMakeFiles/optimization.dir/src/main.cpp.s

CMakeFiles/optimization.dir/src/matrix.cpp.o: CMakeFiles/optimization.dir/flags.make
CMakeFiles/optimization.dir/src/matrix.cpp.o: ../src/matrix.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/macbook/Pliki/AGH/optimization/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/optimization.dir/src/matrix.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/optimization.dir/src/matrix.cpp.o -c /Users/macbook/Pliki/AGH/optimization/src/matrix.cpp

CMakeFiles/optimization.dir/src/matrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/optimization.dir/src/matrix.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/macbook/Pliki/AGH/optimization/src/matrix.cpp > CMakeFiles/optimization.dir/src/matrix.cpp.i

CMakeFiles/optimization.dir/src/matrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/optimization.dir/src/matrix.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/macbook/Pliki/AGH/optimization/src/matrix.cpp -o CMakeFiles/optimization.dir/src/matrix.cpp.s

CMakeFiles/optimization.dir/src/ode_solver.cpp.o: CMakeFiles/optimization.dir/flags.make
CMakeFiles/optimization.dir/src/ode_solver.cpp.o: ../src/ode_solver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/macbook/Pliki/AGH/optimization/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/optimization.dir/src/ode_solver.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/optimization.dir/src/ode_solver.cpp.o -c /Users/macbook/Pliki/AGH/optimization/src/ode_solver.cpp

CMakeFiles/optimization.dir/src/ode_solver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/optimization.dir/src/ode_solver.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/macbook/Pliki/AGH/optimization/src/ode_solver.cpp > CMakeFiles/optimization.dir/src/ode_solver.cpp.i

CMakeFiles/optimization.dir/src/ode_solver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/optimization.dir/src/ode_solver.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/macbook/Pliki/AGH/optimization/src/ode_solver.cpp -o CMakeFiles/optimization.dir/src/ode_solver.cpp.s

CMakeFiles/optimization.dir/src/opt_alg.cpp.o: CMakeFiles/optimization.dir/flags.make
CMakeFiles/optimization.dir/src/opt_alg.cpp.o: ../src/opt_alg.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/macbook/Pliki/AGH/optimization/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/optimization.dir/src/opt_alg.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/optimization.dir/src/opt_alg.cpp.o -c /Users/macbook/Pliki/AGH/optimization/src/opt_alg.cpp

CMakeFiles/optimization.dir/src/opt_alg.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/optimization.dir/src/opt_alg.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/macbook/Pliki/AGH/optimization/src/opt_alg.cpp > CMakeFiles/optimization.dir/src/opt_alg.cpp.i

CMakeFiles/optimization.dir/src/opt_alg.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/optimization.dir/src/opt_alg.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/macbook/Pliki/AGH/optimization/src/opt_alg.cpp -o CMakeFiles/optimization.dir/src/opt_alg.cpp.s

CMakeFiles/optimization.dir/src/solution.cpp.o: CMakeFiles/optimization.dir/flags.make
CMakeFiles/optimization.dir/src/solution.cpp.o: ../src/solution.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/macbook/Pliki/AGH/optimization/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/optimization.dir/src/solution.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/optimization.dir/src/solution.cpp.o -c /Users/macbook/Pliki/AGH/optimization/src/solution.cpp

CMakeFiles/optimization.dir/src/solution.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/optimization.dir/src/solution.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/macbook/Pliki/AGH/optimization/src/solution.cpp > CMakeFiles/optimization.dir/src/solution.cpp.i

CMakeFiles/optimization.dir/src/solution.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/optimization.dir/src/solution.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/macbook/Pliki/AGH/optimization/src/solution.cpp -o CMakeFiles/optimization.dir/src/solution.cpp.s

# Object files for target optimization
optimization_OBJECTS = \
"CMakeFiles/optimization.dir/src/main.cpp.o" \
"CMakeFiles/optimization.dir/src/matrix.cpp.o" \
"CMakeFiles/optimization.dir/src/ode_solver.cpp.o" \
"CMakeFiles/optimization.dir/src/opt_alg.cpp.o" \
"CMakeFiles/optimization.dir/src/solution.cpp.o"

# External object files for target optimization
optimization_EXTERNAL_OBJECTS =

optimization: CMakeFiles/optimization.dir/src/main.cpp.o
optimization: CMakeFiles/optimization.dir/src/matrix.cpp.o
optimization: CMakeFiles/optimization.dir/src/ode_solver.cpp.o
optimization: CMakeFiles/optimization.dir/src/opt_alg.cpp.o
optimization: CMakeFiles/optimization.dir/src/solution.cpp.o
optimization: CMakeFiles/optimization.dir/build.make
optimization: CMakeFiles/optimization.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/macbook/Pliki/AGH/optimization/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable optimization"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/optimization.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/optimization.dir/build: optimization
.PHONY : CMakeFiles/optimization.dir/build

CMakeFiles/optimization.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/optimization.dir/cmake_clean.cmake
.PHONY : CMakeFiles/optimization.dir/clean

CMakeFiles/optimization.dir/depend:
	cd /Users/macbook/Pliki/AGH/optimization/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/macbook/Pliki/AGH/optimization /Users/macbook/Pliki/AGH/optimization /Users/macbook/Pliki/AGH/optimization/cmake-build-debug /Users/macbook/Pliki/AGH/optimization/cmake-build-debug /Users/macbook/Pliki/AGH/optimization/cmake-build-debug/CMakeFiles/optimization.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/optimization.dir/depend

# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.12.3/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.12.3/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/ziqwang/Documents/GitHub/nanogui

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ziqwang/Documents/GitHub/nanogui/build

# Include any dependencies generated for this target.
include CMakeFiles/nanogui-python.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/nanogui-python.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/nanogui-python.dir/flags.make

# Object files for target nanogui-python
nanogui__python_OBJECTS =

# External object files for target nanogui-python
nanogui__python_EXTERNAL_OBJECTS = \
"/Users/ziqwang/Documents/GitHub/nanogui/build/CMakeFiles/nanogui-python-obj.dir/python/main.cpp.o" \
"/Users/ziqwang/Documents/GitHub/nanogui/build/CMakeFiles/nanogui-python-obj.dir/python/constants_glfw.cpp.o" \
"/Users/ziqwang/Documents/GitHub/nanogui/build/CMakeFiles/nanogui-python-obj.dir/python/constants_entypo.cpp.o" \
"/Users/ziqwang/Documents/GitHub/nanogui/build/CMakeFiles/nanogui-python-obj.dir/python/eigen.cpp.o" \
"/Users/ziqwang/Documents/GitHub/nanogui/build/CMakeFiles/nanogui-python-obj.dir/python/widget.cpp.o" \
"/Users/ziqwang/Documents/GitHub/nanogui/build/CMakeFiles/nanogui-python-obj.dir/python/layout.cpp.o" \
"/Users/ziqwang/Documents/GitHub/nanogui/build/CMakeFiles/nanogui-python-obj.dir/python/basics.cpp.o" \
"/Users/ziqwang/Documents/GitHub/nanogui/build/CMakeFiles/nanogui-python-obj.dir/python/button.cpp.o" \
"/Users/ziqwang/Documents/GitHub/nanogui/build/CMakeFiles/nanogui-python-obj.dir/python/tabs.cpp.o" \
"/Users/ziqwang/Documents/GitHub/nanogui/build/CMakeFiles/nanogui-python-obj.dir/python/textbox.cpp.o" \
"/Users/ziqwang/Documents/GitHub/nanogui/build/CMakeFiles/nanogui-python-obj.dir/python/glcanvas.cpp.o" \
"/Users/ziqwang/Documents/GitHub/nanogui/build/CMakeFiles/nanogui-python-obj.dir/python/formhelper.cpp.o" \
"/Users/ziqwang/Documents/GitHub/nanogui/build/CMakeFiles/nanogui-python-obj.dir/python/misc.cpp.o" \
"/Users/ziqwang/Documents/GitHub/nanogui/build/CMakeFiles/nanogui-python-obj.dir/python/glutil.cpp.o" \
"/Users/ziqwang/Documents/GitHub/nanogui/build/CMakeFiles/nanogui-python-obj.dir/python/nanovg.cpp.o" \
"/Users/ziqwang/Documents/GitHub/nanogui/build/CMakeFiles/nanogui-python-obj.dir/ext/coro/coro.c.o"

python/nanogui.cpython-36m-darwin.so: CMakeFiles/nanogui-python-obj.dir/python/main.cpp.o
python/nanogui.cpython-36m-darwin.so: CMakeFiles/nanogui-python-obj.dir/python/constants_glfw.cpp.o
python/nanogui.cpython-36m-darwin.so: CMakeFiles/nanogui-python-obj.dir/python/constants_entypo.cpp.o
python/nanogui.cpython-36m-darwin.so: CMakeFiles/nanogui-python-obj.dir/python/eigen.cpp.o
python/nanogui.cpython-36m-darwin.so: CMakeFiles/nanogui-python-obj.dir/python/widget.cpp.o
python/nanogui.cpython-36m-darwin.so: CMakeFiles/nanogui-python-obj.dir/python/layout.cpp.o
python/nanogui.cpython-36m-darwin.so: CMakeFiles/nanogui-python-obj.dir/python/basics.cpp.o
python/nanogui.cpython-36m-darwin.so: CMakeFiles/nanogui-python-obj.dir/python/button.cpp.o
python/nanogui.cpython-36m-darwin.so: CMakeFiles/nanogui-python-obj.dir/python/tabs.cpp.o
python/nanogui.cpython-36m-darwin.so: CMakeFiles/nanogui-python-obj.dir/python/textbox.cpp.o
python/nanogui.cpython-36m-darwin.so: CMakeFiles/nanogui-python-obj.dir/python/glcanvas.cpp.o
python/nanogui.cpython-36m-darwin.so: CMakeFiles/nanogui-python-obj.dir/python/formhelper.cpp.o
python/nanogui.cpython-36m-darwin.so: CMakeFiles/nanogui-python-obj.dir/python/misc.cpp.o
python/nanogui.cpython-36m-darwin.so: CMakeFiles/nanogui-python-obj.dir/python/glutil.cpp.o
python/nanogui.cpython-36m-darwin.so: CMakeFiles/nanogui-python-obj.dir/python/nanovg.cpp.o
python/nanogui.cpython-36m-darwin.so: CMakeFiles/nanogui-python-obj.dir/ext/coro/coro.c.o
python/nanogui.cpython-36m-darwin.so: CMakeFiles/nanogui-python.dir/build.make
python/nanogui.cpython-36m-darwin.so: libnanogui.dylib
python/nanogui.cpython-36m-darwin.so: CMakeFiles/nanogui-python.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/ziqwang/Documents/GitHub/nanogui/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Linking CXX shared library python/nanogui.cpython-36m-darwin.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/nanogui-python.dir/link.txt --verbose=$(VERBOSE)
	strip -u -r /Users/ziqwang/Documents/GitHub/nanogui/build/python/nanogui.cpython-36m-darwin.so

# Rule to build all files generated by this target.
CMakeFiles/nanogui-python.dir/build: python/nanogui.cpython-36m-darwin.so

.PHONY : CMakeFiles/nanogui-python.dir/build

CMakeFiles/nanogui-python.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/nanogui-python.dir/cmake_clean.cmake
.PHONY : CMakeFiles/nanogui-python.dir/clean

CMakeFiles/nanogui-python.dir/depend:
	cd /Users/ziqwang/Documents/GitHub/nanogui/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ziqwang/Documents/GitHub/nanogui /Users/ziqwang/Documents/GitHub/nanogui /Users/ziqwang/Documents/GitHub/nanogui/build /Users/ziqwang/Documents/GitHub/nanogui/build /Users/ziqwang/Documents/GitHub/nanogui/build/CMakeFiles/nanogui-python.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/nanogui-python.dir/depend


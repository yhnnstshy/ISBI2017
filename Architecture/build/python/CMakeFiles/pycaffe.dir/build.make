# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build

# Include any dependencies generated for this target.
include python/CMakeFiles/pycaffe.dir/depend.make

# Include the progress variables for this target.
include python/CMakeFiles/pycaffe.dir/progress.make

# Include the compile flags for this target's objects.
include python/CMakeFiles/pycaffe.dir/flags.make

python/CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.o: python/CMakeFiles/pycaffe.dir/flags.make
python/CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.o: ../python/caffe/_caffe.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object python/CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.o"
	cd /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/python && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.o -c /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/python/caffe/_caffe.cpp

python/CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.i"
	cd /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/python && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/python/caffe/_caffe.cpp > CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.i

python/CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.s"
	cd /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/python && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/python/caffe/_caffe.cpp -o CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.s

python/CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.o.requires:
.PHONY : python/CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.o.requires

python/CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.o.provides: python/CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.o.requires
	$(MAKE) -f python/CMakeFiles/pycaffe.dir/build.make python/CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.o.provides.build
.PHONY : python/CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.o.provides

python/CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.o.provides.build: python/CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.o

# Object files for target pycaffe
pycaffe_OBJECTS = \
"CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.o"

# External object files for target pycaffe
pycaffe_EXTERNAL_OBJECTS =

lib/_caffe.so: python/CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.o
lib/_caffe.so: python/CMakeFiles/pycaffe.dir/build.make
lib/_caffe.so: lib/libcaffe.so
lib/_caffe.so: /usr/lib/x86_64-linux-gnu/libpython2.7.so
lib/_caffe.so: /usr/lib/x86_64-linux-gnu/libboost_python.so
lib/_caffe.so: lib/libproto.a
lib/_caffe.so: /usr/lib/x86_64-linux-gnu/libboost_system.so
lib/_caffe.so: /usr/lib/x86_64-linux-gnu/libboost_thread.so
lib/_caffe.so: /usr/lib/x86_64-linux-gnu/libpthread.so
lib/_caffe.so: /usr/local/lib/libglog.so
lib/_caffe.so: /usr/local/lib/libgflags.a
lib/_caffe.so: /usr/lib/x86_64-linux-gnu/libprotobuf.so
lib/_caffe.so: /usr/local/lib/libglog.so
lib/_caffe.so: /usr/local/lib/libgflags.a
lib/_caffe.so: /usr/lib/x86_64-linux-gnu/libprotobuf.so
lib/_caffe.so: /usr/lib/x86_64-linux-gnu/libhdf5_hl.so
lib/_caffe.so: /usr/lib/x86_64-linux-gnu/libhdf5.so
lib/_caffe.so: /usr/local/lib/liblmdb.so
lib/_caffe.so: /usr/lib/x86_64-linux-gnu/libleveldb.so
lib/_caffe.so: /usr/lib/libsnappy.so
lib/_caffe.so: /usr/local/cuda-7.0/lib64/libcudart.so
lib/_caffe.so: /usr/local/cuda-7.0/lib64/libcurand.so
lib/_caffe.so: /usr/local/cuda-7.0/lib64/libcublas.so
lib/_caffe.so: /usr/lib/x86_64-linux-gnu/libopencv_highgui.so.2.4.8
lib/_caffe.so: /usr/lib/x86_64-linux-gnu/libopencv_imgproc.so.2.4.8
lib/_caffe.so: /usr/lib/x86_64-linux-gnu/libopencv_core.so.2.4.8
lib/_caffe.so: /usr/lib/liblapack_atlas.so
lib/_caffe.so: /usr/lib/libcblas.so
lib/_caffe.so: /usr/lib/libatlas.so
lib/_caffe.so: /usr/lib/x86_64-linux-gnu/libpython2.7.so
lib/_caffe.so: /usr/lib/x86_64-linux-gnu/libboost_python.so
lib/_caffe.so: python/CMakeFiles/pycaffe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared library ../lib/_caffe.so"
	cd /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/python && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pycaffe.dir/link.txt --verbose=$(VERBOSE)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Creating symlink /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/python/caffe/_caffe.so -> /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/lib/_caffe.so"
	cd /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/python && ln -sf /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/lib/_caffe.so /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/python/caffe/_caffe.so
	cd /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/python && /usr/bin/cmake -E make_directory /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/python/caffe/proto
	cd /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/python && touch /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/python/caffe/proto/__init__.py
	cd /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/python && cp /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/include/caffe/proto/*.py /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/python/caffe/proto/

# Rule to build all files generated by this target.
python/CMakeFiles/pycaffe.dir/build: lib/_caffe.so
.PHONY : python/CMakeFiles/pycaffe.dir/build

python/CMakeFiles/pycaffe.dir/requires: python/CMakeFiles/pycaffe.dir/caffe/_caffe.cpp.o.requires
.PHONY : python/CMakeFiles/pycaffe.dir/requires

python/CMakeFiles/pycaffe.dir/clean:
	cd /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/python && $(CMAKE_COMMAND) -P CMakeFiles/pycaffe.dir/cmake_clean.cmake
.PHONY : python/CMakeFiles/pycaffe.dir/clean

python/CMakeFiles/pycaffe.dir/depend:
	cd /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/python /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/python /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/python/CMakeFiles/pycaffe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : python/CMakeFiles/pycaffe.dir/depend


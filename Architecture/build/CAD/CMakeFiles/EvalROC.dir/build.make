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
include CAD/CMakeFiles/EvalROC.dir/depend.make

# Include the progress variables for this target.
include CAD/CMakeFiles/EvalROC.dir/progress.make

# Include the compile flags for this target's objects.
include CAD/CMakeFiles/EvalROC.dir/flags.make

CAD/CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.o: CAD/CMakeFiles/EvalROC.dir/flags.make
CAD/CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.o: ../CAD/EvalProstateCAD.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CAD/CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.o"
	cd /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/CAD && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.o -c /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/CAD/EvalProstateCAD.cpp

CAD/CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.i"
	cd /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/CAD && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/CAD/EvalProstateCAD.cpp > CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.i

CAD/CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.s"
	cd /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/CAD && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/CAD/EvalProstateCAD.cpp -o CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.s

CAD/CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.o.requires:
.PHONY : CAD/CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.o.requires

CAD/CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.o.provides: CAD/CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.o.requires
	$(MAKE) -f CAD/CMakeFiles/EvalROC.dir/build.make CAD/CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.o.provides.build
.PHONY : CAD/CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.o.provides

CAD/CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.o.provides.build: CAD/CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.o

CAD/CMakeFiles/EvalROC.dir/bsdgetopt.c.o: CAD/CMakeFiles/EvalROC.dir/flags.make
CAD/CMakeFiles/EvalROC.dir/bsdgetopt.c.o: ../CAD/bsdgetopt.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CAD/CMakeFiles/EvalROC.dir/bsdgetopt.c.o"
	cd /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/CAD && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/EvalROC.dir/bsdgetopt.c.o   -c /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/CAD/bsdgetopt.c

CAD/CMakeFiles/EvalROC.dir/bsdgetopt.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/EvalROC.dir/bsdgetopt.c.i"
	cd /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/CAD && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/CAD/bsdgetopt.c > CMakeFiles/EvalROC.dir/bsdgetopt.c.i

CAD/CMakeFiles/EvalROC.dir/bsdgetopt.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/EvalROC.dir/bsdgetopt.c.s"
	cd /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/CAD && /usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/CAD/bsdgetopt.c -o CMakeFiles/EvalROC.dir/bsdgetopt.c.s

CAD/CMakeFiles/EvalROC.dir/bsdgetopt.c.o.requires:
.PHONY : CAD/CMakeFiles/EvalROC.dir/bsdgetopt.c.o.requires

CAD/CMakeFiles/EvalROC.dir/bsdgetopt.c.o.provides: CAD/CMakeFiles/EvalROC.dir/bsdgetopt.c.o.requires
	$(MAKE) -f CAD/CMakeFiles/EvalROC.dir/build.make CAD/CMakeFiles/EvalROC.dir/bsdgetopt.c.o.provides.build
.PHONY : CAD/CMakeFiles/EvalROC.dir/bsdgetopt.c.o.provides

CAD/CMakeFiles/EvalROC.dir/bsdgetopt.c.o.provides.build: CAD/CMakeFiles/EvalROC.dir/bsdgetopt.c.o

# Object files for target EvalROC
EvalROC_OBJECTS = \
"CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.o" \
"CMakeFiles/EvalROC.dir/bsdgetopt.c.o"

# External object files for target EvalROC
EvalROC_EXTERNAL_OBJECTS =

CAD/EvalROC: CAD/CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.o
CAD/EvalROC: CAD/CMakeFiles/EvalROC.dir/bsdgetopt.c.o
CAD/EvalROC: CAD/CMakeFiles/EvalROC.dir/build.make
CAD/EvalROC: CAD/libCADLib.so
CAD/EvalROC: /usr/local/lib/libITKCommon-4.11.a
CAD/EvalROC: /usr/local/lib/libitksys-4.11.a
CAD/EvalROC: /usr/local/lib/libitkvnl_algo-4.11.a
CAD/EvalROC: /usr/local/lib/libitkvnl-4.11.a
CAD/EvalROC: /usr/local/lib/libitkv3p_netlib-4.11.a
CAD/EvalROC: /usr/local/lib/libitknetlib-4.11.a
CAD/EvalROC: /usr/local/lib/libitkvcl-4.11.a
CAD/EvalROC: /usr/local/lib/libITKIOImageBase-4.11.a
CAD/EvalROC: /usr/local/lib/libitkgdcmDICT-4.11.a
CAD/EvalROC: /usr/local/lib/libitkgdcmMSFF-4.11.a
CAD/EvalROC: /usr/local/lib/libITKEXPAT-4.11.a
CAD/EvalROC: /usr/local/lib/libitkzlib-4.11.a
CAD/EvalROC: /usr/local/lib/libITKIOGDCM-4.11.a
CAD/EvalROC: /usr/local/lib/libITKIOMeta-4.11.a
CAD/EvalROC: /usr/local/lib/libITKMetaIO-4.11.a
CAD/EvalROC: /usr/local/lib/libITKIONIFTI-4.11.a
CAD/EvalROC: /usr/local/lib/libITKIONRRD-4.11.a
CAD/EvalROC: /usr/local/lib/libITKIOJPEG-4.11.a
CAD/EvalROC: /usr/local/lib/libITKIOPNG-4.11.a
CAD/EvalROC: /usr/local/lib/libITKMesh-4.11.a
CAD/EvalROC: /usr/local/lib/libITKQuadEdgeMesh-4.11.a
CAD/EvalROC: /usr/local/lib/libitkgdcmMSFF-4.11.a
CAD/EvalROC: /usr/local/lib/libitkgdcmDICT-4.11.a
CAD/EvalROC: /usr/local/lib/libitkgdcmIOD-4.11.a
CAD/EvalROC: /usr/local/lib/libITKEXPAT-4.11.a
CAD/EvalROC: /usr/local/lib/libitkgdcmDSED-4.11.a
CAD/EvalROC: /usr/local/lib/libitkgdcmCommon-4.11.a
CAD/EvalROC: /usr/local/lib/libitkgdcmjpeg8-4.11.a
CAD/EvalROC: /usr/local/lib/libitkgdcmjpeg12-4.11.a
CAD/EvalROC: /usr/local/lib/libitkgdcmjpeg16-4.11.a
CAD/EvalROC: /usr/local/lib/libitkgdcmopenjpeg-4.11.a
CAD/EvalROC: /usr/local/lib/libitkgdcmcharls-4.11.a
CAD/EvalROC: /usr/local/lib/libitkgdcmuuid-4.11.a
CAD/EvalROC: /usr/local/lib/libITKniftiio-4.11.a
CAD/EvalROC: /usr/local/lib/libITKznz-4.11.a
CAD/EvalROC: /usr/local/lib/libITKNrrdIO-4.11.a
CAD/EvalROC: /usr/local/lib/libitkjpeg-4.11.a
CAD/EvalROC: /usr/local/lib/libITKIOImageBase-4.11.a
CAD/EvalROC: /usr/local/lib/libitkpng-4.11.a
CAD/EvalROC: /usr/local/lib/libitkzlib-4.11.a
CAD/EvalROC: /usr/local/lib/libITKMesh-4.11.a
CAD/EvalROC: /usr/local/lib/libITKTransform-4.11.a
CAD/EvalROC: /usr/local/lib/libITKCommon-4.11.a
CAD/EvalROC: /usr/local/lib/libitksys-4.11.a
CAD/EvalROC: /usr/local/lib/libitkdouble-conversion-4.11.a
CAD/EvalROC: /usr/local/lib/libITKVNLInstantiation-4.11.a
CAD/EvalROC: /usr/local/lib/libitkvnl_algo-4.11.a
CAD/EvalROC: /usr/local/lib/libitkvnl-4.11.a
CAD/EvalROC: /usr/local/lib/libitkv3p_netlib-4.11.a
CAD/EvalROC: /usr/local/lib/libitknetlib-4.11.a
CAD/EvalROC: /usr/local/lib/libitkvcl-4.11.a
CAD/EvalROC: CAD/CMakeFiles/EvalROC.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable EvalROC"
	cd /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/CAD && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/EvalROC.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CAD/CMakeFiles/EvalROC.dir/build: CAD/EvalROC
.PHONY : CAD/CMakeFiles/EvalROC.dir/build

CAD/CMakeFiles/EvalROC.dir/requires: CAD/CMakeFiles/EvalROC.dir/EvalProstateCAD.cpp.o.requires
CAD/CMakeFiles/EvalROC.dir/requires: CAD/CMakeFiles/EvalROC.dir/bsdgetopt.c.o.requires
.PHONY : CAD/CMakeFiles/EvalROC.dir/requires

CAD/CMakeFiles/EvalROC.dir/clean:
	cd /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/CAD && $(CMAKE_COMMAND) -P CMakeFiles/EvalROC.dir/cmake_clean.cmake
.PHONY : CAD/CMakeFiles/EvalROC.dir/clean

CAD/CMakeFiles/EvalROC.dir/depend:
	cd /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/CAD /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/CAD /home/tsehayyk/Work/ProstateCAD/ISBI2017/Architecture/build/CAD/CMakeFiles/EvalROC.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CAD/CMakeFiles/EvalROC.dir/depend


cmake_minimum_required(VERSION 3.6)
project(SynImageGen)

set(CMAKE_CXX_STANDARD 14)

# Find the ArrayFire package.
FIND_PACKAGE(ArrayFire REQUIRED)
# If ArrayFire is found, the following variables will be defined:
#
# ArrayFire_INCLUDE_DIRS    - Location of ArrayFire's include directory.
# ArrayFire_LIBRARIES       - Location of ArrayFire's libraries. This will default
#                             to a GPU backend if one is found.
# ArrayFire_FOUND           - True if ArrayFire has been located
#
# You may provide a hint to where ArrayFire's root directory may be located
# by setting ArrayFire_DIR.
#
# ----------------------------------------------------------------------------
#
# ArrayFire_CPU_FOUND        - True of the ArrayFire CPU library has been found.
# ArrayFire_CPU_LIBRARIES    - Location of ArrayFire's CPU library, if found
# ArrayFire_CUDA_FOUND       - True of the ArrayFire CUDA library has been found.
# ArrayFire_CUDA_LIBRARIES   - Location of ArrayFire's CUDA library, if found
# ArrayFire_OpenCL_FOUND     - True of the ArrayFire OpenCL library has been found.
# ArrayFire_OpenCL_LIBRARIES - Location of ArrayFire's OpenCL library, if found

set(SOURCE_FILES main.cpp src/MeshDataAF.cpp src/MeshDataAF.h)

# Include the ArrayFire hreaders
INCLUDE_DIRECTORIES(${ArrayFire_INCLUDE_DIRS})

# ArrayFire OpenCL backend
FIND_PACKAGE(OpenCL)
IF(${ArrayFire_OpenCL_FOUND} AND ${OpenCL_FOUND})
    # We need to find OpenCL as transitive linking is disabled on some OSes
    MESSAGE(STATUS "ArrayFire OpenCL backend found.")
    ADD_EXECUTABLE(SynImageGen ${SOURCE_FILES} )
    TARGET_LINK_LIBRARIES(SynImageGen ${ArrayFire_OpenCL_LIBRARIES}
            ${OpenCL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
ENDIF()

if(${ArrayFire_CPU_FOUND})
    MESSAGE(STATUS "ArrayFire CPU backend found.")
    ADD_EXECUTABLE(SynImageGen_cpu ${SOURCE_FILES})
    TARGET_LINK_LIBRARIES(SynImageGen_cpu ${ArrayFire_CPU_LIBRARIES}
            ${CMAKE_THREAD_LIBS_INIT})
ENDIF()
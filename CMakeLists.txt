cmake_minimum_required(VERSION 3.2)

project(homeworks)



message(2)
message(STATUS "CMAKE_SYSTEM_NAME: ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_SYSTEM_VERSION: ${CMAKE_SYSTEM_VERSION}")

if(APPLE)
    message(STATUS "Detected Apple System")
else()
    message(STATUS "Not an Apple System")
endif()
if(APPLE)
    FIND_LIBRARY(COCOA_LIBRARY Cocoa)
    FIND_LIBRARY(OpenGL_LIBRARY OpenGL)
    FIND_LIBRARY(IOKit_LIBRARY IOKit)
    FIND_LIBRARY(CoreVideo_LIBRARY CoreVideo)
    MARK_AS_ADVANCED(COCOA_LIBRARY OpenGL_LIBRARY)
    SET(APPLE_LIBS ${COCOA_LIBRARY} ${IOKit_LIBRARY} ${OpenGL_LIBRARY} ${CoreVideo_LIBRARY})
endif()

# include(${CMAKE_TOOLCHAIN_FILE})

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
include(libigl)

igl_include(glfw)
igl_include(imgui)
igl_include(opengl)
# /Users/liuxing/Documents/WASM-TOOL/emsdk/upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake



set(DATA_FOLDER ${CMAKE_CURRENT_SOURCE_DIR}/data)
set(DEST_DATA_FOLDER ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${DATA_FOLDER} DESTINATION ${DEST_DATA_FOLDER})

add_subdirectory(a1-mass-spring-1d)
add_subdirectory(a2-mass-spring-3d)
add_subdirectory(a3-finite-elements-3d)
add_subdirectory(a4-cloth-simulation)
add_subdirectory(a5-rigid-bodies)
add_subdirectory(a6-rigid-bodies-contact)

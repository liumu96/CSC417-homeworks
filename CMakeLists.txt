cmake_minimum_required(VERSION 3.2)

project(homeworks)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
include(libigl)

igl_include(glfw)
igl_include(imgui)
igl_include(opengl)

if(APPLE)
    FIND_LIBRARY(COCOA_LIBRARY Cocoa)
    FIND_LIBRARY(OpenGL_LIBRARY OpenGL)
    FIND_LIBRARY(IOKit_LIBRARY IOKit)
    FIND_LIBRARY(CoreVideo_LIBRARY CoreVideo)
    MARK_AS_ADVANCED(COCOA_LIBRARY OpenGL_LIBRARY)
    SET(APPLE_LIBS ${COCOA_LIBRARY} ${IOKit_LIBRARY} ${OpenGL_LIBRARY} ${CoreVideo_LIBRARY})
endif()

set(DATA_FOLDER ${CMAKE_CURRENT_SOURCE_DIR}/data)
set(DEST_DATA_FOLDER ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${DATA_FOLDER} DESTINATION ${DEST_DATA_FOLDER})

add_subdirectory(a1-mass-spring-1d)


cmake_minimum_required(VERSION 3.2)

project(mass-spring-1d)

#include directories
# include_directories(${PROJECT_SOURCE_DIR}/shared/include)
# include_directories(${PROJECT_SOURCE_DIR}/include)

#include source code
# file(GLOB H_FILES ${PROJECT_SOURCE_DIR}/include/*.h)
# file(GLOB H_SHARED_FILES ${PROJECT_SOURCE_DIR}/shared/include/*.h)
file(GLOB CPP_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)
# file(GLOB CPP_SHARED_FILES ${PROJECT_SOURCE_DIR}/shared/src/*.cpp)

add_executable(${PROJECT_NAME} 
main.cpp 
# ${H_FILES} 
# ${H_SHARED_FILES} 
${CPP_FILES} 

# ${CPP_SHARED_FILES}
)
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl igl::glfw igl::imgui )
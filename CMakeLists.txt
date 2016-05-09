cmake_minimum_required(VERSION 3.0)
project(Glitter)

find_package(CUDA REQUIRED)

#CUDA NVCC SETTINGS
list(APPEND CUDA_NVCC_FLAGS "-arch=sm_20;-std=c++11;-O2;-DVERBOSE")
SET(CUDA_PROPAGATE_HOST_FLAGS OFF)

option(GLFW_BUILD_DOCS OFF)
option(GLFW_BUILD_EXAMPLES OFF)
option(GLFW_BUILD_TESTS OFF)
add_subdirectory(Glitter/Vendor/glfw)

if(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -std=c++11")
    if(NOT WIN32)
        set(GLAD_LIBRARIES dl)
    endif()
endif()

link_directories(libs/)

include_directories(Glitter/Headers/
                    Glitter/Vendor/glad/include/
                    Glitter/Vendor/glfw/include/
                    Glitter/Vendor/glm/)

file(GLOB VENDORS_SOURCES Glitter/Vendor/glad/src/glad.c)
file(GLOB PROJECT_HEADERS Glitter/Headers/*.hpp
			  Glitter/Headers/*.h)
file(GLOB PROJECT_SOURCES Glitter/Sources/*.cpp
                          Glitter/Sources/*.cu)
file(GLOB PROJECT_SHADERS Glitter/Shaders/*.cs
                          Glitter/Shaders/*.fs
                          Glitter/Shaders/*.gs
                          Glitter/Shaders/*.vs)
file(GLOB PROJECT_CONFIGS CMakeLists.txt
                          Readme.md
                         .gitattributes
                         .gitignore
                         .gitmodules)

source_group("Headers" FILES ${PROJECT_HEADERS})
source_group("Shaders" FILES ${PROJECT_SHADERS})
source_group("Sources" FILES ${PROJECT_SOURCES})
source_group("Vendors" FILES ${VENDORS_SOURCES})

add_definitions(-DGLFW_INCLUDE_NONE
                -DPROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\")
cuda_add_executable(${PROJECT_NAME} ${PROJECT_SOURCES} ${PROJECT_HEADERS}
                               ${PROJECT_SHADERS} ${PROJECT_CONFIGS}
                               ${VENDORS_SOURCES})
target_link_libraries(${PROJECT_NAME} glfw
                      ${GLFW_LIBRARIES} ${GLAD_LIBRARIES} soil2 GL)
set_target_properties(${PROJECT_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${PROJECT_NAME})

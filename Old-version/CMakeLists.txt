cmake_minimum_required(VERSION 3.6)
project(myFemProject)

set(EXEC_NAME myFem)

# OpenGL
find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})

# GLFW
add_subdirectory(glfw)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/glfw/include)

# GMSH
if(MSVC OR MSYS OR MINGW)
    set(GMSH_LIBRARY ${CMAKE_SOURCE_DIR}/gmsh/gmsh-4.11.1-Windows64-sdk/lib/gmsh.lib)
    set(GMSH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/gmsh/gmsh-4.11.1-Windows64-sdk/include)
    #Copy .dll into Debug
    set(GMSH_DLL ${CMAKE_SOURCE_DIR}/gmsh/gmsh-4.11.1-Windows64-sdk/lib/gmsh-4.11.dll)
    file(COPY ${GMSH_DLL} DESTINATION ${CMAKE_BINARY_DIR}/${CMAKE_BUILD_TYPE})
endif(MSVC OR MSYS OR MINGW)

if(APPLE)
    if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm64")
        set(GMSH_LIBRARY ${CMAKE_SOURCE_DIR}/gmsh/gmsh-4.11.1-MacOSARM-sdk/lib/libgmsh.dylib)
        set(GMSH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/gmsh/gmsh-4.11.1-MacOSARM-sdk/include)
    else()    
        set(GMSH_LIBRARY ${CMAKE_SOURCE_DIR}/gmsh/gmsh-4.11.1-MacOSX-sdk/lib/libgmsh.dylib)
        set(GMSH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/gmsh/gmsh-4.11.1-MacOSX-sdk/include)
    endif()
endif(APPLE)

if(UNIX AND NOT APPLE)
    set(GMSH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/gmsh/gmsh-4.11.1-Linux64-sdk/include)
    set(GMSH_LIBRARY ${CMAKE_SOURCE_DIR}/gmsh/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so)
endif(UNIX AND NOT APPLE)

message("GMSH_INCLUDE_DIR : ${GMSH_INCLUDE_DIR}")
message("GMSH_LIBRARY : ${GMSH_LIBRARY}")
include_directories(${GMSH_INCLUDE_DIR})


# link the aforementioned libraries
add_executable(${EXEC_NAME} src/fem.c src/glfem.c src/homework.c src/main.c)
target_link_libraries(${EXEC_NAME} glfw ${OPENGL_gl_LIBRARY} ${GMSH_LIBRARY})

if(APPLE)
    target_compile_options(${EXEC_NAME} PRIVATE "-Wdeprecated-declarations")
    target_compile_definitions(${EXEC_NAME} PRIVATE -D GL_SILENCE_DEPRECATION)
endif()

# set exec as the startup project in visual studio
if (MSVC)
    set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${EXEC_NAME})
    target_compile_options(${EXEC_NAME} PRIVATE "/W2")
endif()

target_compile_definitions(${EXEC_NAME} PRIVATE graphic)

# add a custom target so we can do make run directly
add_custom_target(run COMMAND ${EXEC_NAME} DEPENDS ${EXEC_NAME} COMMENT "./${EXEC_NAME}")


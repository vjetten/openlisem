cmake_minimum_required(VERSION 3.16)

project(Lisem)

# Instruct CMake to run moc, uic, and rcc automatically when needed.
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find the Qt 6 libraries
find_package(Qt6 REQUIRED COMPONENTS Widgets Core Gui Network)

# Include directories for Qt 6
include_directories(${Qt6Widgets_INCLUDE_DIRS})
include_directories(${Qt6Core_INCLUDE_DIRS})
include_directories(${Qt6Gui_INCLUDE_DIRS})
include_directories(${Qt6Network_INCLUDE_DIRS})

# Add definitions for Qt 6
add_definitions(${Qt6Widgets_DEFINITIONS})
add_definitions(${Qt6Core_DEFINITIONS})
add_definitions(${Qt6Gui_DEFINITIONS})
add_definitions(${Qt6Network_DEFINITIONS})

# Add compile flags for Qt 6
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt6Widgets_EXECUTABLE_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt6Core_EXECUTABLE_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt6Gui_EXECUTABLE_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt6Network_EXECUTABLE_COMPILE_FLAGS}")

set(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_RELEASE} -fprofile-arcs -ftest-coverage")

# Update target properties for Qt 6
set_target_properties(Qt6::Core PROPERTIES MAP_IMPORTED_CONFIG_COVERAGE "RELEASE")

SET(CMAKE_MODULE_PATH
    $ENV{CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
)

# Enable ccache if available
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
endif()

# Include application makefile
INCLUDE(Application)

target_sources(Lisem
  PRIVATE
  AUTHORS
  README.txt
  changes.txt  # changes is a short list with latest changes for the installer
  compile_lisem_windows.txt
  compile_lisem_linux.txt
)


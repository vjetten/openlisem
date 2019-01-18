cmake_minimum_required(VERSION 2.8.11)

project(Lisem)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

set_property(SOURCE ui_lisemqt.h PROPERTY SKIP_AUTOMOC ON)
set_property(SOURCE ui_lisemqt.h PROPERTY SKIP_AUTOUIC ON)
set_property(SOURCE ui_lisemqt.h PROPERTY SKIP_AUTOGEN ON)
set_property(SOURCE qrc_openlisem.cpp PROPERTY SKIP_AUTOMOC ON)
set_property(SOURCE qrc_openlisem.cpp PROPERTY SKIP_AUTOUIC ON)
set_property(SOURCE qrc_openlisem.cpp PROPERTY SKIP_AUTOGEN ON)

# Find the QtWidgets library
find_package(Qt5Widgets)
find_package(Qt5Core)
find_package(Qt5Gui)
find_package(Qt5OpenGL)
find_package(OpenGL)

include_directories(${Qt5Widgets_INCLUDE_DIRS})
include_directories(${Qt5Core_INCLUDE_DIRS})
include_directories(${Qt5Gui_INCLUDE_DIRS})
include_directories(${Qt5OpenGL_INCLUDE_DIRS})
include_directories( ${OPENGL_INCLUDE_DIRS} )

add_definitions(${Qt5Widgets_DEFINITIONS})
add_definitions(${Qt5Core_DEFINITIONS})
add_definitions(${Qt5Gui_DEFINITIONS})
add_definitions(${Qt5Widgets_DEFINITIONS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Core_EXECUTABLE_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Gui_EXECUTABLE_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5OpenGL_EXECUTABLE_COMPILE_FLAGS}")

set(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_RELEASE} -fprofile-arcs -ftest-coverage")

set_target_properties(Qt5::Core PROPERTIES MAP_IMPORTED_CONFIG_COVERAGE "RELEASE")

SET(CMAKE_MODULE_PATH
    $ENV{CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
)


# Include application makefile
INCLUDE(Application)



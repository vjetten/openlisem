# USER: specify locations of all libraries on the local machine

IF(UNIX AND NOT CYGWIN)
    SET(CMAKE_SKIP_BUILD_RPATH FALSE)
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

    SET(QWT_LIBRARIES "${QWT_BUILD_DIR}/lib/libqwt.so")
    SET(QWT_INCLUDE_DIRS "${QWT_BUILD_DIR}/include/")
ENDIF()

# assuming MSYS instalation: https://wiki.qt.io/MSYS2
IF(WIN32)
    SET(GDAL_INCLUDE_DIRS "${GDAL_BUILD_DIR}/include")
    SET(GDAL_LIBRARIES "${GDAL_BUILD_DIR}/lib/libgdal.dll.a")

    #SET(QWT_INCLUDE_DIRS "${QWT_BUILD_DIR}/include/qwt") # qwt build outside msys
    #SET(QWT_LIBRARIES "${QWT_BUILD_DIR}/lib/libqwt.dll.a")

    SET(QWT_INCLUDE_DIRS "${QWT_BUILD_DIR}/include")
    SET(QWT_LIBRARIES "${QWT_BUILD_DIR}/lib/libqwt.dll.a")

    SET(OMP_INCLUDE_DIRS "C:/msys/mingw64/lib/gcc/x86_64-w64-mingw32/11.2.0/include")
ENDIF()


INCLUDE(CheckCXXCompilerFlag)

IF(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wcast-qual -Wwrite-strings -Wno-sign-conversion -Werror=strict-aliasing -std=c++11 -fopenmp")
    IF(UNIX)
       SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -Wl,-rpath=${ORIGIN}./lib")
       # extra flags for thread and looking for so libs in ./lib
    ENDIF()
ENDIF()


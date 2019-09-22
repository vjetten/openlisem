
# Find packages. ---------------------------------------------------------------
# this calls fndq
FIND_PACKAGE(GDAL REQUIRED)
FIND_PACKAGE(Qwt REQUIRED)
FIND_PACKAGE(PCRasterRasterFormat REQUIRED)

# Fixup GDAL variable.
IF(WIN32 AND MINGW)
    # On Windows, find_package always selects the release dll named
    # libgdal.dll. This is a symbolic links which trips the linker. Also,
    # for debug builds, we want to link against the debug library.
    # Before fixing this, verify that the variable contains a single string
    # and that the name is libgdal.dll. The CMake module may be fixed at
    # some point.
    LIST(LENGTH GDAL_LIBRARIES NR_GDAL_LIBRARIES)
    IF(NR_GDAL_LIBRARIES EQUAL 1)
        GET_FILENAME_COMPONENT(GDAL_LIBRARY_NAME ${GDAL_LIBRARIES} NAME)

        IF(GDAL_LIBRARY_NAME STREQUAL "libgdal.dll")
            SET(GDAL_VERSION "1.11.1")
            SET(OPTIMIZED_GDAL_LIBRARY ${GDAL_LIBRARIES}.${GDAL_VERSION})
            SET(DEBUG_GDAL_LIBRARY ${GDAL_LIBRARIES}.${GDAL_VERSION})
            STRING(REPLACE gdal gdald DEBUG_GDAL_LIBRARY ${DEBUG_GDAL_LIBRARY})
            SET(GDAL_LIBRARIES
                optimized ${OPTIMIZED_GDAL_LIBRARY}
                debug ${DEBUG_GDAL_LIBRARY})
        ENDIF()
    ENDIF()
ENDIF()



# Configure project. -----------------------------------------------------------

INCLUDE_DIRECTORIES(
    SYSTEM
    ${GDAL_INCLUDE_DIRS}
    ${QWT_INCLUDE_DIRS}
    ${PCRASTER_RASTER_FORMAT_INCLUDE_DIRS}
)


SET(LISEM_EXTERNAL_LIBRARIES
    ${GDAL_LIBRARIES}
    ${QWT_LIBRARIES}
    ${PCRASTER_RASTER_FORMAT_LIBRARIES}
)

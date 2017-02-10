IF((NOT DEFINED LISEM2016_3RD_PARTY_ROOT) AND (NOT DEFINED
        ENV{LISEM2016_3RD_PARTY_ROOT}))
    # Third party software can be built using Peacock. In that case the
    # PEACOCK_PREFIX environment or CMake variable should be set. Given its
    # value the LISEM2016_3RD_PARTY_ROOT variable is set.
    # PEACOCK_PREFIX is only used if LISEM2016_3RD_PARTY_ROOT is not set.
    IF((NOT DEFINED PEACOCK_PREFIX) AND (DEFINED ENV{PEACOCK_PREFIX}))
        SET(PEACOCK_PREFIX $ENV{PEACOCK_PREFIX})
    ENDIF()


    IF(DEFINED PEACOCK_PREFIX)
        IF(EXISTS ${PEACOCK_PREFIX}/lisem)
            SET(PEACOCK_PREFIX ${PEACOCK_PREFIX}/lisem)
        ENDIF()

        SET(LISEM2016_3RD_PARTY_ROOT ${PEACOCK_PREFIX}/${peacock_target_platform})
    ENDIF()
ENDIF()


IF(NOT DEFINED LISEM2016_3RD_PARTY_ROOT AND DEFINED ENV{LISEM2016_3RD_PARTY_ROOT})
    SET(LISEM2016_3RD_PARTY_ROOT $ENV{LISEM2016_3RD_PARTY_ROOT})
ENDIF()


IF(DEFINED LISEM2016_3RD_PARTY_ROOT)
    # Assume that all 3rd party software that LISEM depends on is rooted at
    # $LISEM2016_3RD_PARTY_ROOT. If not, the user can edit the cache variables.

    # Boost.
    FILE(GLOB DEFAULT_PATH ${LISEM2016_3RD_PARTY_ROOT}/boost-*)
    IF(NOT DEFAULT_PATH)
        SET(DEFAULT_PATH ${LISEM2016_3RD_PARTY_ROOT})
    ENDIF()
    SET(LISEM_BOOST_ROOT
        ${DEFAULT_PATH}
        CACHE PATH "Path to root of Boost software"
    )

    # Qt.
#    FILE(GLOB DEFAULT_PATH ${LISEM2016_3RD_PARTY_ROOT}/qt-*)
#    IF(NOT DEFAULT_PATH)
#        SET(DEFAULT_PATH ${LISEM2016_3RD_PARTY_ROOT})
#    ENDIF()
#    SET(LISEM_QT_ROOT
#        ${DEFAULT_PATH}
#        CACHE PATH "Path to root of Qt software"
#    )

    # Qwt.
    FILE(GLOB DEFAULT_PATH ${LISEM2016_3RD_PARTY_ROOT}/qwt-*)
    IF(NOT DEFAULT_PATH)
        SET(DEFAULT_PATH ${LISEM2016_3RD_PARTY_ROOT})
    ENDIF()
    SET(LISEM_QWT_ROOT
        ${DEFAULT_PATH}
        CACHE PATH "Path to root of Qwt software"
    )

    # PCRaster raster format.
    FILE(GLOB DEFAULT_PATH ${LISEM2016_3RD_PARTY_ROOT}/pcraster_raster_format-*)
    IF(NOT DEFAULT_PATH)
        SET(DEFAULT_PATH ${LISEM2016_3RD_PARTY_ROOT})
    ENDIF()
    SET(LISEM_RASTER_FORMAT_ROOT
        ${DEFAULT_PATH}
        CACHE PATH "Path to root of PCRaster raster format software"
    )

    # GDAL.
    FILE(GLOB DEFAULT_PATH ${LISEM2016_3RD_PARTY_ROOT}/gdal-*)
    IF(NOT DEFAULT_PATH)
        SET(DEFAULT_PATH ${LISEM2016_3RD_PARTY_ROOT})
    ENDIF()
    SET(LISEM_GDAL_ROOT
        ${DEFAULT_PATH}
        CACHE PATH "Path to root of GDAL software"
    )

    # Fern.
    FILE(GLOB DEFAULT_PATH ${LISEM2016_3RD_PARTY_ROOT}/fern-*)
    IF(NOT DEFAULT_PATH)
        SET(DEFAULT_PATH ${LISEM2016_3RD_PARTY_ROOT})
    ENDIF()
    SET(LISEM_FERN_ROOT
        ${DEFAULT_PATH}
        CACHE PATH "Path to root of Fern software"
    )
ELSE()
    # Add cache variables so the user can tell us where the various 3rd party
    # software packages are located.
    SET(LISEM_BOOST_ROOT
        ""
        CACHE PATH "Path to root of Boost software"
    )
#    SET(LISEM_QT_ROOT
#        ""
#        CACHE PATH "Path to root of Qt software"
#    )
    SET(LISEM_QWT_ROOT
        ""
        CACHE PATH "Path to root of Qwt software"
    )
    SET(LISEM_RASTER_FORMAT_ROOT
        ""
        CACHE PATH "Path to root of PCRaster raster format software"
    )
    SET(LISEM_GDAL_ROOT
        ""
        CACHE PATH "Path to root of GDAL software"
    )
    SET(LISEM_FERN_ROOT
        ""
        CACHE PATH "Path to root of Fern software"
    )
ENDIF()


# In case we are provided with paths to 3rd party software, then we need to
# tell cmake to look into these locations.
IF(LISEM_BOOST_ROOT)
    LIST(APPEND CMAKE_PREFIX_PATH
        ${LISEM_BOOST_ROOT}
    )
ENDIF()
#IF(LISEM_QT_ROOT)
#    LIST(APPEND CMAKE_PREFIX_PATH
#        ${LISEM_QT_ROOT}
#    )
#ENDIF()
IF(LISEM_QWT_ROOT)
    LIST(APPEND CMAKE_PREFIX_PATH
        ${LISEM_QWT_ROOT}
    )
ENDIF()
IF(LISEM_RASTER_FORMAT_ROOT)
    LIST(APPEND CMAKE_PREFIX_PATH
        ${LISEM_RASTER_FORMAT_ROOT}
    )
ENDIF()
IF(LISEM_GDAL_ROOT)
    LIST(APPEND CMAKE_PREFIX_PATH
        ${LISEM_GDAL_ROOT}
    )
ENDIF()
IF(LISEM_FERN_ROOT)
    LIST(APPEND CMAKE_PREFIX_PATH
        ${LISEM_FERN_ROOT}
    )
ENDIF()

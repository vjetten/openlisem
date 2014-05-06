IF(DEFINED ENV{LISEM_3RD_PARTY_ROOT})
    # Assume that all 3rd party software that LISEM depends on is rooted at
    # $LISEM_3RD_PARTY_ROOT. If not, the user can edit the cache variables.
    SET(LISEM_3RD_PARTY_ROOT
        $ENV{LISEM_3RD_PARTY_ROOT}
    )

    # Qt.
    FILE(GLOB DEFAULT_PATH ${LISEM_3RD_PARTY_ROOT}/qt-*)
    IF(NOT DEFAULT_PATH)
        SET(DEFAULT_PATH ${LISEM_3RD_PARTY_ROOT})
    ENDIF()
    SET(LISEM_QT_ROOT
        ${DEFAULT_PATH}
        CACHE PATH "Path to root of Qt software"
    )

    # Qwt.
    FILE(GLOB DEFAULT_PATH ${LISEM_3RD_PARTY_ROOT}/qwt-*)
    IF(NOT DEFAULT_PATH)
        SET(DEFAULT_PATH ${LISEM_3RD_PARTY_ROOT})
    ENDIF()
    SET(LISEM_QWT_ROOT
        ${DEFAULT_PATH}
        CACHE PATH "Path to root of Qwt software"
    )

    # PCRaster raster format.
    FILE(GLOB DEFAULT_PATH ${LISEM_3RD_PARTY_ROOT}/pcraster_raster_format-*)
    IF(NOT DEFAULT_PATH)
        SET(DEFAULT_PATH ${LISEM_3RD_PARTY_ROOT})
    ENDIF()
    SET(LISEM_RASTER_FORMAT_ROOT
        ${DEFAULT_PATH}
        CACHE PATH "Path to root of PCRaster raster format software"
    )
ELSE()
    # Add cache variables so the user can tell us where the various 3rd party
    # software packages are located.
    SET(LISEM_QT_ROOT
        ""
        CACHE PATH "Path to root of Qt software"
    )
    SET(LISEM_QWT_ROOT
        ""
        CACHE PATH "Path to root of Qwt software"
    )
    SET(LISEM_RASTER_FORMAT_ROOT
        ""
        CACHE PATH "Path to root of PCRaster raster format software"
    )
ENDIF()


# In case we are provided with paths to 3rd party software, then we need to
# tell cmake to look into these locations.
IF(LISEM_QT_ROOT)
    LIST(APPEND CMAKE_PREFIX_PATH
        ${LISEM_QT_ROOT}
    )
ENDIF()
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

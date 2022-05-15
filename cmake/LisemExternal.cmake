
# Find packages. ---------------------------------------------------------------


#===== OMP =====
if(OpenMP_CXX_FOUND)
    target_link_libraries(OpenMP::OpenMP_CXX)
endif()

MARK_AS_ADVANCED(
    OMP_INCLUDE_DIRS
)

FIND_PATH(OMP_INCLUDE_DIRS
    NAMES omp.h
)

#===== QWT =====

FIND_LIBRARY(QWT_LIBRARIES
    NAMES qwt
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Qwt
    REQUIRED_VARS
        QWT_LIBRARIES
        QWT_INCLUDE_DIRS
)


MARK_AS_ADVANCED(
    QWT_LIBRARIES
    QWT_INCLUDE_DIRS
)

FIND_PATH(QWT_INCLUDE_DIRS
    NAMES qwt.h
)

FIND_LIBRARY(QWT_LIBRARIES
    NAMES qwt
)


# Configure project. -----------------------------------------------------------

INCLUDE_DIRECTORIES(
    SYSTEM
    ${GDAL_INCLUDE_DIRS}
    ${QWT_INCLUDE_DIRS}
    ${OMP_INCLUDE_DIRS}
)


SET(LISEM_EXTERNAL_LIBRARIES
    ${GDAL_LIBRARIES}
    ${QWT_LIBRARIES}
)

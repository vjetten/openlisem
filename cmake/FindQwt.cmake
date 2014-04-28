# Qt Widgets for Technical Applications is available at
# http://qwt.sourceforge.net/
#
# This module defines the following variables:
#  QWT_FOUND
#  QWT_INCLUDE_DIRS
#  QWT_LIBRARIES
FIND_PATH(QWT_INCLUDE_DIRS
    NAMES qwt_plot.h
)


FIND_LIBRARY(QWT_LIBRARIES
    NAMES qwt
)
IF(WIN32)
    FIND_LIBRARY(QWT_DEBUG_LIBRARY
        NAMES qwtd
    )
    SET(QWT_LIBRARIES
        optimized ${QWT_LIBRARIES}
        debug ${QWT_DEBUG_LIBRARY}
    )
ENDIF()


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

INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/ui_full
    ${CMAKE_CURRENT_BINARY_DIR}/.
)

SET(LIB_SOURCES
    CsfMap
    error
    fixture
    io
)
ADD_LIBRARY(liblisem STATIC
    ${LIB_SOURCES}
)
SET_TARGET_PROPERTIES(liblisem
    PROPERTIES
        OUTPUT_NAME lisem)

SET(APP_SOURCES
    lisTotalsMB
    ui_full/LisUItreecheck
    ui_full/LisUIModel
    ui_full/LisUIrunfile
    ui_full/LisUImapnames
    ui_full/LisUItreeitem
    ui_full/LisUItreemodel
    ui_full/LisUIDefaultNames
    ui_full/lisemqt
    ui_full/LisUIplot
    ui_full/LisUImapplot
    lisDataInit
    lisErosion
    lisInfiltration
    lisModel
    lisRainintc
    lisUnifiedFlow
    lisUnifiedFlowCommon
    lisUnifiedFlowMomentum2
    lisUnifiedFlowMomentumBalance
    lisUnifiedFlowAdvection2
    lisUnifiedFlowSediment
    lisUnifiedFlowEntrainmentDeposition
    lisUnifiedFlowTimestep
    lisUnifiedFlowMath
    lisUnifiedFlowMuscle
    lisUnifiedFlowConnection
    lisUnifiedFlowBoundary
    lisUnifiedFlowBarriers
    lisUnifiedFlowInitialize
    lisUnifiedFlowInitialize
    include/lisUnifiedFlowThreadPool.h
    include/lisUnifiedFlowThread.h
    lisUnifiedFlowThreadPool.cpp
    lisUnifiedFlowThread.cpp
    lisSlopeStability
    lisReportfile
    lisRunfile
    lisSurfstor
    lisSnowmelt
    main
    operation
    swatre/swatstep
    swatre/swatinit
    swatre/soillut
    swatre/lutio
    swatre/lookup
    swatre/swatinp
    lisInterception
    include/version.h
    include/model.h
    include/UFmodel.h
    include/TMmapVariables.h
    include/LisUIoutput.h
    openlisemico.rc
)
QT4_WRAP_CPP(MOC_SOURCES
    include/lisUnifiedFlowThreadPool.h
    include/lisUnifiedFlowThread.h
    include/model.h
    ui_full/lisemqt.h
    ui_full/LisUItreemodel.h
)
QT4_WRAP_UI(UI_SOURCES
    ui_full/lisemqt.ui
)
QT4_ADD_RESOURCES(RCC_SOURCES
    resources/openlisem.qrc
)
ADD_EXECUTABLE(lisem WIN32
    ${MOC_SOURCES}
    ${UI_SOURCES}
    ${RCC_SOURCES}
    ${APP_SOURCES}
)
TARGET_LINK_LIBRARIES(lisem
    liblisem
    ${LISEM_EXTERNAL_LIBRARIES}
    stdc++
)

# TODO CONFIG += exceptions
# TODO CONFIG += precompile_header
# TODO PRECOMPILED_HEADER = include/stable.h
# TODO RC_FILE = openlisemico.rc

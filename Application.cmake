cmake_minimum_required(VERSION 2.8.11)

INCLUDE(LISEM)

# GET ALL THE APPLICATION SOURCES
INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/ui_full
    ${CMAKE_CURRENT_SOURCE_DIR}/ui_full/3D
    ${CMAKE_CURRENT_SOURCE_DIR}/ui_full/3D/World
    ${CMAKE_CURRENT_SOURCE_DIR}/ui_full/3D/Objects
    ${CMAKE_CURRENT_SOURCE_DIR}/ui_full/3D/Graphics
    ${CMAKE_CURRENT_BINARY_DIR}/.
) 

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
    ui_full/3D/Common
    ui_full/3D/GL3DWidget
    ui_full/3D/GL3DDrawFunctions
    ui_full/3D/World/GL3DWorld
    ui_full/3D/World/GL3DCamera
    ui_full/3D/Objects/GL3DBuilding
    ui_full/3D/Objects/GL3DFlowSurface
    ui_full/3D/Objects/GL3DGrass
    ui_full/3D/Objects/GL3DOcean
    ui_full/3D/Objects/GL3DSkyBox
    ui_full/3D/Objects/GL3DSurface
    ui_full/3D/Objects/GL3DTree
    ui_full/3D/Objects/GL3DSkyBox
    ui_full/3D/Objects/GL3DBuilding
    ui_full/3D/Objects/GL3DObject
    ui_full/3D/Graphics/GL3DMath
    ui_full/3D/Graphics/GL3DModels
    ui_full/3D/Graphics/GL3DShaders
    ui_full/3D/Graphics/GL3DTextures
    ui_full/3D/Graphics/GL3DGeometry
    ui_full/3D/Graphics/GL3DColorRamp
    ui_full/3D/Graphics/GL3DMapMath
    ui_full/3D/World/GL3DCameraController
    ui_full/3D/GL3DInput
    ui_full/3D/GL3DWorldCreator
    lisUnifiedFlowThreadPool.cpp
    lisUnifiedFlowThread.cpp
    lisSlopeStability
    lisReportfile
    lisRunfile
    lisSurfstor
    lisSnowmelt
    main
    CsfMap.cpp
    io.cpp
    fixture.cpp
    error.cpp
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
    include/TMmapVariables.h
    include/LisUIoutput.h
    include/Csfmap.h
    include/array.h
    include/error.h
    include/fixture.h
    include/global.h
    include/io.h
    include/LisUIoutput.h
    include/masked_raster.h
    include/mmath.h
    include/model.h
    include/operation.h
    include/option.h
    include/raster.h
    include/stable.h
    include/swatre_g.h
    include/swatre_p.h
    include/swatreLookup.h
    include/swatremisc.h
    include/swatresoillut.h
    openlisemico.rc
)

QT5_WRAP_UI(UI_SOURCES
    ui_full/lisemqt.ui
)
QT5_ADD_RESOURCES(RCC_SOURCES
    resources/openlisem.qrc
)

# Tell CMake to create the executable
add_executable(Lisem WIN32
    ${UI_SOURCES}
    ${RCC_SOURCES}
    ${APP_SOURCES}
    )

# Use the Widgets module from Qt 5.
target_link_libraries(Lisem ${OPENGL_gl_LIBRARY} Qt5::Widgets Qt5::Gui Qt5::Core Qt5::OpenGL ${LISEM_EXTERNAL_LIBRARIES})

TARGET_LINK_LIBRARIES(Lisem ${QWT_LIBRARIES})

# CMakeLists.txt for BMI wrapper

cmake_minimum_required(VERSION 3.16)

# Find Python and pybind11
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

# BMI C++ library
add_library(bmi_openlisem_core SHARED
    bmi_openlisem.cpp
)

target_include_directories(bmi_openlisem_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/include
    ${Qt6Core_INCLUDE_DIRS}
    ${GDAL_INCLUDE_DIRS}
)

target_link_libraries(bmi_openlisem_core
    Qt6::Core
    Qt6::Gui
    ${GDAL_LIBRARIES}
    ${OPENSSL_LIBRARIES}
    ${CURL_LIBRARIES}
    OpenMP::OpenMP_CXX
)

# Link against the main lisem library sources
# We need to compile the core model files without GUI
set(LISEM_CORE_SOURCES
    ${CMAKE_SOURCE_DIR}/maps/CsfMap.cpp
    ${CMAKE_SOURCE_DIR}/maps/CsfRGBMap.cpp
    ${CMAKE_SOURCE_DIR}/maps/error.cpp
    ${CMAKE_SOURCE_DIR}/maps/fixture.cpp
    ${CMAKE_SOURCE_DIR}/maps/io.cpp
    ${CMAKE_SOURCE_DIR}/maps/operation.cpp
    ${CMAKE_SOURCE_DIR}/model/lisDataInit.cpp
    ${CMAKE_SOURCE_DIR}/model/lisRunfile.cpp
    ${CMAKE_SOURCE_DIR}/model/lisModel.cpp
    ${CMAKE_SOURCE_DIR}/model/lisTotalsMB.cpp
    ${CMAKE_SOURCE_DIR}/model/lisReportfile.cpp
    ${CMAKE_SOURCE_DIR}/model/lisReportmaps.cpp
    ${CMAKE_SOURCE_DIR}/hydro/lisInfiltration.cpp
    ${CMAKE_SOURCE_DIR}/hydro/lisInterception.cpp
    ${CMAKE_SOURCE_DIR}/hydro/lisPercolation.cpp
    ${CMAKE_SOURCE_DIR}/hydro/lisSoilmoisture.cpp
    ${CMAKE_SOURCE_DIR}/hydro/lisSurfstor.cpp
    ${CMAKE_SOURCE_DIR}/flow/lisOverlandflow.cpp
    ${CMAKE_SOURCE_DIR}/flow/lisKinematic.cpp
    ${CMAKE_SOURCE_DIR}/flow/lisSWOF2DopenMUSCL.cpp
    ${CMAKE_SOURCE_DIR}/flow/lisGWflow.cpp
    ${CMAKE_SOURCE_DIR}/flow/lisTiledrainflow.cpp
    ${CMAKE_SOURCE_DIR}/flow/lisFlowBarriers.cpp
    ${CMAKE_SOURCE_DIR}/channel/lisChannelflow.cpp
    ${CMAKE_SOURCE_DIR}/channel/lisChannelflood.cpp
    ${CMAKE_SOURCE_DIR}/channel/lisChannelAux.cpp
    ${CMAKE_SOURCE_DIR}/erosion/lisErosion.cpp
    ${CMAKE_SOURCE_DIR}/erosion/lisErosionSplash.cpp
    ${CMAKE_SOURCE_DIR}/erosion/lisChannelErosion.cpp
    ${CMAKE_SOURCE_DIR}/erosion/lisSWOF2DSediment.cpp
    ${CMAKE_SOURCE_DIR}/swatre/lisInfilSwatre.cpp
    ${CMAKE_SOURCE_DIR}/swatre/swatstep.cpp
    ${CMAKE_SOURCE_DIR}/swatre/swatinit.cpp
    ${CMAKE_SOURCE_DIR}/swatre/lookup.cpp
    ${CMAKE_SOURCE_DIR}/swatre/swatinp.cpp
    ${CMAKE_SOURCE_DIR}/meteo/lisRainfall.cpp
    ${CMAKE_SOURCE_DIR}/meteo/lisEvaporation.cpp
    ${CMAKE_SOURCE_DIR}/meteo/lisSnowmelt.cpp
    ${CMAKE_SOURCE_DIR}/PCR/attravail.c
    ${CMAKE_SOURCE_DIR}/PCR/attrsize.c
    ${CMAKE_SOURCE_DIR}/PCR/csfglob.c
    ${CMAKE_SOURCE_DIR}/PCR/csfsup.c
    ${CMAKE_SOURCE_DIR}/PCR/gattravl.c
    ${CMAKE_SOURCE_DIR}/PCR/gattrblk.c
    ${CMAKE_SOURCE_DIR}/PCR/gcellrep.c
    ${CMAKE_SOURCE_DIR}/PCR/gdattype.c
    ${CMAKE_SOURCE_DIR}/PCR/ggisfid.c
    ${CMAKE_SOURCE_DIR}/PCR/getattr.c
    ${CMAKE_SOURCE_DIR}/PCR/getproj.c
    ${CMAKE_SOURCE_DIR}/PCR/gputproj.c
    ${CMAKE_SOURCE_DIR}/PCR/gvalscal.c
    ${CMAKE_SOURCE_DIR}/PCR/gvartype.c
    ${CMAKE_SOURCE_DIR}/PCR/gversion.c
    ${CMAKE_SOURCE_DIR}/PCR/mclose.c
    ${CMAKE_SOURCE_DIR}/PCR/mopen.c
    ${CMAKE_SOURCE_DIR}/PCR/mperror.c
    ${CMAKE_SOURCE_DIR}/PCR/putattr.c
    ${CMAKE_SOURCE_DIR}/PCR/putx0.c
    ${CMAKE_SOURCE_DIR}/PCR/rcoords.c
    ${CMAKE_SOURCE_DIR}/PCR/rgetrow.c
    ${CMAKE_SOURCE_DIR}/PCR/rmalloc.c
    ${CMAKE_SOURCE_DIR}/PCR/rputrow.c
    ${CMAKE_SOURCE_DIR}/PCR/strpad.c
    ${CMAKE_SOURCE_DIR}/PCR/trackmm.c
    ${CMAKE_SOURCE_DIR}/PCR/vsdef.c
    ${CMAKE_SOURCE_DIR}/PCR/vscatchs.c
)

target_sources(bmi_openlisem_core PRIVATE ${LISEM_CORE_SOURCES})

# Python module using pybind11
pybind11_add_module(bmi_openlisem MODULE
    bmi_python.cpp
)

target_link_libraries(bmi_openlisem PRIVATE
    bmi_openlisem_core
)

target_include_directories(bmi_openlisem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/include
)

# Install targets
install(TARGETS bmi_openlisem_core bmi_openlisem
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

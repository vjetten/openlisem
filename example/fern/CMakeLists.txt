# ADD_SUBDIRECTORY(test)
ADD_EXECUTABLE(add_two_rasters
    ../../CsfMap
    add_two_rasters
)
TARGET_LINK_LIBRARIES(add_two_rasters
    ${LISEM_EXTERNAL_LIBRARIES}
)

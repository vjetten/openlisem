add_subdirectory(data)

set(TEST_NAMES
    array_test
    io_test
    masked_raster_test
    raster_test
)


foreach(NAME ${TEST_NAMES})
    add_unit_test(
        NAME ${NAME}
        LINK_LIBRARIES liblisem ${LISEM_EXTERNAL_LIBRARIES}
    )
endforeach()

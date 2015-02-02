add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/default.map
    COMMAND
        ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/default.map
            ${CMAKE_CURRENT_BINARY_DIR}
)
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/wrong_format.map
    COMMAND
        ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/wrong_format.map
            ${CMAKE_CURRENT_BINARY_DIR}
)
add_custom_target(copy_data ALL
    DEPENDS
        ${CMAKE_CURRENT_BINARY_DIR}/default.map
        ${CMAKE_CURRENT_BINARY_DIR}/wrong_format.map
)

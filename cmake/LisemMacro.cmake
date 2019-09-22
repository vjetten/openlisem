# Add a test target.
# Also configures the environment to point to the location of shared libs.
# The idea of this is to keep the dev's shell as clean as possible. Use
# ctest command to run unit tests.
# NAME : Name of test module, without extension.
# LINK_LIBRARIES: Libraries to link against.
FUNCTION(add_unit_test)
    SET(options "")
    SET(one_value_arguments SCOPE NAME)
    SET(multi_value_arguments LINK_LIBRARIES)

    CMAKE_PARSE_ARGUMENTS(add_unit_test "${options}" "${one_value_arguments}"
        "${multi_value_arguments}" ${ARGN})

    IF(add_unit_test_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR
            "Macro called with unrecognized arguments: "
            "${add_unit_test_UNPARSED_ARGUMENTS}"
        )
    ENDIF()

    SET(test_module_name ${add_unit_test_NAME})
    SET(test_exe_name ${test_module_name})

    ADD_EXECUTABLE(${test_exe_name} ${test_module_name})
    TARGET_LINK_LIBRARIES(${test_exe_name}
        ${add_unit_test_LINK_LIBRARIES}
  #      ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
        stdc++)
    ADD_TEST(NAME ${test_exe_name}
        # catch_system_errors: Prevent UTF to detect system errors. This
        #     messes things up when doing system calls to Python unit tests.
        #     See also: http://lists.boost.org/boost-users/2009/12/55048.php
        COMMAND ${test_exe_name} --catch_system_errors=no)

    # Maybe add ${EXECUTABLE_OUTPUT_PATH} in the future. If needed.
    SET(path_list $ENV{PATH})
  #  LIST(INSERT path_list 0 ${Boost_LIBRARY_DIRS})
    SET(path_string "${path_list}")

    IF(${host_system_name} STREQUAL "windows")
        STRING(REPLACE "\\" "/" path_string "${path_string}")
        STRING(REPLACE ";" "\\;" path_string "${path_string}")
    ELSE()
        STRING(REPLACE ";" ":" path_string "${path_string}")
    ENDIF()

    SET_PROPERTY(TEST ${test_exe_name}
        PROPERTY ENVIRONMENT "PATH=${path_string}")
endfunction()

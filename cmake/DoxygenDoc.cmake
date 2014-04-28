SET(DOXYGEN_TEMPLATE "
    QUIET                   = YES
    WARNINGS                = NO  # Turn on when writing docs.
    WARN_IF_DOC_ERROR       = NO  # Idem.
    WARN_NO_PARAMDOC        = NO  # Idem.
    ALWAYS_DETAILED_SEC     = NO
    INLINE_INHERITED_MEMB   = NO
    INHERIT_DOCS            = YES
    EXTRACT_ALL             = YES
    EXTRACT_PRIVATE         = NO
    MULTILINE_CPP_IS_BRIEF  = YES
    FULL_PATH_NAMES         = YES
    STRIP_FROM_PATH         = ${CMAKE_SOURCE_DIR}
    STRIP_FROM_INC_PATH     = ${CMAKE_SOURCE_DIR}
    EXTRACT_STATIC          = YES
    SOURCE_BROWSER          = YES
    FILE_PATTERNS           = *.h *.hpp *.hxx *.c *.cc *.cpp *.cxx *.dox *.md openlisem.main
    EXCLUDE_PATTERNS        = *Test.h *Test.cc *_test.cc
    EXPAND_ONLY_PREDEF      = NO
    GENERATE_LATEX          = NO
    RECURSIVE               = YES
    SORT_MEMBER_DOCS        = NO
    DOXYFILE_ENCODING       = UTF-8
    USE_MATHJAX             = YES
    # MATHJAX_EXTENSIONS      = TeX/AMSmath TeX/AMSsymbols
")

CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
)

ADD_CUSTOM_TARGET(cpp_doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
)

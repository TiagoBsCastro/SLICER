# FindCCFITS.cmake

if(NOT CCFITS_FOUND)
    find_path(CCFITS_INCLUDE_DIR NAMES CCfits CCfits.h
        HINTS ${CCFITS_ROOT}/include)

    find_library(CCFITS_LIBRARY NAMES CCfits
        HINTS ${CCFITS_ROOT}/lib)

    mark_as_advanced(CCFITS_INCLUDE_DIR CCFITS_LIBRARY)

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(CCFITS DEFAULT_MSG
        CCFITS_LIBRARY CCFITS_INCLUDE_DIR)

    if(CCFITS_FOUND)
        set(CCFITS_INCLUDE_DIRS ${CCFITS_INCLUDE_DIR})
        set(CCFITS_LIBRARIES ${CCFITS_LIBRARY})
    endif()
endif()

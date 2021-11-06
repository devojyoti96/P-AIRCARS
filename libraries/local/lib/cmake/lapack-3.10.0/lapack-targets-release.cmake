#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "blas" for configuration "Release"
set_property(TARGET blas APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(blas PROPERTIES
  IMPORTED_LOCATION_RELEASE "/data1/devojyoti/PhD/P-AIRCARS/libraries/local/lib/libblas.so.3.10.0"
  IMPORTED_SONAME_RELEASE "libblas.so.3"
  )

list(APPEND _IMPORT_CHECK_TARGETS blas )
list(APPEND _IMPORT_CHECK_FILES_FOR_blas "/data1/devojyoti/PhD/P-AIRCARS/libraries/local/lib/libblas.so.3.10.0" )

# Import target "lapack" for configuration "Release"
set_property(TARGET lapack APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(lapack PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "blas"
  IMPORTED_LOCATION_RELEASE "/data1/devojyoti/PhD/P-AIRCARS/libraries/local/lib/liblapack.so.3.10.0"
  IMPORTED_SONAME_RELEASE "liblapack.so.3"
  )

list(APPEND _IMPORT_CHECK_TARGETS lapack )
list(APPEND _IMPORT_CHECK_FILES_FOR_lapack "/data1/devojyoti/PhD/P-AIRCARS/libraries/local/lib/liblapack.so.3.10.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)

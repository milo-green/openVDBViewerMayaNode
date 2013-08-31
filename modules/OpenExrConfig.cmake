# - Try to find OPENEXR
# Once done this will define
#
#  OpenExr_FOUND - system has OpenExr
#  OpenExr_INCLUDE_DIR - the OpenExr include directory


IF( OpenExr_INCLUDE_DIR )
  SET( OpenExr_FIND_QUIETLY TRUE )
ENDIF( OpenExr_INCLUDE_DIR )

FIND_PATH( OpenExr_INCLUDE_DIR OpenEXR/half.h ${OPENEXR_DIR}/include )



# all listed variables are TRUE
INCLUDE( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( OpenExr DEFAULT_MSG OpenExr_INCLUDE_DIR )

MARK_AS_ADVANCED( OpenExr_INCLUDE_DIR )

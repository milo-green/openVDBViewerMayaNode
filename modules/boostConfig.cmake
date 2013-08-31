# - Try to find BOOST INCLUDES
# Once done this will define
#
#  boost_FOUND - system has boost
#  boost_INCLUDE_DIR - the boost include directory




IF( boost_INCLUDE_DIR  )
  SET( boost_FIND_QUIETLY TRUE )
ENDIF( boost_INCLUDE_DIR  )

FIND_PATH( boost_INCLUDE_DIR boost/thread/thread.hpp ${BOOST_DIR}/include )


# handle the QUIETLY and REQUIRED arguments and set boost_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( boost DEFAULT_MSG boost_INCLUDE_DIR )

MARK_AS_ADVANCED( boost_INCLUDE_DIR)
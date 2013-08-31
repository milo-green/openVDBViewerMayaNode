# - Try to find OPENVDB INCLUDES AND LIB
# Once done this will define
#
#  Vdb_FOUND - system has Vdb
#  Vdb_INCLUDE_DIR - the Vdb include directory
#  Vdb_LIBRARY - Link this to use Vdb


IF( Vdb_INCLUDE_DIR AND Vdb_LIBRARY )
  SET( Vdb_FIND_QUIETLY TRUE )
ENDIF( Vdb_INCLUDE_DIR AND Vdb_LIBRARY )

FIND_PATH( Vdb_INCLUDE_DIR openvdb/openvdb.h ${VDB_DIR}/include )

FIND_LIBRARY( Vdb_LIBRARY NAMES openvdb PATHS ${VDB_DIR}/lib )


# all listed variables are TRUE
INCLUDE( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( Vdb DEFAULT_MSG Vdb_LIBRARY Vdb_INCLUDE_DIR )

MARK_AS_ADVANCED( Vdb_INCLUDE_DIR Vdb_LIBRARY )

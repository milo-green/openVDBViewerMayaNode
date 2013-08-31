########### install locations ###########
if( "$ENV{MAYA_PLUGIN_DIR}" STREQUAL "" )

    set( MAYA_PLUGIN_DIR ${CMAKE_INSTALL_PREFIX}/bin/plug-ins )

endif()


if( "$ENV{MAYA_SCRIPT_DIR}" STREQUAL "" )
    set( MAYA_SCRIPT_DIR ${CMAKE_INSTALL_PREFIX}/scripts/others )

endif()


if( CMAKE_SYSTEM_NAME STREQUAL "Darwin" )        # MacOS

    add_definitions( -Dfoo ) # TODO

elseif( CMAKE_SYSTEM_NAME STREQUAL "Linux" )     # Linux

    add_definitions(
        -DBits64_ -D_BOOL -DLINUX -DFUNCPROTO -D_GNU_SOURCE -DLINUX_64
        -DREQUIRE_IOSTREAM -DUNIX -m64 -fPIC -fno-strict-aliasing -O3
        -funsigned-char -ftemplate-depth-25 -pthread
        -Wno-deprecated -Wall -Wno-multichar -Wno-comment
        -Wno-sign-compare -Wno-reorder -fno-gnu-keywords )
    set( CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined" )

elseif( CMAKE_SYSTEM_NAME STREQUAL "Windows" )   # Windows

    add_definitions( -Dfoo ) # TODO

endif()


#find_package( OpenGL REQUIRED ) # seems to work w/o this
# libGLU is provided by Maya

# SETUP MAYA DIRS/LIBS AND INCLUDES -REQUIRED
set( Maya_DIR ../modules )
find_package( Maya REQUIRED )


# SETUP BOOST INCLUDE_DIR-REQUIRED
set( boost_DIR ../modules )
find_package( boost REQUIRED )
set( BOOST_INCLUDE_DIR ${boost_INCLUDE_DIR} )

# SETUP VDB INC_DIRS + LIBS  -REQUIRED
set( Vdb_DIR ../modules )
find_package( Vdb REQUIRED )
set( VDB_LIBRARY ${Vdb_LIBRARY} )
set( VDB_INCLUDE_DIR ${Vdb_INCLUDE_DIR} )

# SETUP OPENEXR INC_ DIRS -REQUIRED
set( OpenExr_DIR ../modules )
find_package( OpenExr REQUIRED )
set( OPENEXR_INCLUDE_DIR ${OpenExr_INCLUDE_DIR} )





include_directories(
    ${Maya_INCLUDE_DIR}
    ${BOOST_INCLUDE_DIR} 
    ${VDB_INCLUDE_DIR}
    ${OPENEXR_INCLUDE_DIR})

########### next target #################

set( vdbViewerNode_LIB_SRCS
    vdbViewerNode.cpp )


add_library( vdbViewerNode MODULE ${vdbViewerNode_LIB_SRCS} )
# omit "lib" prefix from output
set_target_properties( vdbViewerNode PROPERTIES PREFIX "" )

target_link_libraries( vdbViewerNode ${Maya_LIBRARIES} ${VDB_LIBRARY} )

########### install files ###############

install( TARGETS vdbViewerNode
    LIBRARY
        DESTINATION bin/plug-ins )

install( FILES AEvdbViewerNodeTemplate.mel
    DESTINATION scripts/others )



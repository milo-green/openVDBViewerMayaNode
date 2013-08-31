########### install locations ###########

if( "$ENV{MAYA_ICON_DIR}" STREQUAL "" )

    set( MAYA_ICON_DIR ${CMAKE_INSTALL_PREFIX}/icons )

endif()


########### install files ###############

set( vdbViewerNode_LIB_ICONS
    out_vdbViewerNode.xpm vdbViewerNode.xpm )

install( FILES ${vdbViewerNode_LIB_ICONS}
    DESTINATION ${MAYA_ICON_DIR} )

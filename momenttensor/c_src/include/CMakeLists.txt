CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/compearth_config.h)
SET(headers compearth_config.h compearth.h)
FILE(RELATIVE_PATH _relPath "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
INSTALL(FILES ${headers} DESTINATION ${_relPath} COMPONENT headers)

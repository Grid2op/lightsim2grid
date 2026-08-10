list(APPEND SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/BaseBatchSolverSynch.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/YbusPolicy.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/SbusPolicy.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/BaseBatchSweep.cpp"
)
set(SOURCES ${SOURCES} PARENT_SCOPE)

list(APPEND SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/BaseBatchSolverSynch.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/YbusPolicy.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/SbusPolicy.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/BaseBatchSweep.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/BatchAdjoint.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ContinuationSweep.cpp"
)
set(SOURCES ${SOURCES} PARENT_SCOPE)

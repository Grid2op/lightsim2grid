file (RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}")
list(APPEND SOURCES 
    "${_relPath}/BaseBatchSolverSynch.cpp"
    "${_relPath}/TimeSeries.cpp"
    "${_relPath}/ContingencyAnalysis.cpp"
)
set (SOURCES ${SOURCES} PARENT_SCOPE)

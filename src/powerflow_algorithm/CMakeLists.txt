file (RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}")
list(APPEND SOURCES 
    "${_relPath}/GaussSeidelAlgo.cpp"
    "${_relPath}/GaussSeidelSynchAlgo.cpp"
    "${_relPath}/BaseAlgo.cpp"
)
set (SOURCES ${SOURCES} PARENT_SCOPE)

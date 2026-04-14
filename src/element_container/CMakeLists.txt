file (RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}")
list(APPEND SOURCES 
    "${_relPath}/GenericContainer.cpp"
    "${_relPath}/ShuntContainer.cpp"
    "${_relPath}/TrafoContainer.cpp"
    "${_relPath}/LoadContainer.cpp"
    "${_relPath}/GeneratorContainer.cpp"
    "${_relPath}/SGenContainer.cpp"
    "${_relPath}/DCLineContainer.cpp"
    "${_relPath}/LineContainer.cpp"
)
set (SOURCES ${SOURCES} PARENT_SCOPE)

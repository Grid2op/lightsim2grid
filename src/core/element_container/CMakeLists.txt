list(APPEND SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/GenericContainer.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ShuntContainer.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/TrafoContainer.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/LoadContainer.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/GeneratorContainer.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/SGenContainer.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/DCLineContainer.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/LineContainer.cpp"
)
set(SOURCES ${SOURCES} PARENT_SCOPE)

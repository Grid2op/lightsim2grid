file (RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}")

LIST(APPEND 
    SOURCES 
    "${_relPath}/SparseLUSolver.cpp"
)

if(LIB_CKTSO)
    LIST(APPEND 
        SOURCES 
        "${_relPath}/CKTSOSolver.cpp"
    )
endif()

set (SOURCES ${SOURCES} PARENT_SCOPE)

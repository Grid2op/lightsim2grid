file (RELATIVE_PATH _relPath "${PROJECT_SOURCE_DIR}/src" "${CMAKE_CURRENT_SOURCE_DIR}")

LIST(APPEND 
    SOURCES 
    "${_relPath}/SparseLUSolver.cpp"
)

if(KLU_SOLVER_FOUND)
    LIST(APPEND
        SOURCES
        "${_relPath}/KLUSolver.cpp"
    )
endif()

if(LIB_CKTSO)
    LIST(APPEND 
        SOURCES 
        "${_relPath}/CKTSOSolver.cpp"
    )
endif()

if(LIB_NICSLU)
    LIST(APPEND 
        SOURCES 
        "${_relPath}/NICSLUSolver.cpp"
    )
endif()

set (SOURCES ${SOURCES} PARENT_SCOPE)

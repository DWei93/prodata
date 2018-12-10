Set (SRCS 
    Main.cpp )

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wno-unknown-pragmas -DPARALLEL=1 " )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 " )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")
add_executable(prodata ${SRCS})
add_definitions( -DPARALLEL=1)

INCLUDE_DIRECTORIES(../include)    
target_link_libraries(prodata m pthread)
Install(TARGETS prodata
        RUNTIME DESTINATION ${PARADIS_SOURCE_DIR}/bin)


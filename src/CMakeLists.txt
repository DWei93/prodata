Set (SRCS 
    Main.cpp  ReadData.cpp  Util.cpp Math.cpp HandleExtendedDislocation.cpp
    WriteData.cpp Parse.cpp Custom.cpp)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wno-unknown-pragmas -DPARALLEL=1 " )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 " )

if(OPENMP_FOUND)  
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")
endif(OPENMP_FOUND)

add_executable(prodata ${SRCS})
add_definitions( -DPARALLEL=1)

INCLUDE_DIRECTORIES(../include)    
target_link_libraries(prodata m pthread)
Install(TARGETS prodata
        RUNTIME DESTINATION bin)


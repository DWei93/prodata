Set (SRCS 
    Main.cpp  ReadData.cpp  Util.cpp Math.cpp HandleExtendedDislocation.cpp
    WriteData.cpp Parse.cpp GenerateDislocation.cpp 
    ../custom/Custom.cpp)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wno-unknown-pragmas -std=c++11" )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11" )

if(OPENMP_FOUND)     
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")
endif(OPENMP_FOUND)

add_executable(prodata ${SRCS})
add_definitions( -DPARALLEL=1)

INCLUDE_DIRECTORIES(../include $ENV{HOME}/LIB/include)
LINK_DIRECTORIES($ENV{HOME}/LIB/lib)   
If(GSL_FOUND)
    target_link_libraries(prodata m pthread gsl gslcblas)
    add_definitions( -DGSL)
else
    target_link_libraries(prodata m pthread)
endif 
Install(TARGETS prodata
        RUNTIME DESTINATION bin)


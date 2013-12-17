project("DEC")
cmake_minimum_required(VERSION 2.8)

#find_package(AMDIS REQUIRED COMPONENTS umfpack )
find_package(AMDIS REQUIRED)

if(AMDIS_FOUND)
	message("amdis was found\n")
	include(${AMDIS_USE_FILE})
	set(BASIS_LIBS ${AMDIS_LIBRARIES})
endif(AMDIS_FOUND)

FILE(GLOB allsrcs *.cc)
set(dec_src ${allsrcs})
add_library(dec SHARED ${dec_src})
target_link_libraries(dec ${BASIS_LIBS})

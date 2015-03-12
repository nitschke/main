project("DEC")
cmake_minimum_required(VERSION 2.8)

#set(CMAKE_MODULE_PATH /home/ingo/Software/Snapshots/2015-02-19/AMDiS/sequentiell/share/amdis/)

set(AMDIS_DIR /home/ingo/Software/Snapshots/2015-02-19/AMDiS/sequentiell/share/amdis/)

find_package(AMDIS REQUIRED COMPONENTS umfpack )
#find_package(AMDIS REQUIRED)

if(AMDIS_FOUND)
	message("amdis was found\n")
  message(STATUS "Libs found in: " ${AMDIS_LIBRARIES})
	include(${AMDIS_USE_FILE})
	set(BASIS_LIBS ${AMDIS_LIBRARIES})
endif(AMDIS_FOUND)

add_subdirectory(src/dec)
include_directories(src/dec)



add_executable(sphereSharp src/sphereSharp.cc)
target_link_libraries(sphereSharp ${BASIS_LIBS} dec)


project("DEC")
cmake_minimum_required(VERSION 2.8)

#find_package(AMDIS REQUIRED COMPONENTS umfpack )
find_package(AMDIS REQUIRED)

if(AMDIS_FOUND)
	message("amdis was found\n")
	include(${AMDIS_USE_FILE})
	set(BASIS_LIBS ${AMDIS_LIBRARIES})
endif(AMDIS_FOUND)

add_subdirectory(src/dec)
include_directories(src/dec)

add_executable(laplaceSphere src/sphereLaplaceDEC.cc)
target_link_libraries(laplaceSphere dec ${BASIS_LIBS})

add_executable(laplaceTorus src/torusLaplaceDEC.cc)
target_link_libraries(laplaceTorus dec ${BASIS_LIBS})

#create the output dir
#file(MAKE_DIRECTORY output)

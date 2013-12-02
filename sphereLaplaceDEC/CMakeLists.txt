project("DEC")
cmake_minimum_required(VERSION 2.8)

#find_package(AMDIS REQUIRED COMPONENTS umfpack )
find_package(AMDIS REQUIRED)

if(AMDIS_FOUND)
	message("amdis was found\n")
	include(${AMDIS_USE_FILE})
	set(BASIS_LIBS ${AMDIS_LIBRARIES})
endif(AMDIS_FOUND)

file(GLOB SRC src/*.cc)
add_executable(laplaceDEC ${SRC})
target_link_libraries(laplaceDEC ${BASIS_LIBS})

#create the output dir
#file(MAKE_DIRECTORY output)

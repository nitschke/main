project("DEC")
cmake_minimum_required(VERSION 2.8)

#set(CMAKE_MODULE_PATH /home/ingo/Software/Snapshots/2015-02-19/AMDiS/sequentiell/share/amdis/)

set(AMDIS_DIR /home/ingo/Software/Snapshots/2015-02-19/AMDiS/sequentiell/share/amdis/)
#zuhause
#set(AMDIS_DIR /home/ingo/Software/Snapshots/2015-03-06/AMDiS/sequentiell/share/amdis/)

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

add_executable(sphereLaplace src/sphereLaplace.cc)
target_link_libraries(sphereLaplace ${BASIS_LIBS} dec)

add_executable(sphereLaplaceTest src/sphereLaplaceTest.cc)
target_link_libraries(sphereLaplaceTest ${BASIS_LIBS} dec)

add_executable(sphereSimple src/sphereSimple.cc)
target_link_libraries(sphereSimple ${BASIS_LIBS} dec)

add_executable(sphereDeRhamDiffusion src/sphereDeRhamDiffusion.cc)
target_link_libraries(sphereDeRhamDiffusion ${BASIS_LIBS} dec)

add_executable(deRhamDecompSphere src/deRhamDecompSphere.cc)
target_link_libraries(deRhamDecompSphere ${BASIS_LIBS} dec)

add_executable(deRhamDecompTorus src/deRhamDecompTorus.cc)
target_link_libraries(deRhamDecompTorus ${BASIS_LIBS} dec)

add_executable(hodgeDeRhamDecompTorus src/HodgeDeRhamDecompTorus.cc)
target_link_libraries(hodgeDeRhamDecompTorus ${BASIS_LIBS} dec)

add_executable(normSphere src/NormSphere.cc)
target_link_libraries(normSphere ${BASIS_LIBS} dec)

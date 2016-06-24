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

add_executable(ellipsoidLaplaceTest src/ellipsoidLaplaceTest.cc)
target_link_libraries(ellipsoidLaplaceTest ${BASIS_LIBS} dec)

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

add_executable(energySphere src/EnergySphere.cc)
target_link_libraries(energySphere ${BASIS_LIBS} dec)

add_executable(directorFieldSphere src/directorFieldSphere.cc)
target_link_libraries(directorFieldSphere ${BASIS_LIBS} dec)

add_executable(directorFieldEllipsoid src/directorFieldEllipsoid.cc)
target_link_libraries(directorFieldEllipsoid ${BASIS_LIBS} dec)

add_executable(directorFieldQuartic src/directorFieldQuartic.cc)
target_link_libraries(directorFieldQuartic ${BASIS_LIBS} dec)

add_executable(directorFieldOctic src/directorFieldOctic.cc)
target_link_libraries(directorFieldOctic ${BASIS_LIBS} dec)

add_executable(directorFieldTorus src/directorFieldTorus.cc)
target_link_libraries(directorFieldTorus ${BASIS_LIBS} dec)

add_executable(directorField src/directorField.cc)
target_link_libraries(directorField ${BASIS_LIBS} dec)

add_executable(directorFieldOcticCurvTest src/directorFieldOcticCurvTest.cc)
target_link_libraries(directorFieldOcticCurvTest ${BASIS_LIBS} dec)

add_executable(dirCosPhi src/dirCosPhi.cc)
target_link_libraries(dirCosPhi ${BASIS_LIBS} dec)

add_executable(sphereSystemTest src/sphereSystemTest.cc)
target_link_libraries(sphereSystemTest ${BASIS_LIBS} dec)

add_executable(movingSurfaceScalar src/movingSurfaceScalar.cc)
target_link_libraries(movingSurfaceScalar ${BASIS_LIBS} dec)

add_executable(spherePDTest src/spherePDTest.cc)
target_link_libraries(spherePDTest ${BASIS_LIBS} dec)

add_executable(geoPot src/geoPot.cc)
target_link_libraries(geoPot ${BASIS_LIBS} dec)

add_executable(interpolTestSphere src/InterpolTestSphere.cc)
target_link_libraries(interpolTestSphere ${BASIS_LIBS} dec)

add_executable(linMapping src/linMapping.cc)
target_link_libraries(linMapping ${BASIS_LIBS} dec)

add_executable(laplaceConvergenceTest src/LaplaceConvergenceTest.cc)
target_link_libraries(laplaceConvergenceTest ${BASIS_LIBS} dec)

add_executable(directionDerivationTest src/DirectionDerivationTest.cc)
target_link_libraries(directionDerivationTest ${BASIS_LIBS} dec)

add_executable(gradDivProblemTest src/GradDivProblemTest.cc)
target_link_libraries(gradDivProblemTest ${BASIS_LIBS} dec)

add_executable(rotFunTest src/RotFunTest.cc)
target_link_libraries(rotFunTest ${BASIS_LIBS} dec)

add_executable(phiFlatTest src/phiFlatTest.cc)
target_link_libraries(phiFlatTest ${BASIS_LIBS} dec)

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

add_executable(meshCorrector src/correctMesh.cc)
target_link_libraries(meshCorrector ${BASIS_LIBS} dec)

add_executable(laplaceSphere src/sphereLaplaceDEC.cc)
target_link_libraries(laplaceSphere ${BASIS_LIBS} dec)

add_executable(laplaceWeightedSphere src/sphereWeightedLaplace.cc)
target_link_libraries(laplaceWeightedSphere ${BASIS_LIBS} dec)

add_executable(laplaceTorus src/torusLaplaceDEC.cc)
target_link_libraries(laplaceTorus ${BASIS_LIBS} dec)

add_executable(meanCurvTorus src/torusMeanCurvature.cc)
target_link_libraries(meanCurvTorus ${BASIS_LIBS} dec)

add_executable(weingartenTorus src/torusWeingarten.cc)
target_link_libraries(weingartenTorus ${BASIS_LIBS} dec)

add_executable(simpleWeingartenTorus src/torusWeingartenSimple.cc)
target_link_libraries(simpleWeingartenTorus ${BASIS_LIBS} dec)

add_executable(platteTorus src/torusPlatteDEC.cc)
target_link_libraries(platteTorus ${BASIS_LIBS} dec)

add_executable(platteSphere src/spherePlatteDEC.cc)
target_link_libraries(platteSphere ${BASIS_LIBS} dec)

add_executable(heatSphere src/sphereHeat.cc)
target_link_libraries(heatSphere ${BASIS_LIBS} dec)

add_executable(jacobianSphere src/sphereJacobian.cc)
target_link_libraries(jacobianSphere ${BASIS_LIBS} dec)

add_executable(gradSphere src/sphereGrad.cc)
target_link_libraries(gradSphere ${BASIS_LIBS} dec)

add_executable(gradTorus src/torusGrad.cc)
target_link_libraries(gradTorus ${BASIS_LIBS} dec)

add_executable(curvSphere src/sphereCurv.cc)
target_link_libraries(curvSphere ${BASIS_LIBS} dec)

add_executable(defSphere src/defSphere.cc)
target_link_libraries(defSphere ${BASIS_LIBS} dec)

add_executable(heine51b src/heine51b.cc)
target_link_libraries(heine51b ${BASIS_LIBS} dec)

add_executable(heine51c src/heine51c.cc)
target_link_libraries(heine51c ${BASIS_LIBS} dec)

add_executable(heine51cWeingarten src/heine51cWeingarten.cc)
target_link_libraries(heine51cWeingarten ${BASIS_LIBS} dec)

add_executable(heine51bWeingarten src/heine51bWeingarten.cc)
target_link_libraries(heine51bWeingarten ${BASIS_LIBS} dec)

add_executable(sphereWeingarten src/sphereWeingarten.cc)
target_link_libraries(sphereWeingarten ${BASIS_LIBS} dec)

#create the output dir
#create the output dir
#file(MAKE_DIRECTORY output)

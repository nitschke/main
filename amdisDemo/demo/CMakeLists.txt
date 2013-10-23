project("amdis_demo")
cmake_minimum_required(VERSION 2.8)

#find_package(AMDIS REQUIRED COMPONENTS umfpack )
find_package(AMDIS REQUIRED)

if(AMDIS_FOUND)
	message(STATUS "AMDiS was found.")
	include(${AMDIS_USE_FILE})
	SET(BASIS_LIBS ${AMDIS_LIBRARIES})
endif(AMDIS_FOUND)

set(ball src/ball.cc)
set(bunny src/bunny.cc)
set(ellipt src/ellipt.cc)
set(heat src/heat.cc)
set(neumann src/neumann.cc)
set(nonlin src/nonlin.cc)
set(parametric src/parametric.cc)
set(periodic src/periodic.cc)
set(sphere src/sphere.cc)
set(torus src/torus.cc)
set(vecellipt src/vecellipt.cc)
set(vecheat src/vecheat.cc)

add_executable("ball" ${ball})
target_link_libraries("ball" ${BASIS_LIBS})
install(TARGETS ball RUNTIME DESTINATION bin)

add_executable("bunny" ${bunny})
target_link_libraries("bunny" ${BASIS_LIBS})
install(TARGETS bunny RUNTIME DESTINATION bin)

add_executable("ellipt" ${ellipt})
target_link_libraries("ellipt" ${BASIS_LIBS})
install(TARGETS ellipt RUNTIME DESTINATION bin)

add_executable("heat" ${heat})
target_link_libraries("heat" ${BASIS_LIBS})
install(TARGETS heat RUNTIME DESTINATION bin)

add_executable("neumann" ${neumann})
target_link_libraries("neumann" ${BASIS_LIBS})
install(TARGETS neumann RUNTIME DESTINATION bin)

add_executable("nonlin" ${nonlin})
target_link_libraries("nonlin" ${BASIS_LIBS})
install(TARGETS nonlin ball RUNTIME DESTINATION bin)

add_executable("parametric" ${parametric})
target_link_libraries("parametric" ${BASIS_LIBS})
install(TARGETS parametric RUNTIME DESTINATION bin)

add_executable("periodic" ${periodic})
target_link_libraries("periodic" ${BASIS_LIBS})
install(TARGETS periodic RUNTIME DESTINATION bin)

add_executable("sphere" ${sphere})
target_link_libraries("sphere" ${BASIS_LIBS})
install(TARGETS sphere RUNTIME DESTINATION bin)

add_executable("torus" ${torus})
target_link_libraries("torus" ${BASIS_LIBS})
install(TARGETS torus RUNTIME DESTINATION bin)

add_executable("vecellipt" ${vecellipt})
target_link_libraries("vecellipt" ${BASIS_LIBS})
install(TARGETS vecellipt RUNTIME DESTINATION bin)

add_executable("vecheat" ${vecheat})
target_link_libraries("vecheat" ${BASIS_LIBS})
install(TARGETS vecheat RUNTIME DESTINATION bin)

#create the output dir
# file(MAKE_DIRECTORY output)

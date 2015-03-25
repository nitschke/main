#####################################################################################################################
##
## SCRIPT TO INSTALL 
##   - SUITESPARSE
##   - PETSC
##   - SCOTCH
##   - ZOLTAN
##   - AMDIS SEQUENTIELL
##   - AMDIS SEQUENTIELL DEBUG
##   - AMDIS PARALLEL
##   - AMDIS PARALLEL DEBUG
##
## REQUIERED EXTERNAL LIBARIES AND PACKAGES (TODO)
##   - openmpi
##   - cmake
##   - libboost
##   - libblas
##   - liblapack
##   - subversion
##   - gfortran
##   - gcc
##   - build-essential
##   - t.b.d.
##   - bison and flex (for PTScotch within PETSc)
##
#####################################################################################################################

## MACHINE CONFIGURATION VARIABLES
COMP_ON_TAURUS=0
COMP_ON_JUROPA=0
PROCESSORS=5

# BUILD AND DOWNLOAD SUITESPARSE
DOWNLOAD_SUITESPARSE=1
BUILD_SUITESPARES=1

# BUILD AND DOWNLOAD PETSC
DOWNLOAD_PETSC=1
BUILD_PETSC=1

# BUILD AND DOWNLOAD SCOTCH
DOWNLOAD_SCOTCH=0
BUILD_SCOTCH=0
SCOTCH_VERSION="6.0"

# BUILD AND DOWNLOAD ZOLTAN
DOWNLOAD_ZOLTAN=0
BUILD_ZOLTAN=0
ZOLTAN_VERSION="3.81"

# BUILD AND DOWNLOAD AMDIS
DOWNLOAD_AMDIS=1
BUILD_AMDIS_SEQUENTIELL=1
BUILD_AMDIS_SEQUENTIELL_DEBUG=1
BUILD_AMDIS_PARALLEL=1
BUILD_AMDIS_PARALLEL_WITH_DEBUG_SYMBOLS=1
BUILD_AMDIS_PARALLEL_DEBUG=1
ENABLE_ZOLTAN=0
SET_MANUAL_ZOLTAN_FILE_PATHS=0
ZOLTAN_HEADER_FILE=/scratch/wir/software/Zoltan_v3.8/include/zoltan.h
ZOLTAN_HEADER_DIR=/scratch/wir/software/Zoltan_v3.8/include
AMDIS_BRANCH="https://fusionforge.zih.tu-dresden.de/svn/amdis/trunk"
DOWNLOAD_ALTERNATE_REVISION=0                 # CHECK IF OLDER AMDIS REVISION IS NEEDED AND SPECIFY THE FOLLOWING
ALTERNATE_REVISION_NUMBER=2949                # NUMBER OF THE AMDIS REVISION
ALTERNATE_REVISION_DATE="2014-06-16"          # ALTERNATE REVISION DATE. ONLY USED IN FOLDER STRUCTURE (FORMAT: yyyy-mm-dd)

# FOLDER INFORMATION
HOME_FOLDER=$HOME
if [ "$DOWNLOAD_ALTERNATE_REVISION" -eq 1 ] ; then
	DATE=${ALTERNATE_REVISION_DATE}       
else
	DATE=`date +"%Y-%m-%d"`
fi
DESTBASE=${HOME_FOLDER}/Software/Snapshots/${DATE}
BUILDDIR=${HOME_FOLDER}/Software/Snapshots/${DATE}/build
EXTDIR=${BUILDDIR}/extensions

#########################################################################################
### -------------------------------------------------------------------------------------
### DO NOT EDIT AFTER THIS LINE
### -------------------------------------------------------------------------------------
#########################################################################################

if [ "$COMP_ON_JUROPA" -eq 1 ] ; then
	MPI_HEADER_DIR="${MPIHOME}/include"
elif [ "$COMP_ON_TAURUS" -eq 1 ] ; then
	MPIHOME="/opt/mpi/bullxmpi/1.2.4.1"
	MPI_HEADER_DIR="${MPIHOME}/include"
else
	MPI_HEADER_DIR="/usr/include/mpi/"
fi

if [ "$DOWNLOAD_ALTERNATE_REVISION" -eq 1 ] ; then
        AMDIS_BRANCH="-r ${ALTERNATE_REVISION_NUMBER}  ${AMDIS_BRANCH}"
else
        DATE=`date +"%Y-%m-%d"`
fi

## LOAD DEPENDENCIES, IF COMPILATION ON ONE OF THE HPC's
if [ "$COMP_ON_JUROPA" -eq 1 ] ; then
	module unload intel
	module unload mkl
	module load intel/12.1.4
	module load boost
	module load cmake
	export CC="icc"
	export CXX="icpc"
	export CXXFLAGS="${CXXFLAGS} -w"
elif [ "$COMP_ON_TAURUS" -eq 1 ] ; then
	module load cmake
	module load boost
fi

## MAKE BUILD DIRECTORY 
if test -d ${BUILDDIR}; then
  	cd ${BUILDDIR}
else
	mkdir -p ${BUILDDIR}
	cd ${BUILDDIR}
fi

## UNSET PETSC_DIR VARIABLE: THIS IS ONLY NEEDED IF A PRIOR PETSC INSTALLATION EXISTS 
unset PETSC_DIR

## BUILD SUITESPARSE
if [ "$BUILD_SUITESPARES" -eq 1 ] ; then
	echo "BUILD SUITESPARSE"
	if [ "$DOWNLOAD_SUITESPARSE" -eq 1 ] ; then
		wget http://www.cise.ufl.edu/research/sparse/SuiteSparse/SuiteSparse-3.7.1.tar.gz
	fi
	tar xzf SuiteSparse-3.7.1.tar.gz
	cd SuiteSparse
	DEST=${DESTBASE}/SuiteSparse-3.7.1/
	mv UFconfig/UFconfig.mk UFconfig/UFconfig_old.mk
	mkdir -p ${DEST}/lib
	mkdir -p ${DEST}/include
	sed "s:INSTALL_LIB = /usr/local/lib:INSTALL_LIB = ${DEST}/lib:" UFconfig/UFconfig_old.mk > UFconfig/UFconfig.mk
	sed -i "s:INSTALL_INCLUDE = /usr/local/include:INSTALL_INCLUDE =  ${DEST}/include/:" UFconfig/UFconfig.mk
	sed -i "s:# CHOLMOD_CONFIG = -DNPARTITION:CHOLMOD_CONFIG = -DNPARTITION:" UFconfig/UFconfig.mk
	make -j${PROCESSORS} library && make -j${PROCESSORS} install > logfile_make
	if test 0 -eq $?; then
		echo "SUITESPARSE WAS BUILT SUCCESSFULLY"
	else
		echo "COULD NOT BUILD SUITESPARSE"
		echo " > SEE ${BUILDDIR}/SuiteSparse/logfile_make FOR DETAILS"
		exit
	fi
	export SUITESPARSE_LIB="${DEST}/lib/"
	export SUITESPARSE_INC="${DEST}/include/"
	export CHOLMOD_LIBRARIES="${SUITESPARSE_LIB}/libcholmod.a,${SUITESPARSE_LIB}/libcolamd.a"
	export UMFPACK_LIBRARIES="${SUITESPARSE_LIB}/libumfpack.a,${SUITESPARSE_LIB}/libamd.a,${CHOLMOD_LIBRARIES}"
	export PATH=${SUITESPARSE_INC}:$PATH
	cd ..
else
	DEST=${DESTBASE}/SuiteSparse-3.7.1/
	export SUITESPARSE_LIB="${DEST}/lib/"
	export SUITESPARSE_INC="${DEST}/include/"
	export CHOLMOD_LIBRARIES="${SUITESPARSE_LIB}/libcholmod.a,${SUITESPARSE_LIB}/libcolamd.a"
	export UMFPACK_LIBRARIES="${SUITESPARSE_LIB}/libumfpack.a,${SUITESPARSE_LIB}/libamd.a,${CHOLMOD_LIBRARIES}"
	export PATH=${SUITESPARSE_INC}:$PATH
fi

## BUILD PETSC
if [ "$COMP_ON_TAURUS" -eq 1 ] ; then
	module load petsc/3.3-p6
	export PATH="$PATH:$PETSC_DIR/include:$PETSC_LIB"
else
	echo "BUILD PETSC"
	PETSC_VERSION="3.4.3" # 3.3-p1 or 3.4.3
	PETSC_PATCH="" # -p6
	PETSC_STRING="petsc-lite" # petsc-lite
	if [ "$BUILD_PETSC" -eq 1 ] ; then
		if [ "$DOWNLOAD_PETSC" -eq 1 ] ; then
			wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/${PETSC_STRING}-${PETSC_VERSION}${PETSC_PATCH}.tar.gz
		fi
		tar xzf ${PETSC_STRING}-${PETSC_VERSION}${PETSC_PATCH}.tar.gz
		DEST=${DESTBASE}/petsc-${PETSC_VERSION}${PETSC_PATCH}
		cd petsc-${PETSC_VERSION}${PETSC_PATCH}
		ADDITIONAL_CONFIG_FLAGS=""
		
		# DOES NOT WORK -> DON'T KNOW WHY
		#if [ "$ENABLE_ZOLTAN" -eq 1 ] ; then
			#ADDITIONAL_CONFIG_FLAGS="--with-zoltan=1 --download-zoltan=yes" # --download-zoltan=yes" # --with-x=1 --download-x"
		#fi
		
		#./configure --with-pic=1 \
		#	--with-clanguage=C++ \
                #       --prefix=${DEST} \
                #       --with-umfpack=1 \
                #       --with-umfpack-include=$SUITESPARSE_INC \
                #       --with-umfpack-lib=[${UMFPACK_LIBRARIES}] \
                #       --with-blas-lapack-lib=${ACML_LIB}/libacml.so \
                #       --with-metis=1 \
                #       --with-metis-include=${METIS_INC} \
                #       --with-metis-lib=${METIS_LIBRARIES} \
                #       --with-parmetis=1 \
                #       --with-parmetis-include=[${METIS_INC},${PARMETIS_INC}] \
                #       --with-parmetis-lib=[${PARMETIS_LIBRARIES}] \
                #       --with-hypre=1 \
                #       --download-hypre=yes \
                #       --with-superlu=1 \
                #       --download-superlu=yes \
                #       --with-superlu_dist=1 \
                #       --download-superlu_dist=yes \
                #       --with-scalapack=1 \
                #       --download-scalapack=yes \
                #       --with-blacs=1 \
                #       --download-blacs=yes \
                #       --with-mumps=1 \
                #       --download-mumps=yes \
                #       --with-ml=1 \
                #       --download-ml=yes \
                #       --with-cholmod=1 \
                #       --with-cholmod-include=$SUITESPARSE_INC \
                #       --with-cholmod-lib=[${CHOLMOD_LIBRARIES}] \
                #       --with-debugging=0 > logfile_configure
		#./configure --with-pic=1 \
		#	--with-clanguage=C++ \
                #       --prefix=${DEST} \
                #       --with-umfpack=1 \
                #       --with-umfpack-include=$SUITESPARSE_INC \
                #       --with-umfpack-lib=[${UMFPACK_LIBRARIES}] \
                #       --with-blas-lapack-lib=${ACML_LIB}/libacml.so \
                #       --with-metis=1 \
                #       --download-metis=1 \
                #       --with-parmetis=1 \
                #       --download-parmetis=1 \
                #       --with-hypre=1 \
                #       --download-hypre=yes \
                #       --with-superlu=1 \
                #       --download-superlu=yes \
                #       --with-superlu_dist=1 \
                #       --download-superlu_dist=yes \
                #       --with-scalapack=1 \
                #       --download-scalapack=yes \
                #       --with-blacs=1 \
                #       --download-blacs=yes \
                #       --with-mumps=1 \
                #       --download-mumps=yes \
                #       --with-ml=1 \
                #       --download-ml=yes \
                #       --with-cholmod=1 \
                #       --with-cholmod-include=$SUITESPARSE_INC \
                #       --with-cholmod-lib=[${CHOLMOD_LIBRARIES}] \
                #       --with-debugging=0
		if [ "$COMP_ON_JUROPA" -eq 1 ] ; then
			./configure --with-pic=1 \
				--with-clanguage=C++ \
                                --prefix=${DEST} \
                                --with-umfpack=1 \
                                --with-umfpack-include=$SUITESPARSE_INC \
                                --with-umfpack-lib=[${UMFPACK_LIBRARIES}] \
                                --with-ptscotch=1 \
                                --download-ptscotch=yes \
                                --with-metis=1 \
                                --download-metis=1 \
                                --with-parmetis=1 \
                                --download-parmetis=1 \
                                --with-hypre=1 \
                                --download-hypre=yes \
                                --with-superlu=1 \
                                --download-superlu=yes \
                                --with-superlu_dist=1 \
                                --download-superlu_dist=yes \
                                --with-scalapack=1 \
                                --download-scalapack=yes \
                                --with-blacs=1 \
                                --download-blacs=yes \
                                --with-mumps=1 \
                                --download-mumps=yes \
                                --with-ml=1 \
                                --download-ml=yes \
                                --with-cholmod=1 \
                                --with-cholmod-include=$SUITESPARSE_INC \
                                --with-cholmod-lib=[${CHOLMOD_LIBRARIES}] \
                                --with-debugging=0 \
                                ${ADDITIONAL_CONFIG_FLAGS} > logfile_configure
		else
			# PTScotch needs bison and flex to be installed
			./configure --with-pic=1 \
                                --with-clanguage=C++ \
                                --prefix=${DEST} \
                                --with-umfpack=1 \
                                --with-umfpack-include=$SUITESPARSE_INC \
                                --with-umfpack-lib=[${UMFPACK_LIBRARIES}] \
                                --with-ptscotch \
                                --download-ptscotch=yes \
                                --with-metis=1 \
                                --download-metis=1 \
                                --with-parmetis=1 \
                                --download-parmetis=1 \
                                --with-hypre=1 \
                                --download-hypre=yes \
                                --with-superlu=1 \
                                --download-superlu=yes \
                                --with-superlu_dist=1 \
                                --download-superlu_dist=yes \
                                --with-scalapack=1 \
                                --download-scalapack=yes \
                                --with-blacs=1 \
                                --download-blacs=yes \
                                --with-mumps=1 \
                                --download-mumps=yes \
                                --with-ml=1 \
                                --download-ml=yes \
                                --with-cholmod=1 \
                                --with-cholmod-include=$SUITESPARSE_INC \
                                --with-cholmod-lib=[${CHOLMOD_LIBRARIES}] \
                                --with-debugging=0 \
                                ${ADDITIONAL_CONFIG_FLAGS} # > logfile_configure
                                # --download-f2cblaslapack
		fi
		if test 0 -eq $?; then
			echo "PETSC WAS CONFIGURED SUCCESSFULLY"
		else
			echo "COULD NOT CONFIGURE PETSC"
			echo " > SEE ${BUILDDIR}/petsc-${PETSC_VERSION}${PETSC_PATCH}/logfile_configure FOR DETAILS"
			exit 1
		fi
		#make -j${PROCESSORS} PETSC_DIR=`pwd` PETSC_ARCH=arch-linux2-cxx-opt all
		make PETSC_DIR=`pwd` PETSC_ARCH=arch-linux2-cxx-opt all
		if test 0 -eq $?; then
			echo "PETSC WAS BUILT SUCCESSFULLY"
		else
			echo "COULD NOT BUILD PETSC"
			echo " > SEE ${BUILDDIR}/petsc-${PETSC_VERSION}${PETSC_PATCH}/logfile_configure FOR DETAILS"
			exit 1
		fi
		#make -j${PROCESSORS} PETSC_DIR=`pwd` PETSC_ARCH=arch-linux2-cxx-opt install
		make PETSC_DIR=`pwd` PETSC_ARCH=arch-linux2-cxx-opt install
		if test 0 -eq $?; then
			echo "PETSC WAS INSTALLED SUCCESSFULLY"
		else
			echo "COULD NOT INSTALL PETSC"
			exit 1
		fi
		export PETSC_DIR=${DEST}
		cd ..
	else
		export PETSC_DIR=${DESTBASE}/petsc-${PETSC_VERSION}${PETSC_PATCH}
	fi
fi

if [ "$BUILD_SCOTCH" -eq 1 ] ; then
	if [ "$DOWNLOAD_SCOTCH" -eq 1 ] ; then
		svn checkout --username anonsvn --password anonsvn https://scm.gforge.inria.fr/svn/scotch/scotch_${SCOTCH_VERSION}/ scotch_${SCOTCH_VERSION}
	fi
	cd scotch_${SCOTCH_VERSION}/
	DEST=${DESTBASE}/scotch_${SCOTCH_VERSION}
	if test -d ${DEST}; then
		echo "${DEST} already existing!"
	else
		mkdir -p ${DEST}
	fi
	# CONFIGURATION
	cd src/
	if [ "$COMP_ON_JUROPA" -eq 1 ] ; then
		MAKEFILE_INC=Make.inc/Makefile.inc.x86-64_pc_linux2.icc
	elif [ "$COMP_ON_TAURUS" -eq 1 ] ; then
		MAKEFILE_INC=Make.inc/Makefile.inc.x86-64_pc_linux2.icc
	else
		MAKEFILE_INC=Make.inc/Makefile.inc.x86-64_pc_linux2
		# ln -s ${MAKEFILE_INC} Makefile.inc
		cp ${MAKEFILE_INC} Makefile.inc
		sed -i "s:CCD\t\t= gcc:CCD\t\t= gcc -I${MPI_HEADER_DIR}:" Makefile.inc
	fi
	# COMPILATION
	make scotch
	if test 0 -eq $?; then
		echo "SCOTCH WAS BUILT SUCCESSFULLY"
	else
		echo "COULD NOT BUILD SCOTCH"
		exit 1
	fi
	make ptscotch
	if test 0 -eq $?; then
		echo "PTSCOTCH WAS BUILT SUCCESSFULLY"
	else
		echo "COULD NOT BUILD PTSCOTCH"
		exit 1
	fi
	# INSTALL
	make prefix=${DEST} install
	if test 0 -eq $?; then
		echo "SCOTCH/PTSCOTCH WAS INSTALLED SUCCESSFULLY"
	else
		echo "COULD NOT INSTALL SCOTCH/PTSCOTCH"
		exit 1
	fi
	cd ..
	cd ..
	SCOTCH_DIR=${DEST}
else
	SCOTCH_DIR=${DESTBASE}/scotch_${SCOTCH_VERSION}
fi

# if [ "$BUILD_ZOLTAN" -eq 1 ] ; then
# 	if [ "$DOWNLOAD_ZOLTAN" -eq 1 ] ; then
# 		wget http://ftp.mcs.anl.gov/pub/petsc/externalpackages/zoltan_distrib_v3.8.tar.gz
# 	fi
# 	tar xzf zoltan_distrib_v3.8.tar.gz > tar.log
# 	cd Zoltan_v3.8/
# 		mkdir build
# 		cd build
# 			# configure
# 			DEST=${DESTBASE}/Zoltan_v3.8
# 			if [ "$COMP_ON_JUROPA" -eq 1 ] ; then
# 				../configure --prefix=${DEST} \
# 					--enable-mpi \
# 					--with-mpi-compilers=yes \
# 					--with-mpi-libdir=${MPIHOME}/lib \
# 					--with-mpi-incdir=${MPIHOME}/include \
# 					--with-parmetis \
# 					--with-parmetis-incdir=${PETSC_DIR}/include \
# 					--with-parmetis-libdir=${PETSC_DIR}/lib \
# 					CXXFLAGS=-O3 CPPFLAGS=-O3 CFLAGS=-O3
# 			elif [ "$COMP_ON_TAURUS" -eq 1 ] ; then
# 				../configure --prefix=${DEST} \
# 					--enable-mpi \
# 					--with-mpi-compilers=yes \
# 					--with-mpi-libdir=/opt/mpi/bullxmpi/1.2.4.1/lib \
# 					--with-mpi-incdir=/opt/mpi/bullxmpi/1.2.4.1/include \
# 					--with-parmetis \
# 					--with-parmetis-incdir=${PETSC_DIR}/include \
# 					--with-parmetis-libdir=${PETSC_DIR}/lib \
# 					CXXFLAGS=-O3 CPPFLAGS=-O3 CFLAGS=-O3
# 			else
# 				../configure --prefix=${DEST} \
# 					--enable-mpi \
# 					--with-mpi-compilers=yes \
# 					--with-parmetis \
# 					--with-parmetis-incdir=${PETSC_DIR}/include \
# 					--with-parmetis-libdir=${PETSC_DIR}/lib \
# 					CXXFLAGS=-O3 CPPFLAGS=-O3 CFLAGS=-O3
# 			fi
# # 				--with-scotch \
# # 				--with-scotch-libdir=${PETSC_DIR}/lib \
# # 				--with-scotch-incdir=${PETSC_DIR}/include \
# 			if test 0 -eq $?; then
# 				echo "ZOLTAN WAS CONFIGURED SUCCESSFULLY"
# 			else
# 				echo "COULD NOT CONFIGURE ZOLTAN"
# 				exit 1
# 			fi
# 			make everything 
# 			if test 0 -eq $?; then
# 				echo "ZOLTAN WAS BUILT SUCCESSFULLY"
# 			else
# 				echo "COULD NOT BUILD ZOLTAN"
# 				exit 1
# 			fi
# 			make install
# 			if test 0 -eq $?; then
# 				echo "ZOLTAN WAS INSTALLED SUCCESSFULLY"
# 			else
# 				echo "COULD NOT INSTALL ZOLTAN"
# 				exit 1
# 			fi
# 		cd ..
# 	cd ..
# 	ZOLTAN_DIR=${DEST}
# else
# 	ZOLTAN_DIR=${DESTBASE}/Zoltan_v3.8
# fi
if [ "$BUILD_ZOLTAN" -eq 1 ] ; then
	if [ "$DOWNLOAD_ZOLTAN" -eq 1 ] ; then
		wget http://www.cs.sandia.gov/~kddevin/Zoltan_Distributions/zoltan_distrib_v${ZOLTAN_VERSION}.tar.gz
	fi
	tar xzf zoltan_distrib_v${ZOLTAN_VERSION}.tar.gz
	cd Zoltan_v${ZOLTAN_VERSION}/
		mkdir build
		cd build
			# configure
			DEST=${DESTBASE}/Zoltan_v${ZOLTAN_VERSION}
			if [ "$COMP_ON_JUROPA" -eq 1 || "$COMP_ON_TAURUS" -eq 1 ] ; then
				../configure --prefix=${DEST} \
					--enable-mpi \
					--with-mpi-compilers=yes \
					--with-mpi-libdir=${MPIHOME}/lib \
					--with-mpi-incdir=${MPIHOME}/include \
					--with-parmetis \
					--with-parmetis-incdir=${PETSC_DIR}/include \
					--with-parmetis-libdir=${PETSC_DIR}/lib \
					--with-scotch \
					--with-scotch-libdir=${SCOTCH_DIR}/lib \
					--with-scotch-incdir=${SCOTCH_DIR}/include \
					CXXFLAGS=-O3 CPPFLAGS=-O3 CFLAGS=-O3
			else
				../configure --prefix=${DEST} \
					--enable-mpi \
					--with-mpi-compilers=yes \
					--with-parmetis \
					--with-parmetis-incdir=${PETSC_DIR}/include \
					--with-parmetis-libdir=${PETSC_DIR}/lib \
					--with-scotch \
					--with-scotch-libdir=${SCOTCH_DIR}/lib \
					--with-scotch-incdir=${SCOTCH_DIR}/include \
					CXXFLAGS=-O3 CPPFLAGS=-O3 CFLAGS=-O3
			fi
			if test 0 -eq $?; then
				echo "ZOLTAN WAS CONFIGURED SUCCESSFULLY"
			else
				echo "COULD NOT CONFIGURE ZOLTAN"
				exit 1
			fi
			make everything 
			if test 0 -eq $?; then
				echo "ZOLTAN WAS BUILT SUCCESSFULLY"
			else
				echo "COULD NOT BUILD ZOLTAN"
				exit 1
			fi
			make install
			if test 0 -eq $?; then
				echo "ZOLTAN WAS INSTALLED SUCCESSFULLY"
			else
				echo "COULD NOT INSTALL ZOLTAN"
				exit 1
			fi
		cd ..
	cd ..
	ZOLTAN_DIR=${DEST}
else
	ZOLTAN_DIR=${DESTBASE}/Zoltan_v${ZOLTAN_VERSION}
fi
if [ "$SET_MANUAL_ZOLTAN_FILE_PATHS" -eq 0 ] ; then
	ZOLTAN_HEADER_DIR=${ZOLTAN_DIR}/include
	ZOLTAN_HEADER_FILE=${ZOLTAN_HEADER_DIR}/zoltan.h
fi

## BUILD AMDIS 
if [ "$DOWNLOAD_AMDIS" -eq 1 ] ; then
	if [ "$BUILD_AMDIS_SEQUENTIELL" -eq 1 ] || [ "$BUILD_AMDIS_SEQUENTIELL_DEBUG" -eq 1 ] || [ "$BUILD_AMDIS_PARALLEL" -eq 1 ] || ["$BUILD_AMDIS_PARALLEL_DEBUG" -eq 1 ] ; then
		svn co ${AMDIS_BRANCH}/AMDiS amdis
		svn co ${AMDIS_BRANCH}/demo demo
		svn co ${AMDIS_BRANCH}/doc doc
		svn co ${AMDIS_BRANCH}/test test
		svn co ${AMDIS_BRANCH}/extensions extensions
		svn co ${AMDIS_BRANCH}/tools tools
	fi
fi
cd amdis
# build sequentiell AMDiS
if [ "$BUILD_AMDIS_SEQUENTIELL" -eq 1 ] ; then
	AMDIS_VERSION_STRING="AMDIS SEQUENTIAL"
	echo "BUILD ${AMDIS_VERSION_STRING}"
	mkdir sequentiell
	cd sequentiell
	DEST=${DESTBASE}/AMDiS/sequentiell
	cmake -DCMAKE_INSTALL_PREFIX=${DEST} \
		-DCMAKE_BUILD_TYPE=RELEASE \
                -DENABLE_COMPRESSION=ON \
                -DENABLE_UMFPACK=ON \
                -DENABLE_EXTENSIONS=ON \
                -DENABLE_BASE_PROBLEMS=ON \
                -DEXTENSIONS_DIR=${EXTDIR} ..
	make -j${PROCESSORS} install
	if test 0 -eq $?; then
		echo "${AMDIS_VERSION_STRING} WAS COMPILED SUCCESSFULLY"
	else
		echo "COULD NOT COMPILE ${AMDIS_VERSION_STRING}"
		exit 1
	fi
	cd ..
fi
# build sequentiell DebugAMDiS
if [ "$BUILD_AMDIS_SEQUENTIELL_DEBUG" -eq 1 ] ; then
	AMDIS_VERSION_STRING="AMDIS SEQUENTIAL DEBUG"
	echo "BUILD ${AMDIS_VERSION_STRING}"
	mkdir sequentiellDebug
	cd sequentiellDebug
	DEST=${DESTBASE}/AMDiS/sequentiellDebug
	cmake -DCMAKE_BUILD_TYPE=DEBUG \
                -DCMAKE_INSTALL_PREFIX=${DEST} \
                -DENABLE_COMPRESSION=ON \
                -DENABLE_UMFPACK=ON \
                -DENABLE_EXTENSIONS=ON \
                -DENABLE_BASE_PROBLEMS=ON \
                -DEXTENSIONS_DIR=${EXTDIR} ..
	make -j${PROCESSORS} install
	if test 0 -eq $?; then
		echo "${AMDIS_VERSION_STRING} WAS COMPILED SUCCESSFULLY"
	else
		echo "COULD NOT COMPILE ${AMDIS_VERSION_STRING}"
		exit 1
	fi
	cd ..
fi
# build paralell AMDiS
if [ "$BUILD_AMDIS_PARALLEL" -eq 1 ] ; then
	AMDIS_VERSION_STRING="AMDIS PETSC"
	echo "BUILD ${AMDIS_VERSION_STRING}"
	mkdir petsc
	cd petsc
	DEST=${DESTBASE}/AMDiS/petsc
	PATH=$PATH:${PETSC_DIR}/include:${PETSC_DIR}/lib 
	if [ "$ENABLE_ZOLTAN" -eq "1" ]
        then
		cmake -DCMAKE_BUILD_TYPE=RELEASE \
                	-DCMAKE_INSTALL_PREFIX=${DEST} \
                        -DENABLE_COMPRESSION=ON \
                        -DENABLE_UMFPACK=ON \
                        -DENABLE_PARALLEL_DOMAIN=PETSC \
                        -DENABLE_EXTENSIONS=ON \
                        -DENABLE_BASE_PROBLEMS=ON \
                        -DEXTENSIONS_DIR=${EXTDIR} \
                        -DENABLE_ZOLTAN=ON \
                        -DZOLTAN_HEADER_FILE=${ZOLTAN_HEADER_FILE} \
                        -DZOLTAN_HEADER_DIR:PATH=${ZOLTAN_HEADER_DIR} ..
        else
		cmake -DCMAKE_BUILD_TYPE=RELEASE \
                        -DCMAKE_INSTALL_PREFIX=${DEST} \
                        -DENABLE_COMPRESSION=ON \
                        -DENABLE_UMFPACK=ON \
                        -DENABLE_PARALLEL_DOMAIN=PETSC \
                        -DENABLE_EXTENSIONS=ON \
                        -DENABLE_BASE_PROBLEMS=ON \
                        -DEXTENSIONS_DIR=${EXTDIR} ..
	fi  
	make -j${PROCESSORS} install
	if test 0 -eq $?; then
		echo "${AMDIS_VERSION_STRING} WAS COMPILED SUCCESSFULLY"
	else
		echo "COULD NOT COMPILE ${AMDIS_VERSION_STRING}"
		exit 1
	fi
	cd ..
fi
# build paralell AMDiS with debug symbols
if [ "$BUILD_AMDIS_PARALLEL_WITH_DEBUG_SYMBOLS" -eq 1 ] ; then
	AMDIS_VERSION_STRING="AMDIS PETSC WITH DEBUG SYMBOLS"
	echo "BUILD ${AMDIS_VERSION_STRING}"
	mkdir petsc_with_debug_symbols
	cd petsc_with_debug_symbols
	DEST=${DESTBASE}/AMDiS/petsc_with_debug_symbols
	PATH=$PATH:${PETSC_DIR}/include:${PETSC_DIR}/lib
	if [ "$ENABLE_ZOLTAN" -eq "1" ]
        then
		cmake -DCMAKE_BUILD_TYPE=RELWITHDEBINFO \
                        -DCMAKE_INSTALL_PREFIX=${DEST} \
                        -DENABLE_COMPRESSION=ON \
                        -DENABLE_UMFPACK=ON \
                        -DENABLE_PARALLEL_DOMAIN=PETSC \
                        -DENABLE_EXTENSIONS=ON \
                        -DENABLE_BASE_PROBLEMS=ON \
                        -DEXTENSIONS_DIR=${EXTDIR} \
                        -DENABLE_ZOLTAN=ON \
                        -DZOLTAN_HEADER_FILE=${ZOLTAN_HEADER_FILE} \
                        -DZOLTAN_HEADER_DIR:PATH=${ZOLTAN_HEADER_DIR} ..
        else
		cmake -DCMAKE_BUILD_TYPE=RELWITHDEBINFO \
                        -DCMAKE_INSTALL_PREFIX=${DEST} \
                        -DENABLE_COMPRESSION=ON \
                        -DENABLE_UMFPACK=ON \
                        -DENABLE_PARALLEL_DOMAIN=PETSC \
                        -DENABLE_EXTENSIONS=ON \
                        -DENABLE_BASE_PROBLEMS=ON \
                        -DEXTENSIONS_DIR=${EXTDIR} ..
        fi
	make -j${PROCESSORS} install
	if test 0 -eq $?; then
		echo "${AMDIS_VERSION_STRING} WAS COMPILED SUCCESSFULLY"
	else
		echo "COULD NOT COMPILE ${AMDIS_VERSION_STRING}"
		exit 1
	fi
	cd ..
fi
# build parallel DebugAMDiS
if [ "$BUILD_AMDIS_PARALLEL_DEBUG" -eq 1 ] ; then
	AMDIS_VERSION_STRING="AMDIS PETSC DEBUG"
	echo "BUILD ${AMDIS_VERSION_STRING}"
	mkdir petscDebug
	cd petscDebug
	DEST=${DESTBASE}/AMDiS/petscDebug
	PATH=$PATH:${PETSC_DIR}/include:${PETSC_DIR}/lib 
	if [ "$ENABLE_ZOLTAN" -eq "1" ]
        then
		cmake -DCMAKE_BUILD_TYPE=DEBUG \
                        -DCMAKE_INSTALL_PREFIX=${DEST} \
                        -DENABLE_COMPRESSION=ON \
                        -DENABLE_UMFPACK=ON \
                        -DENABLE_PARALLEL_DOMAIN=PETSC \
                        -DENABLE_EXTENSIONS=ON \
                        -DENABLE_BASE_PROBLEMS=ON \
                        -DEXTENSIONS_DIR=${EXTDIR} \
                        -DENABLE_ZOLTAN=ON \
                        -DZOLTAN_HEADER_FILE=${ZOLTAN_HEADER_FILE} \
                        -DZOLTAN_HEADER_DIR:PATH=${ZOLTAN_HEADER_DIR} ..
        else
		cmake -DCMAKE_BUILD_TYPE=DEBUG \
                        -DCMAKE_INSTALL_PREFIX=${DEST} \
                        -DENABLE_COMPRESSION=ON \
                        -DENABLE_UMFPACK=ON \
                        -DENABLE_PARALLEL_DOMAIN=PETSC \
                        -DENABLE_EXTENSIONS=ON \
                        -DENABLE_BASE_PROBLEMS=ON \
                        -DEXTENSIONS_DIR=${EXTDIR} ..
        fi
	make -j${PROCESSORS} install
	if test 0 -eq $?; then
		echo "${AMDIS_VERSION_STRING} WAS COMPILED SUCCESSFULLY"
	else
		echo "COULD NOT COMPILE ${AMDIS_VERSION_STRING}"
		exit 1
	fi
	cd ..
fi

cd .. # return to BUILDDIR

echo "SUCCESSFULLY BUILT:"
if [ "$BUILD_SUITESPARES" -eq 1 ] ; then
	echo " > SUITESPARSE"
fi
if [ "$BUILD_PETSC" -eq 1 ] ; then
	echo " > PETSC"
fi
if [ "$BUILD_SCOTCH" -eq 1 ] ; then
	echo " > SCOTCH"
fi
if [ "$BUILD_ZOLTAN" -eq 1 ] ; then
	echo " > ZOLTAN"
fi
if [ "$BUILD_AMDIS_SEQUENTIELL" -eq 1 ] ; then
	echo " > AMDIS SEQUENTIELL"
fi
if [ "$BUILD_AMDIS_SEQUENTIELL_DEBUG" -eq 1 ] ; then
	echo " > AMDIS SEQUENTIELL DEBUG"
fi
if [ "$BUILD_AMDIS_PARALLEL" -eq 1 ] ; then
	echo " > AMDIS PERALLEL"
fi
if [ "$BUILD_AMDIS_PARALLEL_DEBUG" -eq 1 ] ; then
	echo " > AMDIS PARALLEL DEBUG"
fi
if [ "$BUILD_AMDIS_PARALLEL_WITH_DEBUG_SYMBOLS" -eq 1 ] ; then
	echo " > AMDIS PARALLEL WITH DEBUG SYMBOLS"
fi


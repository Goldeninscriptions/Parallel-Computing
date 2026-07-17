#################################################################
# Try to find PETSc                                             #
#                                                               #
# Once done this will define:                                   #
#  PETSC_FOUND         - system has PETSc                       #
#  PETSC_DIR           - PETSc directory                        #
#  PETSC_ARCH          - PETSc architecture                     #
#  PETSC_INC           - PETSc include directory                #
#  PETSC_LIB           - PETSc library (static or dynamic)      #
#  PETSC_VARIABLES     - Content of PETSc 'petscvariables' file #
#                                                               #
#  PETSC_MUMPS         - Was PETSc compiled with MUMPS?         #
#  PETSC_MUMPS_INC     - PETSc MUMPS include file               #
#  PETSC_MUMPS_LIB     - PETSc MUMPS libraries                  #
#                                                               #
#  PETSC_SCALAPACK     - Was PETSc compiled with ScaLAPACK?     #
#  PETSC_SCALAPACK_LIB - PETSc ScaLAPACK libraries              #
#                                                               #
#  PETSC_PARMETIS      - Was PETSc compiled with ParMETIS?      #
#  PETSC_PARMETIS_LIB  - PETSc ParMETIS libraries               #
#                                                               #
#  PETSC_METIS         - Was PETSc compiled with Metis?         #
#  PETSC_METIS_LIB     - PETSc METIS libraries                  #
#                                                               #
#  PETSC_MPIUNI        - Was PETSc compiled with MPIUNI?        #
#  PETSC_MPIUNI_INC    - MPIUNI include file                    #
#                                                               #
# Usage:                                                        #
#  find_package(PETSc)                                          #
#                                                               #
# Setting these changes the behavior of the search              #
#  PETSC_DIR           - PETSc directory                        #
#  PETSC_ARCH          - PETSc architecture                     #
#################################################################

if(NOT DEFINED PETSC_DIR)
  set(PETSC_DIR $ENV{PETSC_DIR})
endif()
if(NOT DEFINED PETSC_ARCH)
  set(PETSC_ARCH $ENV{PETSC_ARCH})
endif()

set(PETSC_ROOT_DIR "${PETSC_DIR}")
if(PETSC_ARCH AND NOT PETSC_ARCH STREQUAL ".")
  set(PETSC_ARCH_DIR "${PETSC_DIR}/${PETSC_ARCH}")
else()
  set(PETSC_ARCH_DIR "${PETSC_DIR}")
endif()

if(EXISTS "${PETSC_ROOT_DIR}/include" AND EXISTS "${PETSC_ARCH_DIR}/include")
  set(PETSC_INC "${PETSC_ROOT_DIR}/include" "${PETSC_ARCH_DIR}/include")
  list(REMOVE_DUPLICATES PETSC_INC)
else()
  message(SEND_ERROR "PETSc includes not found")
endif()

if(EXISTS "${PETSC_ARCH_DIR}/lib/libpetsc.so")
  set(PETSC_LIB "${PETSC_ARCH_DIR}/lib/libpetsc.so")
elseif(EXISTS "${PETSC_ARCH_DIR}/lib/libpetsc.a")
  set(PETSC_LIB "${PETSC_ARCH_DIR}/lib/libpetsc.a")
elseif(EXISTS "${PETSC_ARCH_DIR}/lib/libpetsc.dylib")
  set(PETSC_LIB "${PETSC_ARCH_DIR}/lib/libpetsc.dylib")
else()
  message(SEND_ERROR "PETSc library not found")
endif()

if(EXISTS ${PETSC_ARCH_DIR}/conf/petscvariables)
  file(STRINGS ${PETSC_ARCH_DIR}/conf/petscvariables PETSC_VARIABLES NEWLINE_CONSUME)
elseif(EXISTS ${PETSC_ARCH_DIR}/lib/petsc/conf/petscvariables)
  file(STRINGS ${PETSC_ARCH_DIR}/lib/petsc/conf/petscvariables PETSC_VARIABLES NEWLINE_CONSUME)
else()
  message(SEND_ERROR "PETSc variables not found")
endif()

if(EXISTS ${PETSC_ARCH_DIR}/lib/libmumps_common.a)
  set(PETSC_MUMPS TRUE)
  set(PETSC_MUMPS_INC ${PETSC_INC})
  set(PETSC_MUMPS_LIB
    ${PETSC_ARCH_DIR}/lib/libcmumps.a
    ${PETSC_ARCH_DIR}/lib/libdmumps.a
    ${PETSC_ARCH_DIR}/lib/libsmumps.a
    ${PETSC_ARCH_DIR}/lib/libzmumps.a
    ${PETSC_ARCH_DIR}/lib/libmumps_common.a
    ${PETSC_ARCH_DIR}/lib/libpord.a)

  if(EXISTS ${PETSC_ARCH_DIR}/lib/libmpiseq.a)
    set(PETSC_MUMPS_SEQ ${PETSC_ARCH_DIR}/lib/libmpiseq.a)
  else()
    set(PETSC_MUMPS_SEQ "")
  endif()
else()
  set(PETSC_MUMPS FALSE)
  set(PETSC_MUMPS_INC "")
  set(PETSC_MUMPS_LIB "")
  set(PETSC_MUMPS_SEQ "")
endif()

if(EXISTS ${PETSC_ARCH_DIR}/lib/libscalapack.a)
  set(PETSC_SCALAPACK TRUE)
  set(PETSC_SCALAPACK_LIB ${PETSC_ARCH_DIR}/lib/libscalapack.a)
else()
  set(PETSC_SCALAPACK FALSE)
  set(PETSC_SCALAPACK_LIB "")
endif()

if(EXISTS ${PETSC_ARCH_DIR}/lib/libparmetis.so)
  set(PETSC_PARMETIS TRUE)
  set(PETSC_PARMETIS_LIB ${PETSC_ARCH_DIR}/lib/libparmetis.so)
elseif(EXISTS ${PETSC_ARCH_DIR}/lib/libparmetis.dylib)
  set(PETSC_PARMETIS TRUE)
  set(PETSC_PARMETIS_LIB ${PETSC_ARCH_DIR}/lib/libparmetis.dylib)
else()
  set(PETSC_PARMETIS FALSE)
  set(PETSC_PARMETIS_LIB "")
endif()

if(EXISTS ${PETSC_ARCH_DIR}/lib/libmetis.so)
  set(PETSC_METIS TRUE)
  set(PETSC_METIS_LIB ${PETSC_ARCH_DIR}/lib/libmetis.so)
elseif(EXISTS ${PETSC_ARCH_DIR}/lib/libmetis.dylib)
  set(PETSC_METIS TRUE)
  set(PETSC_METIS_LIB ${PETSC_ARCH_DIR}/lib/libmetis.dylib)
else()
  set(PETSC_METIS FALSE)
  set(PETSC_METIS_LIB "")
endif()

if(PETSC_VARIABLES)
  string(REGEX MATCH "MPI_IS_MPIUNI = [^\n\r]*" PETSC_MPIUNI ${PETSC_VARIABLES})
else()
  set(PETSC_MPIUNI "")
endif()

if(PETSC_MPIUNI)
  string(REPLACE "MPI_IS_MPIUNI = " "" PETSC_MPIUNI ${PETSC_MPIUNI})
  string(STRIP ${PETSC_MPIUNI} PETSC_MPIUNI)
  string(COMPARE EQUAL ${PETSC_MPIUNI} "1" PETSC_MPIUNI)

  string(REGEX MATCH "MPI_INCLUDE = -I[^\n\r]*" PETSC_MPIUNI_INC ${PETSC_VARIABLES})
  if(PETSC_MPIUNI_INC)
    string(REPLACE "MPI_INCLUDE = -I" "" PETSC_MPIUNI_INC ${PETSC_MPIUNI_INC})
    string(STRIP ${PETSC_MPIUNI_INC} PETSC_MPIUNI_INC)
  else()
    set(PETSC_MPIUNI_INC "${PETSC_DIR}/include/petsc/mpiuni")
  endif()
else()
  set(PETSC_MPIUNI FALSE)
  set(PETSC_MPIUNI_INC "")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSc
  "PETSc could not be found: be sure to set PETSC_DIR in your environment variables"
  PETSC_LIB PETSC_INC PETSC_DIR)

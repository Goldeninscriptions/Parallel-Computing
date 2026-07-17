#############################################################
# Try to find SLEPc                                         #
#                                                           #
# Once done this will define:                               #
#  SLEPC_FOUND     - system has SLEPc                       #
#  SLEPC_DIR       - SLEPc directory                        #
#  SLEPC_INC       - SLEPc include directory                #
#  SLEPC_LIB       - SLEPc library (static or dynamic)      #
#                                                           #
# Usage:                                                    #
#  find_package(SLEPc)                                      #
#                                                           #
# Setting these changes the behavior of the search          #
#  SLEPC_DIR       - SLEPc directory                        #
#############################################################

if(NOT DEFINED SLEPC_DIR)
  set(SLEPC_DIR $ENV{SLEPC_DIR})
endif()

if(PETSC_ARCH AND NOT PETSC_ARCH STREQUAL ".")
  set(SLEPC_ARCH_DIR "${SLEPC_DIR}/${PETSC_ARCH}")
else()
  set(SLEPC_ARCH_DIR "${SLEPC_DIR}")
endif()

if(EXISTS "${SLEPC_DIR}/include" AND EXISTS "${SLEPC_ARCH_DIR}/include")
  set(SLEPC_INC "${SLEPC_DIR}/include" "${SLEPC_ARCH_DIR}/include")
  list(REMOVE_DUPLICATES SLEPC_INC)
else()
  message(SEND_ERROR "SLEPc includes not found")
endif()

if(EXISTS "${SLEPC_ARCH_DIR}/lib/libslepc.so")
  set(SLEPC_LIB "${SLEPC_ARCH_DIR}/lib/libslepc.so")
elseif(EXISTS "${SLEPC_ARCH_DIR}/lib/libslepc.a")
  set(SLEPC_LIB "${SLEPC_ARCH_DIR}/lib/libslepc.a")
elseif(EXISTS "${SLEPC_ARCH_DIR}/lib/libslepc.dylib")
  set(SLEPC_LIB "${SLEPC_ARCH_DIR}/lib/libslepc.dylib")
else()
  message(SEND_ERROR "SLEPc library not found")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SLEPc
  "SLEPc could not be found: be sure to set SLEPC_DIR in your environment variables"
  SLEPC_LIB SLEPC_INC SLEPC_DIR)

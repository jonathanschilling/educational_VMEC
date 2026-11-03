# Write the mgrid file, run VMEC on input.raised_ring in a fresh WORK_DIR and check its wout file with CHECKER.
file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
configure_file("${INPUT_DIR}/input.raised_ring" "${WORK_DIR}/input.raised_ring" COPYONLY)
execute_process(
  COMMAND "${MGRID_WRITER}" "${WORK_DIR}/mgrid_raised_ring.nc"
  RESULT_VARIABLE result)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "cannot write the mgrid file: ${result}")
endif()
execute_process(
  COMMAND "${VMEC}" "input.raised_ring"
  WORKING_DIRECTORY "${WORK_DIR}"
  RESULT_VARIABLE result
  OUTPUT_FILE "${WORK_DIR}/xvmec.log"
  ERROR_FILE "${WORK_DIR}/xvmec.err"
  TIMEOUT 120)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "VMEC failed: ${result}")
endif()
execute_process(
  COMMAND "${CHECKER}" "${WORK_DIR}/wout_raised_ring.nc"
  RESULT_VARIABLE result)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "the raised ring did not descend, see ${WORK_DIR}/xvmec.log")
endif()

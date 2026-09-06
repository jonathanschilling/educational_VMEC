file(MAKE_DIRECTORY "${WORK_DIR}")
configure_file("${INPUT}" "${WORK_DIR}/input.recovery" COPYONLY)
file(REMOVE "${WORK_DIR}/wout_recovery.nc")
execute_process(
  COMMAND "${VMEC}" input.recovery
  WORKING_DIRECTORY "${WORK_DIR}"
  RESULT_VARIABLE result
  OUTPUT_FILE "${WORK_DIR}/solver.log"
  ERROR_FILE "${WORK_DIR}/solver.err"
  TIMEOUT 120)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "VMEC failed: ${result}")
endif()
file(READ "${WORK_DIR}/solver.log" log)
string(REGEX MATCHALL "NS = +3 +NO[.] FOURIER" coarse_stages "${log}")
list(LENGTH coarse_stages coarse_count)
if(NOT coarse_count EQUAL 1)
  message(FATAL_ERROR "Expected one three-surface recovery stage, got ${coarse_count}")
endif()
execute_process(
  COMMAND "${CHECKER}" "${WORK_DIR}/wout_recovery.nc"
  RESULT_VARIABLE result)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "The recovered equilibrium failed validation")
endif()

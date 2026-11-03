# Run VMEC on the two inputs in a fresh WORK_DIR and compare their wout files with CHECKER.
file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
foreach(extension IN ITEMS vertical_shift_0 vertical_shift_1)
  configure_file("${INPUT_DIR}/input.${extension}" "${WORK_DIR}/input.${extension}" COPYONLY)
  execute_process(
    COMMAND "${VMEC}" "input.${extension}"
    WORKING_DIRECTORY "${WORK_DIR}"
    RESULT_VARIABLE result
    OUTPUT_FILE "${WORK_DIR}/${extension}.log"
    ERROR_FILE "${WORK_DIR}/${extension}.err"
    TIMEOUT 120)
  if(NOT result EQUAL 0)
    message(FATAL_ERROR "VMEC failed on input.${extension}: ${result}")
  endif()
endforeach()
execute_process(
  COMMAND "${CHECKER}" "${WORK_DIR}/wout_vertical_shift_0.nc" "${WORK_DIR}/wout_vertical_shift_1.nc"
  RESULT_VARIABLE result)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "a vertical shift changes the equilibrium, see ${WORK_DIR}")
endif()

# Run VMEC on INPUT in a fresh WORK_DIR and check its wout file with CHECKER.
file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
get_filename_component(input_name "${INPUT}" NAME)
string(REGEX REPLACE "^input[.]" "" extension "${input_name}")
configure_file("${INPUT}" "${WORK_DIR}/${input_name}" COPYONLY)
execute_process(
  COMMAND "${VMEC}" "${input_name}"
  WORKING_DIRECTORY "${WORK_DIR}"
  RESULT_VARIABLE result
  OUTPUT_FILE "${WORK_DIR}/xvmec.log"
  ERROR_FILE "${WORK_DIR}/xvmec.err"
  TIMEOUT 120)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "VMEC failed: ${result}")
endif()
execute_process(
  COMMAND "${CHECKER}" "${WORK_DIR}/wout_${extension}.nc"
  RESULT_VARIABLE result)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "wout_${extension}.nc is not a converged equilibrium, see ${WORK_DIR}/xvmec.log")
endif()

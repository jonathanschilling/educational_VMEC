module dbgout
  use json
  implicit none

  character(len=255) :: dump_filename

  contains

function should_write()
  use vmec_main, only: iter2
  use vmec_input, only: max_dump
  implicit none
  logical :: should_write

  should_write = iter2.le.max_dump

end ! function should_write

subroutine open_dbg_context(context_name)
  use vmec_dim,   only: ns
  use vmec_main,  only: iter2
  use vmec_input, only: input_extension
  implicit none

  character(len=*), intent(in) :: context_name

  write(dump_filename, 995) trim(context_name), ns, iter2, trim(input_extension)
995 format(a,'_',i5.5,'_',i6.6,'.',a,'.json')

  call open_dbg_out(dump_filename)

end ! subroutine open_dbg_context

end module dbgout

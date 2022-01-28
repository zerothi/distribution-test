      subroutine die(str)

      character(len=*), intent(in)  :: str
      
      write(6,'(a)') trim(str)
      write(0,'(a)') trim(str)

      stop
      end subroutine die
      

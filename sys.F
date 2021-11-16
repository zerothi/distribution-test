! ---
! Copyright (C) 1996-2016	The SIESTA group
!  This file is distributed under the terms of the
!  GNU General Public License: see COPYING in the top directory
!  or http://www.gnu.org/copyleft/gpl.txt .
! See Docs/Contributors.txt for a list of contributors.
! ---
!!@LICENSE

      module sys
!
!     Termination and messaging routines, MPI aware
!
      implicit none

      public :: die      ! Prints an error message and calls MPI_Abort
      public :: bye      ! Prints an error message and calls MPI_Finalize
      public :: message  ! Prints a message string if node==0

      private

      CONTAINS

      subroutine message(level,str)

      character(len=*), intent(in)  :: level
      ! One of INFO, WARNING, FATAL
      character(len=*), intent(in)  :: str

      integer  ::  lun

         write(6,'(a)') trim(str)
         write(0,'(a)') trim(str)

      end subroutine message
!
!--------------------------------------------------
      subroutine die(str)

#ifdef MPI
      use mpi_siesta
#endif

      character(len=*), intent(in), optional   :: str

      integer  ::  lun
#ifdef MPI
      integer MPIerror
#endif

! Even though formally (in MPI 1.X), only the master node
! can do I/O, in those systems that allow it having each
! node state its complaint can be useful.

!!                                       if (Node.eq.0) then
      if (present(str)) then
         write(6,'(a)') trim(str)
         write(0,'(a)') trim(str)
      endif
!!                                       endif

#ifdef MPI
      call MPI_Abort(MPI_Comm_World,1,MPIerror)
      stop
#else
#endif
      end subroutine die

!---------------------------------------------------------
      subroutine bye(str)

#ifdef MPI
      use mpi_siesta
#endif

      character(len=*), intent(in), optional   :: str

#ifdef MPI
      integer rc
#endif

         if (present(str)) then
            write(6,'(a)') trim(str)
         endif
         write(6,'(a)') 'Requested End of Run. Bye!!'

#ifdef MPI
      call MPI_Finalize(rc)
      stop
#else
      stop
#endif
      end subroutine bye

      end module sys


!>> some auxilary subroutines for time or others

!>> Convert the current wall-time into a real number
!>  The unit is second.

  subroutine now(time_now)

     use wmpi
     use para, only : dp
     implicit none

     integer   :: time_new(8)
     real(dp), intent(inout)  :: time_now

     call Date_and_time(values=time_new)
     time_now= time_new(3)*24d0*3600d0+time_new(5)*3600d0+&
               time_new(6)*60d0+time_new(7)+time_new(8)/1000d0
  
     return
  end subroutine now

  !>  The unit is second.
  subroutine print_time_cost(time_start, time_end, subname)

     use wmpi
     use para, only : dp, stdout, cpuid
     implicit none
     real(dp), intent(in)   :: time_start
     real(dp), intent(in)   :: time_end
     character(*) :: subname

     if (cpuid==0)write(stdout, 101)'Time cost for ', subname, ' is about ', &
        time_end- time_start, ' s' 

     101 format(1x, 3a, f16.3, a)
     return
  end subroutine print_time_cost

  subroutine printallocationinfo(variablename, ierr)
     use para, only : stdout

     implicit none

     character(len=*), intent(in) :: variablename
     integer, intent(in) :: ierr

     if (ierr/=0) then
        write(*, *)"ERROR: no enough memory for ", variablename, '  STAT=',ierr
#if defined (MPI)
        call mpi_finalize(ierr)
#endif
        stop
     endif

     return
  end subroutine printallocationinfo

  subroutine printerrormsg(errormsg)
     use para, only : stdout
     use wmpi, only : cpuid
     implicit none
     character(*), intent(in) :: errormsg
     integer :: ierr

     if (cpuid==0) then
        write(stdout, *) trim(errormsg)
     endif

#if defined (MPI)
        call mpi_finalize(ierr)
#endif
     stop
     return
  end subroutine printerrormsg




  !> print header
  subroutine header
     use wmpi
     use para, only : stdout, cpuid, version
     implicit none
     if (cpuid==0) then
        write(stdout, '(a)') " ----------------------------------------------------------"
        write(stdout, '(a)') " =========================================================="
        write(stdout, '(a)') " =       ===  ============================================="
        write(stdout, '(a)') " =  ====  ==  ============================================="
        write(stdout, '(a)') " =  ====  ==  ============================================="
        write(stdout, '(a)') " =  ====  ==  ======   ===  = ===  = ===  ===   ===  =   =="
        write(stdout, '(a)') " =       ===    ===     ==     ==     ======  =  ==    =  ="
        write(stdout, '(a)') " =  ========  =  ==  =  ==  =  ==  =  ==  ==     ==  ======"
        write(stdout, '(a)') " =  ========  =  ==  =  ==  =  ==  =  ==  ==  =====  ======"
        write(stdout, '(a)') " =  ========  =  ==  =  ==  =  ==  =  ==  ==  =  ==  ======"
        write(stdout, '(a)') " =  ========  =  ===   ===  =  ==  =  ==  ===   ===  ======"
        write(stdout, '(a)') " =========================================================="
        write(stdout, '(a)') "                                                           "
        write(stdout, '(a)') "               Hi! You are now using Phonnier              "
        write(stdout, '(a,a10)') "                  Version ", version                             
        write(stdout, '(a)') "                      I hope I can help!                   "
        write(stdout, '(a)') "                                                           "
        write(stdout, '(a)') "                  Author: Francesc Ballester              "
        write(stdout, '(a)') "               Email: fballestermacia@gmail.com           "
        write(stdout, '(a)') "               TODO: INCLUDE PAPER AND GITHUB              "
        write(stdout, '(a)') "                                                           "
        write(stdout, '(a)') "         Now, let's start by reading your input file...    "
        write(stdout, '(a)') " =========================================================="
        write(stdout, '(a)') "                                                           "
     endif
  end subroutine header

  !> print footer
  subroutine footer
     use wmpi
     use para
     implicit none
     if (cpuid==0) then
        write(stdout, '(2x,a)') ''
        write(stdout, '(2x,a)') ''
        write(stdout, '(2x,a)') "======================================================================="
        write(stdout, '(2x,a)') 'Hurray! I finished running!'
        write(stdout, '(2x,a)') "Check the rest of the output files to see the results."
        write(stdout, '(2x,a)') "If you think everything is correct and you want to use it for any "
        write(stdout, '(2x,a)') "publication, please acknowledge and cite this program as:"
        write(stdout, '(2x,a)') ''
        write(stdout, '(2x,a)') "Loren ipsum dolor sit amet"
        write(stdout, '(2x,a)') "name name name name"
        write(stdout, '(2x,a)') "where and when"
        write(stdout, '(2x,a)') "doi"
        write(stdout, '(2x,a)') ''
        write(stdout, '(2x,a)') "If you saw any bug or anything odd during execution,        "
        write(stdout, '(2x,a)') "report it at fballestermacia@gmail.com or on github                 "
        write(stdout, '(2x,a)') "Thank you for using Phonnier!           "
        write(stdout, '(2x,a)') ''
        write(stdout, '(2x,a)')  "''And then there's quantum, of course.'  The monk sighed. "
        write(stdout, '(2x,a)')  "'There's always bloody quantum.'' - Sir Terry Pratchett, Night Watch"
        write(stdout, '(2x,a)') "======================================================================="
     endif
  end subroutine footer


  !> Sorting arr in ascending order
   subroutine sortheap(n, arr)
      use para, only : dp
      implicit none
      integer, intent(in) :: n
      real(dp), intent(inout) :: arr(n)

      !* local variables
      integer :: i

      do i=n/2, 1, -1
         call sift_down(i, n)
      enddo

      do i=n, 2, -1
         call swap(arr(1), arr(i))
         call sift_down(1, i-1)
      enddo
      contains
      subroutine sift_down(l, r)
         use para, only : dp
         integer, intent(in) :: l, r
         integer :: j, jold
         real(dp) :: a
         a= arr(l)
         jold= l
         j= l+ l

         do while (j<=r)
            if (j<r) then
               if (arr(j)<arr(j+1))j=j+1
            endif
            if (a>= arr(j))exit
            arr(jold)= arr(j)
            jold= j
            j= j+ j
         enddo
         arr(jold)= a
         return
      end subroutine sift_down
   end subroutine sortheap

   !>> swap two real numbers
   subroutine swap(a, b)
      use para, only : dp
      real(dp), intent(inout) :: a
      real(dp), intent(inout) :: b
      real(dp) :: c
      c=a
      a=b
      b=c
      return
   end subroutine swap



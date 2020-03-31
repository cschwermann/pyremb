program Do_test
   implicit none
   integer, parameter :: DP = Selected_real_kind(15, 307)
   real(kind=DP), allocatable :: array(:), array2(:), temp(:)
   real(kind=DP) :: mysum, tstart, tend, mysum2, value
   integer :: i, n, izero

   n=100000000

   allocate( array(1:n) )
   allocate( array2(1:n) )
   allocate( temp(1:n) )

   array(1:n) = 1.0E-5_DP
   array2(1:n) = 1.0E-5_DP

   !! normal do loop
   mysum = 0.0_DP
   call Cpu_time( tstart )
   do i = 1, n
      mysum = mysum + array(i)
   end do
   call Cpu_time( tend )
   write(*,*) "Normal:"
   write(*,*) "Sum: ", mysum, ";  Timing: ", tend - tstart, " seconds"

   !! concurrent do loop
   mysum = 0.0_DP
   call Cpu_time( tstart )
   do concurrent ( i = 1:n )
      mysum = mysum + array(i)
   end do
   call Cpu_time( tend )
   write(*,*) "Concurrent:"
   write(*,*) "Sum: ", mysum, ";  Timing: ", tend - tstart, " seconds"

   !! intrinsic sum
   mysum = 0.0_DP
   call Cpu_time( tstart )
   mysum = Sum( array(1:n) )
   call Cpu_time( tend )
   write(*,*) "Intrinsic:"
   write(*,*) "Sum: ", mysum, ";  Timing: ", tend - tstart, " seconds"

   !! 1 do loop
   mysum = 0.0_DP
   mysum2 = 0.0_DP
   call Cpu_time( tstart )
   do i = 1, n
      mysum = mysum + array(i)
      mysum2 = mysum2 + array(i)
   end do
   call Cpu_time( tend )
   write(*,*) "One Loop:"
   write(*,*) "Sum: ", mysum, mysum2, ";  Timing: ", tend - tstart, " seconds"

   !! 2 do loops
   mysum = 0.0_DP
   mysum2 = 0.0_DP
   call Cpu_time( tstart )
   do i = 1, n
      mysum = mysum + array(i)
   end do
   do i = 1, n
      mysum2 = mysum2 + array(i)
   end do
   call Cpu_time( tend )
   write(*,*) "Two Loops:"
   write(*,*) "Sum: ", mysum, mysum2, ";  Timing: ", tend - tstart, " seconds"

   !! 2 do loops
   mysum = 0.0_DP
   mysum2 = 0.0_DP
   call Cpu_time( tstart )
   mysum = Sum( array(1:n) )
   mysum2 = Sum( array(1:n) )
   call Cpu_time( tend )
   write(*,*) "2* Intrinsic"
   write(*,*) "Sum: ", mysum, mysum2, ";  Timing: ", tend - tstart, " seconds"

   !! Other test: array maths and sum

   !! 1 do loop
   mysum = 0.0_DP
   call Cpu_time( tstart )
   do i = 1, n
      mysum = mysum + array(i) -array2(i)
   end do
   call Cpu_time( tend )
   write(*,*) "One Loop:"
   write(*,*) "Sum: ", mysum, ";  Timing: ", tend - tstart, " seconds"

   !! intrinsic
   mysum = 0.0_DP
   call Cpu_time( tstart )
   mysum = Sum( array(1:n) - array2(1:n) )
   call Cpu_time( tend )
   write(*,*) "Intrinsic Loop:"
   write(*,*) "Sum: ", mysum, ";  Timing: ", tend - tstart, " seconds"

   !! temp array
   mysum = 0.0_DP
   call Cpu_time( tstart )
   temp(1:n) = array(1:n) - array2(1:n)
   do i = 1, n
      mysum = mysum + temp(i)
   end do
   call Cpu_time( tend )
   write(*,*) "Temp array:"
   write(*,*) "Sum: ", mysum, ";  Timing: ", tend - tstart, " seconds"

   !! Other test: specifically what i do in shepard interpolation
   array(1:n) = 2.0E-5_DP
   array(n/2) = 1.0E-5_DP

   !! temp array
   mysum = 0.0_DP
   call Cpu_time( tstart )
   do concurrent ( i = 1:n )
      value = Abs( array(i) - array2(i) )
      temp(i) = value
      if( value <= 1.0E-10_DP ) izero = i
   end do
   do concurrent ( i = 1:izero )
      mysum = mysum + temp(i)
   end do
   call Cpu_time( tend )
   write(*,*) "Temp array:"
   write(*,*) "Sum: ", mysum, izero, ";  Timing: ", tend - tstart, " seconds"

   !! one loop  
   mysum = 0.0_DP
   call Cpu_time( tstart )
   do concurrent ( i = 1:n )
      value = Abs( array(i) - array2(i) )
      if( value <= 1.0E-10_DP ) izero = i
      mysum = mysum + value
   end do
   call Cpu_time( tend )
   write(*,*) "One loop:"
   write(*,*) "Sum: ", mysum, izero, ";  Timing: ", tend - tstart, " seconds"

   !! one loop with exit 
   mysum = 0.0_DP
   call Cpu_time( tstart )
   do i = 1, n
      value = Abs( array(i) - array2(i) )
      if( value <= 1.0E-10_DP ) then
         izero = i
         exit
      end if
      mysum = mysum + value
   end do
   call Cpu_time( tend )
   write(*,*) "One loop with exit:"
   write(*,*) "Sum: ", mysum, izero, ";  Timing: ", tend - tstart, " seconds"

end program Do_test

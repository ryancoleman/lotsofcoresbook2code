       program stress_test
!======================================================================
!  Simple contrived program to see how systems perform an MPI test 
!  intended to test/stress distributed-memory systems.
! 
!  This test repeatedly "smooths" or "filters" a 3-d (gridded) array,
!  starting either from a random or deterministic data array.
!  
!  This makes a good system "stress test", since the arrays can be as 
!  large or small as required, and the job can be run for as long or 
!  as short as required (by varying the number of filtering passes).
!  The 7-point filter around the centre-point in a 3-d lattice structure
!  ariss quite often numerically, and poses a performance challenge for 
!  scalar processors!
!
!  The program reads 8 input variables, from 6 lines of standard input,
!  or redirected from a file.  See "program parameters" below.
!
!======================================================================
!
       implicit none
!
!======================================================================
!  Considerations for MPI version:
!    "Global" indices have subscript _glob
!    "Local" indices have subscript _loc
!
!    Each "local" array arr_in (over local domain) is extended to 
!    contain an extra "layer", 1-point wide, right around the domain.
!    This extra layer, or "shadow zone", or "halo", is used to hold
!    neighbouring array values in the global (physical) domain, but 
!    which are really part of the local domain of a neighbouring MPI
!    process.  Subroutines "shado_xchange_*" are used to fill this
!    "shadow zone" by exchanging boundary data between MPI processes.
!  
!    The extended domain idea also helps with the numerics of filtering
!    since we can treat the physical (global) boundaries as if there 
!    were another point beyond them.
!
!    Just remember that the (global) array arr_in(nx,ny,nz) in the
!    sequential case becomes a set of (local) arrays 
!       arr_in(nx_loc,ny_loc,nz_loc) 
!    in the domain decomposition, which in turn become the "extended"
!    "local" arrays 
!       arx_in(nx_loc+2,ny_loc_2,nz_loc+2)
!
!    The relationship to keep in mind between these is:
!       arr_in(i,j,k)  ==  arx_in(i+1,j+1,k+1)
!
!======================================================================
!
       integer nx_glob, ny_glob, nz_glob, &
     &         nx_loc,  ny_loc,  nz_loc,  &
     &  i_glob, j_glob, k_glob, i_loc, j_loc, k_loc, i,j,k
       integer iseed, icount, maxcount
       integer iofreq, irkz, indx0, itag1
       integer ierror,myrank,nprocs           ! MPI-related variables
       integer nprocx, nprocy, nprocz         ! domain decomposition!
       integer myrankx, myranky, myrankz      ! relative rank in each direction
       integer iloc, jloc, kloc, iparam(9)
       real    rparam(2)
       include 'mpif.h'
       integer istatus(MPI_STATUS_SIZE)
!
       real xran1, arrt(20), twopi,xdist,ydist,zdist
!      real, allocatable :: arr_in(:,:,:), arr_out(:,:,:) ! Physical domains
       real, allocatable :: arx_in(:,:,:), arx_out(:,:,:) ! Extended domains
       real, allocatable :: ardx(:)
       integer, allocatable :: indx(:)
       real valmax, wght, rwght3,rwght4,rwght5,rwght6
!
!
!=====================================================
!      For an MPI program, initialize MPI first!
!=====================================================
       call MPI_Init(ierror)
             if(ierror.ne.MPI_SUCCESS) then
               print *,'Error in MPI_INIT, ierror=',ierror
               go to 999
             endif
!
         call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierror)
             if(ierror.ne.MPI_SUCCESS) then
               print *,'Error in MPI_COMM_SIZE, ierror=',ierror
               go to 999
             endif
!
         call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierror)
             if(ierror.ne.MPI_SUCCESS) then
               print *,'Error from MPI_COMM_RANK, ierror=',ierror
               go to 999
             endif
!
!debug     write(9+myrank,*) 'Myrank ',myrank,' initialized.'
!
!===================================================================
!      Some program parameters: 
!      (Only one process reads these, then broadcasts to the others)
!===================================================================
      if(myrank.eq.0) then
!     read(*,*) maxcount  ! Max. number of filter iterations (problem length)
      read(*,*) iparam(1) ! Max. number of filter iterations (problem length)
!     read(*,*) nx_glob,ny_glob,nz_glob  ! global grid values
      read(*,*) iparam(2), iparam(3), iparam(4)
!     read(*,*) iofreq    ! No. of iterations between writing output
      read(*,*) iparam(5) ! No. of iterations between writing output
!     read(*,*) wght      ! Filter weights (coefficients)
      read(*,*) rparam(1) ! Filter weights (coefficients)
!     read(*,*) valmax    ! Max. value of 3-d field (for scaling only)
      read(*,*) rparam(2) ! Max. value of 3-d field (for scaling only)
!     read(*,*) iseed     ! Zero for deterministic array; otherwise random seed
      read(*,*) iparam(6) ! Zero for deterministic array; otherwise random seed
      read(*,*) iparam(7), iparam(8), iparam(9) ! domain decomposition 
!
!  Some suggested values:
!      100                 ! maxcount, number of filter iterations
!      400, 400, 800       ! global grid values, to use 1 GB memory altogether.
!      10                  ! iofreq, i.e, output every 10 iterations
!      0.05                ! wght, for relatively "light" filtering
!      100.                ! valmax, just for order-1 output values
!      0                   ! iseed, zero for deterministic array, best for checking correctness
!
! --- Start of MPI-specific section ----
!
!  Check that domain decomp is consistent with nprocs:
        if(nprocs .ne. iparam(7)*iparam(8)*iparam(9)) then
          print *,'Error: domain decomposition inconsistent with total', &
     &   ' process count.  Must quit here.'
          go to 999
        endif
      endif
!====================================================================
!     Distribute these basic control parameters to every process...
!====================================================================
      call MPI_Bcast(iparam(1),9,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      call MPI_Bcast(rparam(1),2,MPI_REAL   ,0,MPI_COMM_WORLD,ierror)
             if(ierror.ne.MPI_SUCCESS) then
               print *,'Error from MPI_Bcast, ierror=',ierror
               go to 999
             endif
! Unpack:
      maxcount = iparam(1)
      nx_glob  = iparam(2)
      ny_glob  = iparam(3)
      nz_glob  = iparam(4)
      iofreq   = iparam(5)
      iseed    = iparam(6)
      nprocx   = iparam(7)
      nprocy   = iparam(8)
      nprocz   = iparam(9)
      wght     = rparam(1)
      valmax   = rparam(2)
!
!    Now we need to get nx_loc, etc. from nx_glob, etc.
      nx_loc = nx_glob/nprocx
      ny_loc = ny_glob/nprocy
      nz_loc = nz_glob/nprocz
!
!
!  Check that nprocx, etc. divide evenly into nx_glob, etc.:
      if( nx_glob.ne.nx_loc*nprocx .or. &
     &    ny_glob.ne.ny_loc*nprocy .or. &
     &    nz_glob.ne.nz_loc*nprocz ) then
          print *,'Error: nprocx, nprocy or nprocz dont divide evenly', & 
     & ' into nx, ny or nz.  Must quit!'
          go to 999
      endif
!
!  Now the subtle part of figuring myrankx, myranky, myrankz.
!  These are the process rankings in each of the 3 dimensions.
!  (ranks (x,y,z) are stored in usual "fortran" (i,j,k) order):
        myrankx = mod(myrank,nprocx)
        myranky = mod(myrank,(nprocx*nprocy))/nprocx
        myrankz = myrank/(nprocx*nprocy)
!debug       write(9+myrank,*) 'Myrank, and myrank in x,y,z dimensions are:'
!debug    &       , myrank, myrankx,myranky, myrankz
!
! --- End of MPI-specific section ----
!
!
!  The beauty of MPI is that we need only allocate memory for local vars:
!     allocate (arr_in(nx_loc,ny_loc,nz_loc))
      allocate (arx_in(nx_loc+2,ny_loc+2,nz_loc+2)) ! extended domains for MPI
!     allocate (arr_out(nx_loc,ny_loc,nz_loc))
      allocate (arx_out(nx_loc+2,ny_loc+2,nz_loc+2))  ! extended MPI domain 
      allocate (indx(0:nprocz-1))
      allocate (ardx(nz_loc))
!
      if(myrank.eq.0) then
          print *, 'Maxcount:                ',maxcount
          print *, 'nx_glob,ny_glob,nz_glob: ',nx_glob,ny_glob,nz_glob
          print *, 'Output interval:         ',iofreq
          print *, 'Weight (wght):           ',wght
          print *, 'Value range:             ',valmax
          print *, 'Random seed:             ',iseed
          if(iseed.eq.0) then
             print *, 'Using deterministic initialization'
          else
             print *, 'Using random initialization'
          endif
          print *, ' '
      endif
!
!
        rwght3 = 1.0 - 3.0*wght !Weight factors derived from fundamental weight
        rwght4 = 1.0 - 4.0*wght
        rwght5 = 1.0 - 5.0*wght
        rwght6 = 1.0 - 6.0*wght
!
        twopi = 8.0*atan(1.0)
!
! Start the OpenMP Parallel region, with the main arrays "shared"
! to avoid message-passing for "halo-exchange" within MPI sub-domains.
! Only scalars or small common arrays should be "private".
!
!$OMP PARALLEL PRIVATE(i,j,k,i_loc,j_loc,k_loc,xdist,ydist,zdist, &
!$OMP&                 icount)
!
!=====================================================
!  Initialize the main data array,
!  either deterministically, or with random numbers:
!=====================================================
      if(iseed.eq.0) then
!       do k=1,nz_glob
!$OMP DO SCHEDULE (DYNAMIC,1)
        do k_loc=1,nz_loc
          k = k_loc + nz_loc*myrankz   ! k is index in global domain
          zdist=twopi*float(k-1)/float(nz_glob-1)
!         do j=1,ny
          do j_loc=1,ny_loc
            j = j_loc + ny_loc*myranky   ! j is index in global domain
            ydist=twopi*float(j-1)/float(ny_glob-1)
!           do i=1,nx
            do i_loc =1,nx_loc
              i = i_loc + nx_loc*myrankx   ! i is index in global domain
              xdist=twopi*float(i-1)/float(nx_glob-1)
!       arr_in(i,j,k) = valmax*Cos(7.*xdist)*cos(13.*ydist)*
!    &              cos(11.*zdist)
!       arr_in(i_loc,j_loc,k_loc)=valmax*Cos(7.*xdist)*cos(13.*ydist)*
!    &              cos(11.*zdist)
!  Fill in extended array instead:
        arx_in(i_loc+1,j_loc+1,k_loc+1)=valmax*Cos(7.*xdist)* &
     &              cos(13.*ydist)*cos(11.*zdist)
            enddo
          enddo
        enddo
!OMP END DO
!
      else
!  First, test random number generator....
!$OMP DO SCHEDULE (DYNAMIC,1)
       do i=1,20
         arrt(i) = valmax* xran1(iseed)
!debug         write(9+myrank,*)'myrank, iseed, arrt(i) are:  ',
!debug     &         myrank,iseed,arrt(i)
       enddo
!$OMP END DO 
!       
!$OMP DO SCHEDULE (DYNAMIC,1)
!      do k=1,nz
       do k=2,nz_loc+1
!        do j=1,ny
         do j=2,ny_loc+1
!          do i=1,nx
           do i=2,nx_loc+1
!            arr_in(i,j,k) = valmax*(-1. + 2.0*xran1(iseed))
!            arr_in(i,j,k) = valmax*(-1. + 2.0*xran1(iseed+myrank))
!  Fill in extended array instead:
             arx_in(i,j,k) = valmax*(-1. + 2.0*xran1(iseed+myrank))
           enddo
         enddo
       enddo
!$OMP END DO 
      endif
!
      if(myrank.eq.0) then   ! Only process 0 need write this
!$OMP SINGLE
        if(iseed.eq.0) then
          write(*,*) 'Initial (deterministic) selected values are:'
        else
          write(*,*) 'Initial (random) selected values are:'
        endif
!$OMP END SINGLE 
      endif
!
! This next "write" should use same data in MPI version as in original.
! So, need to collect it, assuming global (i,j)=(2,2) are on 1st proc...
      if(myrankx+myranky.eq.0) then
!$OMP DO
         do irkz=0,nprocz-1
          indx(irkz) = 0  ! indx(irkz) is index of points to be sent
         enddo
!$OMP END DO
!
!$OMP SINGLE
         do k_glob=2,nz_glob-1,nz_glob/7 ! all pts as in sequential case
           irkz = (k_glob - 1)/nz_loc    ! z-rank this point lives in
           indx(irkz) = indx(irkz) + 1   ! One more point for this proc.
           if(myrankz.eq.irkz)  then
             k_loc = 1 + mod(k_glob-1,nz_loc)  ! Local z-index
! First try: ardx(indx(irkz)) = arr_in(2,2,k_loc)
             ardx(indx(irkz)) = arx_in(3,3,k_loc+1) ! Use extended array 
           endif
         enddo
         if(myrankz.eq.0) indx0 = indx(0)+1
!
!
! So we've run through all the global z-points, identified each with a
! process rank, accumulated as necessary, and collected the real values too.
! Now to send/receive appropriately...
!
         itag1=888
         do irkz=1,nprocz-1    ! N.B. Start irkz at rank 1
           if(indx(irkz).ge.1) then  ! Something to be sent and recvd
             if(myrankz.eq.irkz) call MPI_SEND(ardx(1),indx(irkz), &
     &           MPI_REAL,   0,itag1,MPI_COMM_WORLD,ierror)
             if(myrankz.eq.0) then   ! rank 0 only receives here
!            Note that we must receive from a "global" rank, not irkz:
               call MPI_RECV(ardx(indx0),indx(irkz),MPI_REAL, &
     &          irkz*nprocx*nprocy,itag1,MPI_COMM_WORLD,istatus,ierror)
               indx0 = indx0 + indx(irkz)
             endif
           endif
         enddo
!$OMP END SINGLE
      endif  ! (myrankx and myranky being both 0...)
!
!
!  Next section for master process only:
      if(myrank.eq.0) then 
!$OMP SINGLE
!seq    write(*,100) (arr_in(2,2,k),k=2,nz-1,nz/7)
        write(*,100) (ardx(k),k=1,indx0-1)
!
! Optionally, remove the /tmp/ prefix to open file in run-time directory:
        open (unit=12,file='/tmp/stresstest.dat',form='formatted')
        write(12,*) 'The Parameters of this run are:'
        write(12,*) 'Maxcount:        ',maxcount
        write(12,*) 'nx,ny,nz (glob): ',nx_glob,ny_glob,nz_glob
        write(12,*) 'nx,ny,nz (local):',nx_loc,ny_loc,nz_loc
        write(12,*) 'Output interval: ',iofreq
        write(12,*) 'Weight (wght):   ',wght
        write(12,*) 'Value range:     ',valmax
        write(12,*) 'Random seed:     ',iseed
        if(iseed.eq.0) then
           write(12,*) 'Using deterministic initialization'
        else
           write(12,*) 'Using random initialization'
        endif
        write(12,*) ' '
        write(12,*) 'Initial random values are (i,j,k=1,10):'
!  In most cases these points will be on proc 0. If not, more work is needed
        if(nx_loc.ge.10 .and. ny_loc.ge.10 .and. nz_loc.ge.10) then
!         write (12,100) (arr_in(10,10,k),k=1,10)
          write (12,100) (arx_in(11,11,k),k=2,11)  ! Use extended array instead
        else
          write (12,*) 'Need to collect diagnostic points from remote', &
     &  ' processes... not yet implemented'
        endif
        write (12,*) ' '
!$OMP END SINGLE 
      endif
!
!
!
!=======================================================
!     Start of outermost loop over smoothing iterations:
!=======================================================
!
      do icount=0,maxcount
!
! Fill edges of "extended" array arx_in with data from neighbour procs:
! These MPI messages should be single-threaded (involves only "Shared" data):
       call shado_xch_EW(arx_in,nx_loc,ny_loc,nz_loc,  &
     &    myrankx,myranky,myrankz,nprocx,nprocy,nprocz)
       call shado_xch_NS(arx_in,nx_loc,ny_loc,nz_loc,  &
     &    myrankx,myranky,myrankz,nprocx,nprocy,nprocz)
       call shado_xch_TB(arx_in,nx_loc,ny_loc,nz_loc,  &
     &    myrankx,myranky,myrankz,nprocx,nprocy,nprocz)
!
!======================================================================
!  Main body of (3-d) data:
!
!  Use of the "extended" array arx_in with "shadow zone" (needed for 
!  MPI case) allows us to radically simplify the smoothing numerical
!  code over the entire global domain, since all physical surfaces,
!  edges and corners are now "inside" the numberical ones!!
!
!  No need any more for separate loops for separate domain sections.
!======================================================================
!
!
!$OMP DO  SCHEDULE (DYNAMIC,1)
!     do k=2,nz-1
      do k=2,nz_loc+1
!       do j=2,ny-1
        do j=2,ny_loc+1
!         do i=2,nx-1
          do i=2,nx_loc+1
            arx_out(i,j,k) = rwght6*arx_in(i,j,k) + wght*(  &
     &         arx_in(i-1,j,k) + arx_in(i+1,j,k) +  &
     &         arx_in(i,j-1,k) + arx_in(i,j+1,k) +  &
     &         arx_in(i,j,k-1) + arx_in(i,j,k+1) )
          enddo
        enddo
      enddo
!$OMP END DO
!
!
!  Do a 2nd smoothing pass to update arx_in for next iteration:
!  First, need to fill in "edges" of arx_out:
!
       call shado_xch_EW(arx_out,nx_loc,ny_loc,nz_loc,  &
     &    myrankx,myranky,myrankz,nprocx,nprocy,nprocz)
       call shado_xch_NS(arx_out,nx_loc,ny_loc,nz_loc,  &
     &    myrankx,myranky,myrankz,nprocx,nprocy,nprocz)
       call shado_xch_TB(arx_out,nx_loc,ny_loc,nz_loc,  &
     &    myrankx,myranky,myrankz,nprocx,nprocy,nprocz)
!
!
!    Main body of (3-d) data, for 2nd pass:
!$OMP DO SCHEDULE (DYNAMIC,1)
!     do k=2,nz-1
      do k=2,nz_loc+1
!       do j=2,ny-1
        do j=2,ny_loc+1
!         do i=2,nx-1
          do i=2,nx_loc+1
            arx_in(i,j,k) = rwght6*arx_out(i,j,k) + wght*(  &
     &         arx_out(i-1,j,k) + arx_out(i+1,j,k) +  &
     &         arx_out(i,j-1,k) + arx_out(i,j+1,k) +  &
     &         arx_out(i,j,k-1) + arx_out(i,j,k+1) )
          enddo
        enddo
      enddo
!$OMP END DO 
!
!
!
!
! This next "write" should use same data in MPI version as in original.
! So, need to collect it, assuming (i,j)=(2,2) are on 1st proc...
!$OMP SINGLE
      if(myrank.eq.0 .and. icount.eq.0)  &
     &     write(*,*) 'Successively smoothed values are:'
!
!  Only do diagnostic output at specified intervals:
      if(icount.lt.10 .or. mod(icount,iofreq).eq.0)  then
       if(myrankx.eq.0 .and. myranky.eq.0) then
         do irkz=0,nprocz-1
          indx(irkz) = 0  ! indx(irkz) is index of points to be sent
         enddo
!
         do k_glob=2,nz_glob-1,nz_glob/7 ! all pts as in sequential case
           irkz = (k_glob - 1)/nz_loc    ! z-rank this point lives in
           indx(irkz) = indx(irkz) + 1   ! One more point for this proc.
           if(myrankz.eq.irkz)  then
             k_loc = 1 + mod(k_glob-1,nz_loc)  ! Local z-index
! First try: ardx(indx(irkz)) = arr_in(2,2,k_loc)
             ardx(indx(irkz)) = arx_in(3,3,k_loc+1) ! Use extended array
           endif
         enddo
         if(myrankz.eq.0) indx0 = indx(0) + 1
!
! So we've run through all the global z-points, identified each with a
! process rank, accumulated as necessary, and collected the real values too
! Now to send/receive appropriately...
!
         itag1=itag1+1
         do irkz=1,nprocz-1    ! N.B. Start irkz at rank 1
           if(indx(irkz).ge.1) then
             if(myrankz.eq.irkz) call MPI_SEND(ardx(1),indx(irkz),  &
     &           MPI_REAL,   0,itag1,MPI_COMM_WORLD,ierror)
             if(myrankz.eq.0) then   ! rank 0 only receives here
!            Note that we must receive from a "global" rank, not irkz:
               call MPI_RECV(ardx(indx0),indx(irkz),MPI_REAL,  &
     &          irkz*nprocx*nprocy,itag1,MPI_COMM_WORLD,istatus,ierror)
               indx0 = indx0 + indx(irkz)
             endif
           endif
         enddo
       endif  ! (myrankx and myranky being both 0...)
!
        if(myrank.eq.0)   &
!    &       write(*,100) (arr_in(2,2,k),k=2,nz-1,nz/7)
     &       write(*,100) (ardx(k),k=1,indx0-1)
      endif   ! (iteration for diagnostic output)
!$OMP END SINGLE
!
!
!
      enddo    ! (end of outer-most "iteration" loop over icount)
!
!================================================================
!
!
!$OMP SINGLE
      if(myrank.eq.0) then
      write(12,*) 'Final (smoothed?) values are (i,j,k=1,10):'
!  In most cases these points will be on proc 0. If not, more work is neede
        if(nx_loc.ge.10 .and. ny_loc.ge.10 .and. nz_loc.ge.10) then
          write (12,100) (arx_in(11,11,k),k=2,11)
        else
          write (12,*) 'Need to collect diagnostic points from remote', &
     &  ' processes... not yet implemented'
        endif
      endif
!$OMP END SINGLE
!$OMP END PARALLEL
!
      call MPI_FINALIZE(ierror)
!
      stop 'Normal end, max smoothing iterations completed.'
 999  stop 'Program ended prematurely due to fatal error'
 100  format(6e12.4)
      end




     
      function xran1(idum)
!--------------------------------------------------------------------
!     Routine from Numerical Recipes to return a uniform random
!     deviate between 0.0 and 1.0.  Set idum to any negative number
!     to initialize or reinitialize the sequence.
!--------------------------------------------------------------------
      parameter (nn=97)
      parameter (m1=259200, ia1=7141, ic1=54773, rm1=1./m1)
      parameter (m2=134456, ia2=8121, ic2=28411, rm2=1./m2)
      parameter (m3=243000, ia3=4561, ic3=51349)
      real xran1, r(nn)
      save r, ix1,ix2,ix3
!
      data iff /0/
      if (idum.lt.0 .or. iff.eq.0) then
        iff = 1
        ix1 = mod(ic1-idum,m1)     ! seed the first routine
        ix1 = mod(ia1*ix1+ic1,m1)
        ix2 = mod(ix1,m2)          ! and use it to seed the second
        ix1 = mod(ia1*ix1+ic1,m1)
        ix3 = mod(ix1,m3)          ! and to seed the third
!
        do 11,j=1,nn                ! fill the table with sequential
          ix1 = mod(ia1*ix1+ic1,m1) ! random deviates generated by the
          ix2 = mod(ia2*ix2+ic2,m2) ! first two routines
          r(j) = (float(ix1)+float(ix2)*rm2)*rm1
 11     continue
!
        idum = 1
      endif
!
      ix1 = mod(ia1*ix1+ic1,m1)
      ix2 = mod(ia2*ix2+ic2,m2)
      ix3 = mod(ia3*ix3+ic3,m3)
!
      j = 1+(nn*ix3)/m3  ! use 3rd sequence for no. between 1 and 97
!
      if(j.gt.nn.or.j.lt.1) pause
      xran1 = r(j)        ! return that table entry
      r(j) = (float(ix1)+float(ix2)*rm2)*rm1    ! and refill it
      return
      end



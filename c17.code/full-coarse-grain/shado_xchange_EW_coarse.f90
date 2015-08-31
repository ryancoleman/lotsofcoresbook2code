      subroutine shado_xch_EW(arx_in,nx_loc,ny_loc,nz_loc,  &
     &   myrankx,myranky,myrankz,nprocx,nprocy,nprocz, &
     &   num_threads, mythread, nz_tloc)
!====================================================================
!  Routine to do east-west shadow-zone exchange.
!  
!  If this process has a neighbour to east or west, a shadow-zone
!  swap will occur.
!
!  The eastmost surface of arr_in (i.e., the physical local domain)
!  is sent to the westmost surface of arx_in (i.e., the extended local
!  domain) of the process to the east - if there is one.
!  The westmost surface of arr_in (i.e., the physical local domain)
!  is sent to the eastmost surface of arx_in (i.e., the extended
!  local domain) of the process to the west - if there is one.
!
!  Similarly, the westmost surface of arr_in (i.e., the physical local
!  domain) from any process to the east is received into the eastmost
!  surface of arx_in (i.e., the extended local domain).
!  The eastmost surface of arr_in (i.e., the physical local domain)
!  from any process to the west is received into the westmost surface
!  of arx_in (i.e., the extended local domain).
!
!  To help avoid contention, let even processes send eastwards first,
!  while odd processes send westwards first.  That way, all processes
!  should have data ready to receive when MPI_Recv is called.
!
!  If at eastern or western boundaries of the global domain, just
!  copy the edge of arr_in into the edge of arx_in.
!====================================================================
!
      implicit none
      real,intent(inout) :: arx_in(nx_loc+2,ny_loc+2,nz_loc+2)
      integer,intent(in) :: myrankx,myranky,myrankz,nprocx,nprocy,nprocz
      integer,intent(in) :: nx_loc, ny_loc, nz_loc
      integer,intent(in) :: num_threads,mythread,nz_tloc
!
!  These variables will all be thread-local:
      integer ::  j, k, kmin
!orig real shado_E(ny_loc,nz_loc)       !  For data moving eastwards
!orig real shado_W(ny_loc,nz_loc)       !   "    "   "     westwards
      real shado_E(ny_loc,nz_tloc)       !  For data moving eastwards
      real shado_W(ny_loc,nz_tloc)       !   "    "   "     westwards
      logical least,lwest
!
      include 'mpif.h'
      integer istatus(MPI_STATUS_SIZE),ierror,ideste,idestw
!orig integer, parameter :: itage=123, itagw=321
      integer  itage, itagw
!
! 
!debug      print *,'num_threads, mythread, nz_tloc are:', &
!debug     &         num_threads, mythread, nz_tloc 
! First, western global boundary treatment:
!
      if(myrankx.eq.0) then     ! Am at the western boundary already
!$OMP DO 
        do k=2,nz_loc+1         ! Just copy west sfc of arr_in to arx_in
          do j=2,ny_loc+1
            arx_in(1,j,k) = arx_in(2,j,k)  ! simplifies filter 
          enddo
        enddo
!$OMP END DO
      endif
!
!
! Next, eastern global boundary treatment:
!
      if(myrankx.eq.nprocx-1) then   ! Am at the eastern boundary already
!$OMP DO 
        do k=2,nz_loc+1         ! Just copy east sfc of arr_in to arx_in
          do j=2,ny_loc+1
            arx_in(nx_loc+2,j,k) = arx_in(nx_loc+1,j,k)  ! simplifes filter 
          enddo
        enddo
!$OMP END DO
      endif
!
!
! For max. efficiency, let even procs send east first, then west, while
!                           odd procs send west first, then east.
!
      itage=100 + mythread
      itagw=999 - mythread
      least = .false.
      lwest = .false.
        kmin = mythread*nz_tloc ! thread-local range of k-points
      if( mod(myrankx,2).eq.1 ) go to 102
!
! Eastern local boundary treatment:
!
  101 continue
      if(myrankx.lt.nprocx-1) then   ! need to do a shadow-exchange
!
! === Pack:
!orig !$OMP DO 
!orig   do k=2,nz_loc +1              ! Send to east, receive from east
        do k=2,nz_tloc+1              ! Send to east, receive from east
          do j=2,ny_loc+1
!orig       shado_E(j-1,k-1) = arx_in(nx_loc+1,j,     k)
            shado_E(j-1,k-1) = arx_in(nx_loc+1,j,kmin+k)
          enddo
        enddo
!orig !$OMP END DO
!
! Send/recv to/from process "ideste", with same myranky,z and myrankx+1:
          ideste= myrankz*nprocx*nprocy + myranky*nprocx + myrankx + 1
!
!========================================================================
!  Note: it makes sense to post the "receive" here first, since data from 
!  the "odd" procs is already "sent".  Otherwise, if everything is "sent" 
!  before anything is "received", there could be too much data 
!  simultaneously in flight for some MPI implementations...
!========================================================================
!
!orig !$OMP SINGLE
!Use this if MPI not thread-safe:   !$OMP ORDERED
!orig     call MPI_Recv(shado_W(1,1),ny_loc*nz_loc ,MPI_REAL, &
          call MPI_Recv(shado_W(1,1),ny_loc*nz_tloc,MPI_REAL, &
     &                  ideste,itagw,MPI_COMM_WORLD,istatus,ierror)
!
!orig     call MPI_Send(shado_E(1,1),ny_loc*nz_loc, MPI_REAL, &
          call MPI_Send(shado_E(1,1),ny_loc*nz_tloc,MPI_REAL, &
     &                  ideste,itage,MPI_COMM_WORLD,ierror)
!Use this if MPI not thread-safe:   !$OMP END ORDERED
!orig !$OMP END SINGLE
!
! === Unpack:
!orig !$OMP DO 
!orig     do k=2,nz_loc+1
          do k=2,nz_tloc+1
            do j=2,ny_loc+1
!orig         arx_in(nx_loc+2,j,k)      = shado_W(j-1,k-1)
              arx_in(nx_loc+2,j,k+kmin) = shado_W(j-1,k-1)
            enddo
          enddo
!orig !$OMP END DO 
!
      endif
      least = .true. 
      if(lwest) go to 103
!
!
! Western local boundary treatment:
!
  102 continue
      if(myrankx.gt.0) then   ! need to do a shadow-exchange
!
!orig !$OMP DO 
!orig   do k=2,nz_loc+1       ! Send to west, receive from west 
        do k=2,nz_tloc+1       ! Send to west, receive from west 
          do j=2,ny_loc+1
!orig       shado_W(j-1,k-1) = arx_in(2,j,k)  ! Pack up local array values
            shado_W(j-1,k-1) = arx_in(2,j,k+kmin)  ! Pack up local array values
          enddo
        enddo
!orig !$OMP END DO 
!
! Send to process "idestw", with same myranky,z and myrankx-1:
          idestw= myrankz*nprocx*nprocy + myranky*nprocx + myrankx - 1
!
!Use this if MPI not thread-safe:   !$OMP ORDERED
!orig !$OMP SINGLE
!orig     call MPI_Send(shado_W(1,1),ny_loc*nz_loc, MPI_REAL, &
          call MPI_Send(shado_W(1,1),ny_loc*nz_tloc,MPI_REAL, &
     &                  idestw,itagw,MPI_COMM_WORLD,ierror)
!
! Receive from process "idestw" as well:
!orig     call MPI_Recv(shado_E(1,1),ny_loc*nz_loc, MPI_REAL, &
          call MPI_Recv(shado_E(1,1),ny_loc*nz_tloc,MPI_REAL, &
     &                  idestw,itage,MPI_COMM_WORLD,istatus,ierror)
!Use this if MPI not thread-safe:   !$OMP END ORDERED
!orig !$OMP END SINGLE
!
!orig !$OMP DO 
!orig     do k=2,nz_loc +1
          do k=2,nz_tloc+1
            do j=2,ny_loc+1
!orig         arx_in(1,j,k)      = shado_E(j-1,k-1)
              arx_in(1,j,kmin+k) = shado_E(j-1,k-1)
            enddo
          enddo
!orig !$OMP END DO
!
      endif
      lwest = .true.
      if(.not.least) go to 101
!
!
!
  103 return
      end


      subroutine shado_xch_EW(arx_in,nx_loc,ny_loc,nz_loc,
     &   myrankx,myranky,myrankz,nprocx,nprocy,nprocz)
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
      integer myrankx,myranky,myrankz,nprocx,nprocy,nprocz
      integer nx_loc, ny_loc, nz_loc, j, k
      real arx_in(nx_loc+2,ny_loc+2,nz_loc+2)
      real shado_E(ny_loc,nz_loc)       !  For data moving eastwards
      real shado_W(ny_loc,nz_loc)       !   "    "   "     westwards
      logical least,lwest
!
      include 'mpif.h'
      integer istatus(MPI_STATUS_SIZE),ierror,ideste,idestw
      integer, parameter :: itage=123, itagw=321
!
! 
! 
! First, western global boundary treatment:
!
      if(myrankx.eq.0) then     ! Am at the western boundary already
!$OMP PARALLEL DO PRIVATE(j,k)
        do k=2,nz_loc+1         ! Just copy west sfc of arr_in to arx_in
          do j=2,ny_loc+1
            arx_in(1,j,k) = arx_in(2,j,k)  ! simplifies filter 
          enddo
        enddo
!$OMP END PARALLEL DO
      endif
!
!
! Next, eastern global boundary treatment:
!
      if(myrankx.eq.nprocx-1) then   ! Am at the eastern boundary already
!$OMP PARALLEL DO PRIVATE(j,k)
        do k=2,nz_loc+1         ! Just copy east sfc of arr_in to arx_in
          do j=2,ny_loc+1
            arx_in(nx_loc+2,j,k) = arx_in(nx_loc+1,j,k)  ! simplifes filter 
          enddo
        enddo
!$OMP END PARALLEL DO
      endif
!
!
! For max. efficiency, let even procs send east first, then west, while
!                           odd procs send west first, then east.
!
      least = .false.
      lwest = .false.
      if( mod(myrankx,2).eq.1 ) go to 102
!
! Eastern local boundary treatment:
!
  101 continue
      if(myrankx.lt.nprocx-1) then   ! need to do a shadow-exchange
!$OMP PARALLEL DO PRIVATE(j,k)
        do k=2,nz_loc+1              ! Send to east, receive from east
          do j=2,ny_loc+1
            shado_E(j-1,k-1) = arx_in(nx_loc+1,j,k)
          enddo
        enddo
!$OMP END PARALLEL DO
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
          call MPI_Recv(shado_W(1,1),ny_loc*nz_loc,MPI_REAL,
     &                  ideste,itagw,MPI_COMM_WORLD,istatus,ierror)
!
          call MPI_Send(shado_E(1,1),ny_loc*nz_loc,MPI_REAL,
     &                  ideste,itage,MPI_COMM_WORLD,ierror)
!
!$OMP PARALLEL DO PRIVATE(j,k)
          do k=2,nz_loc+1
            do j=2,ny_loc+1
              arx_in(nx_loc+2,j,k) = shado_W(j-1,k-1)
            enddo
          enddo
!$OMP END PARALLEL DO 
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
!$OMP PARALLEL DO PRIVATE(j,k)
        do k=2,nz_loc+1       ! Send to west, receive from west 
          do j=2,ny_loc+1
            shado_W(j-1,k-1) = arx_in(2,j,k)  ! Pack up local array values
          enddo
        enddo
!$OMP END PARALLEL DO 
!
! Send to process "idestw", with same myranky,z and myrankx-1:
          idestw= myrankz*nprocx*nprocy + myranky*nprocx + myrankx - 1
          call MPI_Send(shado_W(1,1),ny_loc*nz_loc,MPI_REAL,
     &                  idestw,itagw,MPI_COMM_WORLD,ierror)
!
! Receive from process "idestw" as well:
          call MPI_Recv(shado_E(1,1),ny_loc*nz_loc,MPI_REAL,
     &                  idestw,itage,MPI_COMM_WORLD,istatus,ierror)
!
!$OMP PARALLEL DO PRIVATE(j,k)
          do k=2,nz_loc+1
            do j=2,ny_loc+1
              arx_in(1,j,k) = shado_E(j-1,k-1)
            enddo
          enddo
!$OMP END PARALLEL DO 
!
      endif
      lwest = .true.
      if(.not.least) go to 101
!
!
!
  103 return
      end


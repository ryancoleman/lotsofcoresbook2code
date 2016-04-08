      subroutine shado_xch_NS(arx_in,nx_loc,ny_loc,nz_loc, &
     &   myrankx,myranky,myrankz,nprocx,nprocy,nprocz, &
     &   num_threads, mythread, nz_tloc) 
!======================================================================
!  Routine to do north-south shadow-zone exchange.
!
!  If this process has a neighbour to north or south, a shadow-zone
!  swap will occur.
!
!  The northmost surface of arr_in (i.e., the physical local domain)
!  is sent to the southmost surface of arx_in (i.e., the extended local 
!  domain) of the process to the north - if there is one.
!  The southmost surface of arr_in (i.e., the physical local domain)
!  is sent to the northmost surface of arx_in (i.e., the extended  
!  local domain) of the process to the south - if there is one.
!
!  Similarly, the southmost surface of arr_in (i.e., the physical local
!  domain) from any process to the north is received into the northmost
!  surface of arx_in (i.e., the extended local domain). 
!  The northmost surface of arr_in (i.e., the physical local domain)
!  from any process to the south is received into the southmost surface
!  of arx_in (i.e., the extended local domain).
!  
!  To help avoid contention, let even processes send northwards first,
!  while odd processes send southwards first.  That way, all processes
!  should have data ready to receive when MPI_Recv is called.
!
!  If at northern or southern boundaries of the global domain, just 
!  copy the edge of arr_in into the edge of arx_in. 
!======================================================================
!
      implicit none
      real,intent(inout) :: arx_in(nx_loc+2,ny_loc+2,nz_loc+2)
      integer,intent(in) :: myrankx,myranky,myrankz,nprocx,nprocy,nprocz
      integer,intent(in) :: nx_loc, ny_loc, nz_loc
      integer,intent(in) :: num_threads,mythread,nz_tloc
!
!  These arrays will all be thread-local:
      integer ::  i, j, k, kmin
!orig real shado_N(nx_loc,nz_loc)       !  For data moving northwards
!orig real shado_S(nx_loc,nz_loc)       !   "    "   "     southwards
      real shado_N(nx_loc,nz_tloc)       !  For data moving northwards
      real shado_S(nx_loc,nz_tloc)       !   "    "   "     southwards
      logical lnorth,lsouth
!
      include 'mpif.h'
      integer istatus(MPI_STATUS_SIZE),ierror,idestn,idests
!orig integer, parameter :: itagn=425, itags=524
      integer  itagn, itags
!
! 
! 
! First, southern global boundary treatment:
!
      if(myranky.eq.0) then     ! Am at the southern boundary already
!$OMP DO
        do k=2,nz_loc+1         ! Just copy south sfc of arr_in to arx_in
          do i=2,nx_loc+1
            arx_in(i,1,k) = arx_in(i,2,k)  ! simplifies filter 
          enddo
        enddo
!$OMP END DO
      endif
!
!
! Next, northern global boundary treatment:
!
      if(myranky.eq.nprocy-1) then   ! Am at the northern boundary already
!$OMP DO 
        do k=2,nz_loc+1         ! Just copy north sfc of arr_in to arx_in
          do i=2,nx_loc+1
            arx_in(i,ny_loc+2,k) = arx_in(i,ny_loc+1,k)  ! simplifes filter 
          enddo
        enddo
!$OMP END DO
      endif
!
!
! For max. efficiency, let even procs send north first, then south, while
!                           odd procs send south first, then north.
!
      itags = 300 + mythread
      itagn = 600 + mythread
      lnorth = .false.
      lsouth = .false.
        kmin = mythread*nz_tloc ! thread-local range of k-points
      if( mod(myranky,2).eq.1 ) go to 102
!
! Northern local boundary treatment:
!
  101 continue
      if(myranky.lt.nprocy-1) then   ! need to do a shadow-exchange
!
!orig !$OMP DO 
!orig   do k=2,nz_loc+1              ! Send to north, receive from north 
        do k=2,nz_tloc+1              ! Send to north, receive from north 
          do i=2,nx_loc+1
!orig       shado_N(i-1,k-1) = arx_in(i,ny_loc+1,k)
            shado_N(i-1,k-1) = arx_in(i,ny_loc+1,kmin+k)
          enddo
        enddo
!orig !$OMP END DO
!
! Send/recv to/from process "idestn", with same myrankx,z and myranky+1:
          idestn= myrankz*nprocx*nprocy + (myranky+1)*nprocx + myrankx
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
!orig     call MPI_Recv(shado_S(1,1),nx_loc*nz_loc ,MPI_REAL,  &
          call MPI_Recv(shado_S(1,1),nx_loc*nz_tloc,MPI_REAL,  &
     &                  idestn,itags,MPI_COMM_WORLD,istatus,ierror)
!
!orig     call MPI_Send(shado_N(1,1),nx_loc*nz_loc ,MPI_REAL,  &
          call MPI_Send(shado_N(1,1),nx_loc*nz_tloc,MPI_REAL,  &
     &                  idestn,itagn,MPI_COMM_WORLD,ierror)
!Use this if MPI not thread-safe:   !$OMP END ORDERED
!orig !$OMP END SINGLE
!
!=== Unpack:
!orig !$OMP DO 
!orig     do k=2,nz_loc +1
          do k=2,nz_tloc+1
            do i=2,nx_loc+1
!orig         arx_in(i,ny_loc+2,k)      = shado_S(i-1,k-1)
              arx_in(i,ny_loc+2,k+kmin) = shado_S(i-1,k-1)
            enddo
          enddo
!orig !$OMP END DO
!
      endif
      lnorth = .true. 
      if(lsouth) go to 103
!
!
! Southern local boundary treatment:
!
  102 continue
      if(myranky.gt.0) then   ! need to do a shadow-exchange
!
!orig !$OMP DO 
!orig   do k=2,nz_loc+1       ! Send to south, receive from south
        do k=2,nz_tloc+1      ! Send to south, receive from south
          do i=2,nx_loc+1
!orig       shado_S(i-1,k-1) = arx_in(i,2,k)  ! Pack up local array values
            shado_S(i-1,k-1) = arx_in(i,2,k+kmin)  ! Pack up local array values
          enddo
        enddo
!orig !$OMP END DO
!
! Send to process "idests", with same myrankx,z and myranky-1:
          idests= myrankz*nprocx*nprocy + (myranky-1)*nprocx + myrankx
!
!orig !$OMP SINGLE
!Use this if MPI not thread-safe:   !$OMP ORDERED
!orig     call MPI_Send(shado_S(1,1),nx_loc*nz_loc ,MPI_REAL,  &
          call MPI_Send(shado_S(1,1),nx_loc*nz_tloc,MPI_REAL,  &
     &                  idests,itags,MPI_COMM_WORLD,ierror)
!
! Receive from process "idests" as well:
!orig     call MPI_Recv(shado_N(1,1),nx_loc*nz_loc ,MPI_REAL,  &
          call MPI_Recv(shado_N(1,1),nx_loc*nz_tloc,MPI_REAL,  &
     &                  idests,itagn,MPI_COMM_WORLD,istatus,ierror)
!Use this if MPI not thread-safe:   !$OMP END ORDERED
!orig !$OMP END SINGLE
!
!orig !$OMP DO 
!orig     do k=2,nz_loc +1
          do k=2,nz_tloc+1
            do i=2,nx_loc+1
!orig         arx_in(i,1,k)      = shado_N(i-1,k-1)
              arx_in(i,1,k+kmin) = shado_N(i-1,k-1)
            enddo
          enddo
!orig !$OMP END DO
!
      endif
      lsouth = .true.
      if(.not.lnorth) go to 101
!
!
!
  103 return
      end


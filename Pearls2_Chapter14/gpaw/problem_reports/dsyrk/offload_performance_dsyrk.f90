!dir$ attributes offload:mic :: mydsyrk
SUBROUTINE MYDSYRK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
  implicit none
  !dir$ attributes offload:mic :: dgemm
!     .. Scalar Arguments ..
  DOUBLE PRECISION ALPHA,BETA
  INTEGER K,LDA,LDC,N
  CHARACTER TRANS,UPLO
!     ..
!     .. Array Arguments ..
  DOUBLE PRECISION A(LDA,*),C(LDC,*)

!  call dgemm(transa, 'n', n, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  call dgemm(TRANS, 't', n, n, k, ALPHA, A, n, A, n, BETA, C, LDC)
end subroutine

program offload_test
  use omp_lib
  implicit none

  !dir$ attributes offload:mic :: dgemm, dsyrk, dsyr2k, mydsyrk
  
  integer, parameter :: dp = SELECTED_REAL_KIND(12)

  real(kind=dp), dimension(:,:), allocatable :: a, b, c
  integer :: n, k, repeats, i
  integer :: lda, ldb, ldb2, ldc

  integer ::pad1, pad2

  real(kind=dp) :: alpha, beta
  character :: transa

  real(kind=dp) :: t0, t1

  character(len=85) :: arg

  n = 256
  k = 64**3
  transa = 'n'
  repeats = 10

  pad1 = 0
  pad2 = 0

  select case(command_argument_count())
  case(0) ! No arguments -> default values
  case(1) ! grid points
     call get_command_argument(1, arg)
     read(arg, *) n
  case(2) 
     call get_command_argument(1, arg)
     read(arg, *) n
     call get_command_argument(2, arg)
     read(arg, *) k
  case(3) 
     call get_command_argument(1, arg)
     read(arg, *) n
     call get_command_argument(2, arg)
     read(arg, *) k
     call get_command_argument(3, arg)
     read(arg, *) transa
  case(5) 
     call get_command_argument(1, arg)
     read(arg, *) n
     call get_command_argument(2, arg)
     read(arg, *) k
     call get_command_argument(3, arg)
     read(arg, *) transa
     call get_command_argument(4, arg)
     read(arg, *) pad1
     call get_command_argument(5, arg)
     read(arg, *) pad2

  end select

  write(*,*) 'Matrix size', n, k, 'transpose ', transa

  allocate(a(n + pad1, k + pad2))
  allocate(b(n + pad1, k + pad2))
  allocate(c(n + pad1, n + pad1))

  a = 0.43
  b = 9.1
  c = 0.0

  alpha = 1.0
  beta = 0.0

  if (transa == 'n') then
     lda = n + pad1
     ldb = k + pad2
     ldb2 = n + pad2
     ldc = n + pad1
  else if (transa == 't') then
     lda = k + pad2
     ldb = k + pad2
     ldb2 = k + pad2
     ldc = n + pad1
  else
     write(*,*) 'incorrect transa: ', transa
     stop
  end if

  write(*,*) 'DGEMM'
! warm up
  do i=1,3
     call dgemm(transa, 'n', n, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  end do
  t0 = omp_get_wtime()
  do i=1, repeats
     call dgemm(transa, 'n', n, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  end do
  t1 = omp_get_wtime()
  write(*,*) '  Time host', (t1 - t0) / repeats

!dir$ offload begin target(mic:0) in(a:align(64)) in(b:align(64)) out(c)
! warm up
  do i=1,3
     call dgemm(transa, 'n', n, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  end do

  t0 = omp_get_wtime()
  do i=1, repeats
     call dgemm(transa, 'n', n, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
  end do
  t1 = omp_get_wtime()
!dir$ end offload
  write(*,*) '  Time target', (t1 - t0) / repeats

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test DSYRK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) 'DSYRK'
! warm up
  do i=1,3
     call dsyrk('u', transa, n, k, alpha, a, lda, beta, c, ldc)
  end do
  t0 = omp_get_wtime()
  do i=1, repeats
     call dsyrk('u', transa, n, k, alpha, a, lda, beta, c, ldc)
  end do
  t1 = omp_get_wtime()
  write(*,*) '  Time host', (t1 - t0) / repeats

!dir$ offload begin target(mic:0) in(a:align(64)) out(c)
! warm up
  do i=1,3
     call dsyrk('u', transa, n, k, alpha, a, lda, beta, c, ldc)
  end do

  t0 = omp_get_wtime()
  do i=1, repeats
     call dsyrk('u', transa, n, k, alpha, a, lda, beta, c, ldc)
  end do
  t1 = omp_get_wtime()
!dir$ end offload
  write(*,*) '  Time target', (t1 - t0) / repeats

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test MYDSYRK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*) 'myDSYRK'
! warm up
  do i=1,3
     call mydsyrk('u', transa, n, k, alpha, a, lda, beta, c, ldc)
  end do
  t0 = omp_get_wtime()
  do i=1, repeats
     call mydsyrk('u', transa, n, k, alpha, a, lda, beta, c, ldc)
  end do
  t1 = omp_get_wtime()
  write(*,*) '  Time host', (t1 - t0) / repeats

!dir$ offload begin target(mic:0) in(a:align(64)) out(c)
! warm up
  do i=1,3
     call mydsyrk('u', transa, n, k, alpha, a, lda, beta, c, ldc)
  end do

  t0 = omp_get_wtime()
  do i=1, repeats
     call mydsyrk('u', transa, n, k, alpha, a, lda, beta, c, ldc)
  end do
  t1 = omp_get_wtime()
!dir$ end offload
  write(*,*) '  Time target', (t1 - t0) / repeats

  deallocate(a)
  deallocate(b)
  deallocate(c)
end program offload_test

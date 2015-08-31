      module points

      type :: point
      double precision c0,c1
      end type point

      integer :: steps=512
      integer :: ARRAY_SIZE=2048
      end module points

      module time_info
          interface
                integer (C_LONG) function get_current_time() &
                        BIND(C, name='get_current_time')
                        use, intrinsic :: ISO_C_BINDING
                        implicit none
                end function get_current_time
          end interface
      end module time_info

      program main
      use points
      use time_info

      implicit none
      double precision, external :: func
      type(point), allocatable :: vals(:)
      double precision, allocatable :: src(:), dst(:)

      interface
      double precision function Interpolate(x,vals,nvals) bind(C)
!DEC$ IF DEFINED(_OPENMP)
      !$omp declare simd (Interpolate) uniform(vals,nvals)
!DEC$ ENDIF
      use points
      double precision, intent(in),value :: x
      type(point), intent(in), dimension(nvals) :: vals
      end function Interpolate
      end interface

      double precision :: prev_val = 0.;
      integer :: ind = 0;
      integer :: i = 0;
      double precision :: x, val, c0, c1, num, ref_val
      integer*8 :: start_time = 0
      integer*8 :: end_time = 0
      integer :: iter = 0
      integer :: ITER_COUNT = 5000

      ! Set interpolation array
      double precision :: delta, delta_inv
      delta = 1./steps
      delta_inv = 1./delta

      allocate(vals(steps+1))
      allocate(src(ARRAY_SIZE))
      allocate(dst(ARRAY_SIZE))
      !dir$ attributes align:64 :: vals, src,dst

      vals(1)%c0 = 0.
      vals(1)%c1 = func(vals(1)%c0)

      prev_val = vals(1)%c1;
      do ind = 1, steps, 1
         x = real(ind) * delta
         val = func(x)
         c0 = (val - prev_val)*delta_inv
         c1 = val - c0*x
         vals(ind+1)%c0 = c0
         vals(ind+1)%c1 = c1
         prev_val = val
      end do

!     Initialize input array
      call RANDOM_NUMBER(src)
      
      start_time = get_current_time()
      do iter=1, ITER_COUNT, 1
!DEC$ IF DEFINED(_OPENMP)
         !$omp simd
!DEC$ ENDIF
         do i=1, ARRAY_SIZE, 1
            dst(i) = Interpolate(src(i),vals,steps+1)
         end do
      end do
      end_time = get_current_time()

!     Test results 
      do i= 1, ARRAY_SIZE, 1 
         ref_val = func(src(i))
         
         if ( (ref_val-dst(i))/ref_val > 0.01 )  then
            print *,"Error Found"
            stop(-1)
         end if 
      end do 
       
      print *, "Test completed in ", (end_time - start_time)/ITER_COUNT, " nsec"

      end program main
      
      double precision function func(x)
      double precision, intent(in) :: x
      func = dexp(x)
      end function func

      integer function FindPosition(x)
!DEC$ IF DEFINED(_OPENMP)
!DEC$ IF DEFINED(__MIC__)
      !$omp declare simd(FindPosition) simdlen(8)
!DEC$ ELSE 
      !$omp declare simd(FindPosition) simdlen(4)
!DEC$ ENDIF
!DEC$ ENDIF
      use points
      double precision, intent(in) :: x
      FindPosition = dlog(dexp(x*steps))+1
      end function FindPosition

      double precision function Interpolate(x,vals,nvals) bind(C)
!DEC$ IF DEFINED(_OPENMP)
      !$omp declare simd(Interpolate) uniform(vals,nvals)
!DEC$ ENDIF
      use points
      double precision, intent(in), value :: x
      integer, intent(in) :: nvals
      type(point), intent(in) :: vals(nvals)
      integer FindPosition
      integer :: ind
      ind = FindPosition(x)
      Interpolate = dlog(dexp( vals(ind)%c0 * x + vals(ind)%c1 ) )
      end function Interpolate

  program gsfcgce_driver

    USE module_mp_gsfcgce
    USE module_model_constants

    implicit none			

    logical, parameter :: write_output = .true.
    integer, parameter :: no_runs = 10

    integer ids, ide, jds, jde, kds, kde ! domain dims (unused)
    integer ims, ime, jms, jme, kms, kme ! memory dims
    integer its, ite, jts, jte, kts, kte ! tile dims

    real dt
    integer itimestep, hail, ice2, isee_comp
    logical f_qg

    real, dimension(:, :), allocatable :: ht, rainnc, rainncv, snownc, snowncv, sr, graupelnc, graupelncv

    real, dimension(:, :, :), allocatable :: th
    real, dimension(:, :, :), allocatable :: qv_curr
    real, dimension(:, :, :), allocatable :: qc_curr
    real, dimension(:, :, :), allocatable :: qr_curr
    real, dimension(:, :, :), allocatable :: qi_curr
    real, dimension(:, :, :), allocatable :: qs_curr
    real, dimension(:, :, :), allocatable :: rho
    real, dimension(:, :, :), allocatable :: pi_phy
    real, dimension(:, :, :), allocatable :: p
    real, dimension(:, :, :), allocatable :: z
    real, dimension(:, :, :), allocatable :: dz8w
    real, dimension(:, :, :), allocatable :: qg_curr

    INTEGER hz, clock1, clock0,delta
    real, dimension(no_runs) :: time
    integer :: CHUNK_X, CHUNK_Y, run
    integer :: NUM_TILES_X, NUM_TILES_Y
    integer :: i, j, k, i_start, i_end, j_start, j_end

    character*255 str

    call chdir("data") ! data directory for both input and output files

    CALL read_value(ids, "ids.bin")
    CALL read_value(ide, "ide.bin")
    CALL read_value(jds, "jds.bin")
    CALL read_value(jde, "jde.bin")
    CALL read_value(kds, "kds.bin")
    CALL read_value(kde, "kde.bin")

    CALL read_value(ims, "ims.bin")
    CALL read_value(ime, "ime.bin")
    CALL read_value(jms, "jms.bin")
    CALL read_value(jme, "jme.bin")
    CALL read_value(kms, "kms.bin")
    CALL read_value(kme, "kme.bin")

    CALL read_value(its, "its.bin")
    CALL read_value(ite, "ite.bin")
    CALL read_value(jts, "jts.bin")
    CALL read_value(jte, "jte.bin")
    CALL read_value(kts, "kts.bin")
    CALL read_value(kte, "kte.bin")

    call getenv( 'NUM_TILES_X', str )
    read(str,*) NUM_TILES_X
    call getenv( 'NUM_TILES_Y', str )
    read(str,*) NUM_TILES_Y

    if(NUM_TILES_X .lt. 1 .or. NUM_TILES_Y .lt. 1) then
      print *, 'Error in environmental variables NUM_TILES_X and NUM_TILES_Y'
      call EXIT
    endif

    CHUNK_X = (ite-its+1)/NUM_TILES_X+1
    CHUNK_Y = (jte-jts+1)/NUM_TILES_Y+1

    CALL read_value_real(dt, "dt.bin")
    CALL read_value(itimestep, "itimestep.bin")
    CALL read_value(hail, "hail.bin")
    CALL read_value(ice2, "ice2.bin")
    CALL read_value_logic(f_qg, "f_qg.bin")

    allocate(ht (ims:ime, jms:jme))
    allocate(rainnc (ims:ime, jms:jme))
    allocate(rainncv (ims:ime, jms:jme))
    allocate(snownc (ims:ime, jms:jme))
    allocate(snowncv (ims:ime, jms:jme))
    allocate(sr (ims:ime, jms:jme))
    allocate(graupelnc (ims:ime, jms:jme))
    allocate(graupelncv (ims:ime, jms:jme))
    CALL read_data_2d(ht, "ht.bin", ims, ime, jms, jme)

    allocate(th       (ims:ime, kms:kme, jms:jme))
    allocate(qv_curr  (ims:ime, kms:kme, jms:jme))
    allocate(qc_curr  (ims:ime, kms:kme, jms:jme))
    allocate(qr_curr  (ims:ime, kms:kme, jms:jme))
    allocate(qi_curr  (ims:ime, kms:kme, jms:jme))
    allocate(qs_curr  (ims:ime, kms:kme, jms:jme))
    allocate(rho      (ims:ime, kms:kme, jms:jme))
    allocate(pi_phy   (ims:ime, kms:kme, jms:jme))
    allocate(p        (ims:ime, kms:kme, jms:jme))
    allocate(z        (ims:ime, kms:kme, jms:jme))
    allocate(dz8w     (ims:ime, kms:kme, jms:jme))
    allocate(qg_curr  (ims:ime, kms:kme, jms:jme))
    CALL read_data(rho     , "rho.bin"     , ims, ime, jms, jme, kms, kme)
    CALL read_data(pi_phy  , "pi_phy.bin"  , ims, ime, jms, jme, kms, kme)
    CALL read_data(p       , "p.bin"       , ims, ime, jms, jme, kms, kme)
    CALL read_data(z       , "z.bin"       , ims, ime, jms, jme, kms, kme)
    CALL read_data(dz8w    , "dz8w.bin"    , ims, ime, jms, jme, kms, kme)

    do run = 0, no_runs

      CALL read_data(th      , "th.bin"      , ims, ime, jms, jme, kms, kme)
      CALL read_data(qv_curr , "qv_curr.bin" , ims, ime, jms, jme, kms, kme)
      CALL read_data(qc_curr , "qc_curr.bin" , ims, ime, jms, jme, kms, kme)
      CALL read_data(qr_curr , "qr_curr.bin" , ims, ime, jms, jme, kms, kme)
      CALL read_data(qi_curr , "qi_curr.bin" , ims, ime, jms, jme, kms, kme)
      CALL read_data(qs_curr , "qs_curr.bin" , ims, ime, jms, jme, kms, kme)
      CALL read_data(qg_curr , "qg_curr.bin" , ims, ime, jms, jme, kms, kme)
      CALL read_data_2d(rainnc, "rainnc.bin", ims, ime, jms, jme)
      CALL read_data_2d(rainncv, "rainncv.bin", ims, ime, jms, jme)
      CALL read_data_2d(snownc, "snownc.bin", ims, ime, jms, jme)
      CALL read_data_2d(snowncv, "snowncv.bin", ims, ime, jms, jme)
      CALL read_data_2d(sr, "sr.bin", ims, ime, jms, jme)
      CALL read_data_2d(graupelnc, "graupelnc.bin", ims, ime, jms, jme)
      CALL read_data_2d(graupelncv, "graupelncv.bin", ims, ime, jms, jme)

      CALL system_clock(count_rate=hz)
      CALL system_clock(count=clock0)

!$OMP PARALLEL DO &
!$OMP PRIVATE ( i, j, i_start, i_end, j_start, j_end) collapse(2)
   DO j = jts, jte, CHUNK_Y
    DO i = its, ite, CHUNK_X
     i_start = i
     i_end = min(i_start+(CHUNK_X-1),ite)
     j_start = j
     j_end = min(j_start+(CHUNK_Y-1),jte)

      CALL gsfcgce(                                                 &
               TH=th                                              &
               ,QV=qv_curr                                        &
               ,QL=qc_curr                                        &
               ,QR=qr_curr                                        &
               ,QI=qi_curr                                        &
               ,QS=qs_curr                                        &
               ,RHO=rho, PII=pi_phy, P=p, DT_IN=dt, Z=z           &
               ,HT=ht, DZ8W=dz8w, GRAV=G                          &
               ,RHOWATER=rhowater, RHOSNOW=rhosnow                &
               ,ITIMESTEP=itimestep                               &
               ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde &
               ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme &
               ,ITS=i_start, ITE=i_end, JTS=j_start, JTE=j_end, KTS=kts, KTE=kte &
!               ,ITS=its,ITE=ite, JTS=jts,JTE=jte, KTS=kts,KTE=kte &
               ,RAINNC=rainnc, RAINNCV=rainncv                    &
               ,SNOWNC=snownc, SNOWNCV=snowncv ,SR=sr             &
               ,GRAUPELNC=graupelnc ,GRAUPELNCV=graupelncv        &
               ,F_QG=f_qg                                         &
               ,QG=qg_curr                                        &
               ,IHAIL=hail, ICE2=ice2                             &
      )

    end do
  end do
!$OMP END PARALLEL DO

      call system_clock(count=clock1)
      delta = clock1 - clock0
      if(run .ge. 1) time(run) = real(delta)/(real(hz))
    
    end do
!    print *, time
    PRINT *, "Original Goddard microphysics computing time [ms]"," ",sum(time)/no_runs*1000.0

    if ( write_output ) THEN
    CALL write_data_2d(rainnc, "rainnc_output.bin", ims, ime, jms, jme)
    CALL write_data_2d(rainncv, "rainncv_output.bin", ims, ime, jms, jme)
    CALL write_data_2d(snownc, "snownc_output.bin", ims, ime, jms, jme)
    CALL write_data_2d(snowncv, "snowncv_output.bin", ims, ime, jms, jme)
    CALL write_data_2d(sr, "sr_output.bin", ims, ime, jms, jme)
    CALL write_data_2d(graupelnc, "graupelnc_output.bin", ims, ime, jms, jme)
    CALL write_data_2d(graupelncv, "graupelncv_output.bin", ims, ime, jms, jme)
    CALL write_data(th, "th_output.bin", ims, ime, jms, jme, kms, kme)
    CALL write_data(qv_curr, "qv_curr_output.bin", ims, ime, jms, jme, kms, kme)
    CALL write_data(qc_curr, "qc_curr_output.bin", ims, ime, jms, jme, kms, kme)
    CALL write_data(qr_curr, "qr_curr_output.bin", ims, ime, jms, jme, kms, kme)
    CALL write_data(qi_curr, "qi_curr_output.bin", ims, ime, jms, jme, kms, kme)
    CALL write_data(qs_curr, "qs_curr_output.bin", ims, ime, jms, jme, kms, kme)
    CALL write_data(qg_curr, "qg_curr_output.bin", ims, ime, jms, jme, kms, kme)
    endif

    CALL compare(th, "th_output.bin", ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
    CALL compare(qv_curr, "qv_curr_output.bin", ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
    CALL compare(qc_curr, "qc_curr_output.bin", ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
    CALL compare(qr_curr, "qr_curr_output.bin", ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
    CALL compare(qi_curr, "qi_curr_output.bin", ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
    CALL compare(qs_curr, "qs_curr_output.bin", ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
    CALL compare(qg_curr, "qg_curr_output.bin", ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
    CALL compare_2d(rainnc, "rainnc_output.bin", ims, ime, jms, jme, its, ite, jts, jte)
    CALL compare_2d(rainncv, "rainncv_output.bin", ims, ime, jms, jme, its, ite, jts, jte)
    CALL compare_2d(snownc, "snownc_output.bin", ims, ime, jms, jme, its, ite, jts, jte)
    CALL compare_2d(snowncv, "snowncv_output.bin", ims, ime, jms, jme, its, ite, jts, jte)
    CALL compare_2d(sr, "sr_output.bin", ims, ime, jms, jme, its, ite, jts, jte)
    CALL compare_2d(graupelnc, "graupelnc_output.bin", ims, ime, jms, jme, its, ite, jts, jte)
    CALL compare_2d(graupelncv, "graupelncv_output.bin", ims, ime, jms, jme, its, ite, jts, jte)

    deallocate(ht)
    deallocate(rainnc)
    deallocate(rainncv)
    deallocate(snownc)
    deallocate(snowncv)
    deallocate(sr)
    deallocate(graupelnc)
    deallocate(graupelncv)
    deallocate(th)
    deallocate(qv_curr)
    deallocate(qc_curr)
    deallocate(qr_curr)
    deallocate(qi_curr)
    deallocate(qs_curr)
    deallocate(rho)
    deallocate(pi_phy)
    deallocate(p)
    deallocate(z)
    deallocate(dz8w)
    deallocate(qg_curr)

  end program gsfcgce_driver

  subroutine read_value(value, file_name)
    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: value
    CHARACTER *(*), INTENT(IN) :: file_name

    open(unit = 666, ACCESS="STREAM", file=file_name, convert="big_endian", action="read")
      read (666) value
    close (666)
  end subroutine read_value


  subroutine read_value_logic(value, file_name)
    IMPLICIT NONE

    LOGICAL, INTENT(OUT) :: value
    CHARACTER *(*), INTENT(IN) :: file_name

    open(unit = 666, ACCESS="STREAM", file=file_name, convert="big_endian", action="read")
      read (666) value
    close (666)
  end subroutine read_value_logic


  subroutine read_value_real(value, file_name)
    IMPLICIT NONE

    real, INTENT(OUT) :: value
    CHARACTER *(*), INTENT(IN) :: file_name

    open(unit = 666, ACCESS="STREAM", file=file_name, convert="big_endian", action="read")
      read (666) value
    close (666)            
  end subroutine read_value_real


  subroutine read_data(data, file_name, ims, ime, jms, jme, kms, kme)
    IMPLICIT NONE

    INTEGER,  INTENT(IN) :: ims, ime, jms, jme, kms, kme
    real, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(OUT) :: data
    CHARACTER *(*), INTENT(IN) :: file_name

    open(unit = 666, ACCESS="STREAM", file=file_name, convert="big_endian", action="read")
    read (666) data
    close (666)
  end subroutine read_data

  subroutine read_data_1d(data, file_name, ims, ime)
    IMPLICIT NONE

    INTEGER,  INTENT(IN) :: ims, ime
    real, DIMENSION( ims:ime ), INTENT(OUT) :: data
    CHARACTER *(*), INTENT(IN) :: file_name

    open(unit = 666, ACCESS="STREAM", file=file_name, convert="big_endian", action="read")
    read (666) data
    close (666)
  end subroutine read_data_1d

  subroutine read_data_2d(data, file_name, ims, ime, jms, jme)
    IMPLICIT NONE

    INTEGER,  INTENT(IN) :: ims, ime, jms, jme
    real, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: data
    CHARACTER *(*), INTENT(IN) :: file_name

    open(unit = 666, ACCESS="STREAM", file=file_name, convert="big_endian", action="read")
    read (666) data
    close (666)
  end subroutine read_data_2d

  subroutine read_data_2d_logical(data, file_name, ims, ime, jms, jme)
    IMPLICIT NONE

    INTEGER,  INTENT(IN) :: ims, ime, jms, jme
    logical, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: data
    CHARACTER *(*), INTENT(IN) :: file_name

    open(unit = 666, ACCESS="STREAM", file=file_name, convert="big_endian", action="read")
    read (666) data
    close (666)
  end subroutine read_data_2d_logical

  subroutine read_data_2d_integer(data, file_name, ims, ime, jms, jme)
    IMPLICIT NONE

    INTEGER,  INTENT(IN) :: ims, ime, jms, jme
    integer, DIMENSION( ims:ime, jms:jme ), INTENT(OUT) :: data
    CHARACTER *(*), INTENT(IN) :: file_name

    open(unit = 666, ACCESS="STREAM", file=file_name, convert="big_endian", action="read")
    read (666) data
    close (666)
  end subroutine read_data_2d_integer


  subroutine write_data(data, file_name, ims, ime, jms, jme, kms, kme)
    IMPLICIT NONE

    INTEGER,  INTENT(IN) :: ims, ime, jms, jme, kms, kme
    real, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN) :: data
    CHARACTER *(*), INTENT(IN) :: file_name

    open(unit = 666, ACCESS="STREAM", file=file_name, convert="big_endian", action="write")
    write (666) data
    close (666)
  end subroutine write_data


  subroutine write_data_2d(data, file_name, ims, ime, jms, jme)
    IMPLICIT NONE

    INTEGER,  INTENT(IN) :: ims, ime, jms, jme
    real, DIMENSION( ims:ime, jms:jme ), INTENT(IN) :: data
    CHARACTER *(*), INTENT(IN) :: file_name

    open(unit = 666, ACCESS="STREAM", file=file_name, convert="big_endian", action="write")
    write (666) data
    close (666)
  end subroutine write_data_2d


 subroutine compare_2d(data, file_name, ims, ime, jms, jme, its, ite, jts, jte)
    IMPLICIT NONE

    INTEGER,  INTENT(IN) :: ims, ime, jms, jme
    INTEGER,  INTENT(IN) :: its, ite, jts, jte
    real, DIMENSION( ims:ime, jms:jme ), INTENT(IN) :: data
    CHARACTER *(*), INTENT(IN) :: file_name

    INTEGER i, j
    real max_abs_err
    real max_rel_err
    real rel_err, abs_err
    real rmse
    integer equal_values, non_equal_values

    real, dimension(:, :), allocatable :: data_orig
    allocate(data_orig    (ims:ime, jms:jme))

    open(unit = 666, ACCESS="STREAM", file=file_name, convert="big_endian", action="read")
    read (666) data_orig
    close (666)

    max_abs_err = 0.0
    max_rel_err = 0.0
    rmse = 0.0
    equal_values = 0
    non_equal_values = 0

    do i = its, ite
      do j = jts, jte
        if (isnan(data_orig(i,j))) stop '"data_orig(i,j)" is a NaN'
        if (isnan(data(i,j))) stop '"data(i,j)" is a NaN'

        if( data(i,j) .eq. data_orig(i,j) ) then
          equal_values = equal_values + 1
        else
          non_equal_values = non_equal_values + 1
        endif

        if( abs(data(i,j)) .ne. 0.0 .and. abs(data_orig(i,j)) .ne. 0.0 ) then
          rel_err = (abs(data(i,j) - data_orig(i,j)) ) / MAX(abs(data(i,j)), abs(data_orig(i,j)) ) ;
        else 
          rel_err = MAX( abs(data(i,j)), abs(data_orig(i,j)) )
        end if
        if(rel_err .gt. max_rel_err) then
          max_rel_err = rel_err
        end if

        abs_err = abs(data(i,j) - data_orig(i,j))
        rmse = rmse + abs_err**2.0
        if(abs_err .gt. max_abs_err) then
          max_abs_err = abs_err
        end if

      end do
    end do

    print *
    print *, file_name
    if(non_equal_values .eq. 0) then
      print *, ' values are bit-exact' 
    else
      print *, 'max relative error: ', max_rel_err
      print *, 'max absolute error: ', max_abs_err
      print *, '# of bit-exact values', equal_values
      print *, '# of non bit-exact values', non_equal_values
      print *, 'rmse: ', rmse
    endif
    deallocate(data_orig)
  end subroutine compare_2d

  subroutine compare(data, file_name, ims, ime, jms, jme, kms, kme, its, ite, jts, jte, kts, kte)
    IMPLICIT NONE

    INTEGER,  INTENT(IN) :: ims, ime, kms, kme, jms, jme
    INTEGER,  INTENT(IN) :: its, ite, kts, kte, jts, jte
    real, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN) :: data
    CHARACTER *(*), INTENT(IN) :: file_name

    INTEGER i, k, j
    real max_abs_err
    real max_rel_err
    real rel_err, abs_err
    real rmse
    integer equal_values, non_equal_values
    integer data_nan, orig_data_nan

    real, dimension(:, :, :), allocatable :: data_orig
    allocate(data_orig    (ims:ime, kms:kme, jms:jme))

    open(unit = 666, ACCESS="STREAM", file=file_name, convert="big_endian", action="read")
    read (666) data_orig
    close (666)

    max_abs_err = 0.0
    max_rel_err = 0.0
    rmse = 0.0
    equal_values = 0
    non_equal_values = 0
    orig_data_nan = 0
    data_nan = 0

    do i = its, ite
      do k = kts, kte
        do j = jts, jte

          if (isnan(data_orig(i,k,j))) orig_data_nan = orig_data_nan + 1
          if (isnan(data(i,k,j))) data_nan = data_nan + 1

          if( data(i,k,j) .eq. data_orig(i,k,j) ) then
            equal_values = equal_values + 1
          else
            non_equal_values = non_equal_values + 1
          endif

          if( abs(data(i,k,j)) .ne. 0.0 .and. abs(data_orig(i,k,j)) .ne. 0.0 ) then
            rel_err = (abs(data(i,k,j) - data_orig(i,k,j)) ) / MAX(abs(data(i,k,j)), abs(data_orig(i,k,j)) ) ;
          else 
            rel_err = MAX( abs(data(i,k,j)), abs(data_orig(i,k,j)) )
          end if
          if(rel_err .gt. max_rel_err) then
            max_rel_err = rel_err
          end if

          abs_err = abs(data(i,k,j) - data_orig(i,k,j))
          rmse = rmse + abs_err**2.0
          if(abs_err .gt. max_abs_err) then
            max_abs_err = abs_err
          end if

        end do
      end do
    end do

    print *
    print *, file_name
    if(non_equal_values .eq. 0) then
      print *, ' values are bit-exact' 
    else
      print *, 'max relative error: ', max_rel_err
      print *, 'max absolute error: ', max_abs_err
      print *, '# of bit-exact values', equal_values
      print *, '# of non bit-exact values', non_equal_values
      print *, 'rmse: ', rmse
      if (orig_data_nan .gt. 0) print *, 'Number of nans in the original data: ', orig_data_nan
      if (data_nan .gt. 0) print *, 'Number of nans in the data: ', data_nan
    endif
    deallocate(data_orig)
  end subroutine compare

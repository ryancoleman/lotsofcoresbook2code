MODULE module_mp_gsfcgce

!      IMPLICIT NONE
      INTEGER, PARAMETER, PRIVATE:: CHUNK = 16

!JJS 1/3/2008     vvvvv

!  common /bt/
   REAL,    PRIVATE ::          rd1,  rd2,   al,   cp

!  common /cont/
   REAL,    PRIVATE ::          c38, c358, c610, c149, &
                               c879, c172, c409,  c76, &
                               c218, c580, c141
!  common /b3cs/
   REAL,    PRIVATE ::           ag,   bg,   as,   bs, &
                                 aw,   bw,  bgh,  bgq, &
                                bsh,  bsq,  bwh,  bwq

!  common /size/
   REAL,    PRIVATE ::          tnw,  tns,  tng,       &
                               roqs, roqg, roqr

!  common /bterv/
   REAL,    PRIVATE ::          zrc,  zgc,  zsc,       &
                                vrc,  vgc,  vsc

!  common /bsnw/
   REAL,    PRIVATE ::          alv,  alf,  als,   t0,   t00,     &
                                avc,  afc,  asc,  rn1,  bnd1,     &
                                rn2, bnd2,  rn3,  rn4,   rn5,     &
                                rn6,  rn7,  rn8,  rn9,  rn10,     &
                              rn101,rn10a, rn11,rn11a,  rn12

   REAL,    PRIVATE ::         rn14, rn15,rn15a, rn16,  rn17,     &
                              rn17a,rn17b,rn17c, rn18, rn18a,     &
                               rn19,rn19a,rn19b, rn20, rn20a,     &
                              rn20b, bnd3, rn21, rn22,  rn23,     &
                              rn23a,rn23b, rn25,rn30a, rn30b,     &
                              rn30c, rn31, beta, rn32

   REAL,    PRIVATE, DIMENSION( 31 ) ::    rn12a, rn12b, rn13, rn25a

!  common /rsnw1/
   REAL,    PRIVATE ::         rn10b, rn10c, rnn191, rnn192,  rn30,     &
                             rnn30a,  rn33,  rn331,  rn332

!
   REAL,    PRIVATE, DIMENSION( 31 )  ::      aa1,  aa2
   DATA aa1/.7939e-7, .7841e-6, .3369e-5, .4336e-5, .5285e-5,     &
           .3728e-5, .1852e-5, .2991e-6, .4248e-6, .7434e-6,     &
           .1812e-5, .4394e-5, .9145e-5, .1725e-4, .3348e-4,     &
           .1725e-4, .9175e-5, .4412e-5, .2252e-5, .9115e-6,     &
           .4876e-6, .3473e-6, .4758e-6, .6306e-6, .8573e-6,     &
           .7868e-6, .7192e-6, .6513e-6, .5956e-6, .5333e-6,     &
           .4834e-6/
   DATA aa2/.4006, .4831, .5320, .5307, .5319,      &
           .5249, .4888, .3894, .4047, .4318,      &
           .4771, .5183, .5463, .5651, .5813,      &
           .5655, .5478, .5203, .4906, .4447,      &
           .4126, .3960, .4149, .4320, .4506,      &
           .4483, .4460, .4433, .4413, .4382,      &
           .4361/

!JJS 1/3/2008     ^^^^^

CONTAINS

!-------------------------------------------------------------------
!  NASA/GSFC GCE
!  Tao et al, 2001, Meteo. & Atmos. Phy., 97-137
!-------------------------------------------------------------------
  SUBROUTINE gsfcgce(  th,                                         &
                       qv, ql,                                     &
                       qr, qi,                                     &
                       qs,                                         &
                       rho, pii, p, dt_in, z,                      &
                       ht, dz8w, grav,                             &
                       rhowater, rhosnow,                          &
                       itimestep,                                  &
                       ids,ide, jds,jde, kds,kde,                  & ! domain dims
                       ims,ime, jms,jme, kms,kme,                  & ! memory dims
                       its,ite, jts,jte, kts,kte,                  & ! tile   dims
                       rainnc, rainncv,                            &
                       snownc, snowncv, sr,                        &
                       graupelnc, graupelncv,                      &
!                       f_qg, qg, pgold,                            &
                       f_qg, qg,                                   &
                       ihail, ice2                                 &
                                                                   )

!-------------------------------------------------------------------
  IMPLICIT NONE
!-------------------------------------------------------------------
!
! JJS 2/15/2005
!
  INTEGER,      INTENT(IN   )    ::   ids,ide, jds,jde, kds,kde , &
                                      ims,ime, jms,jme, kms,kme , &
                                      its,ite, jts,jte, kts,kte 
  INTEGER,      INTENT(IN   )    ::   itimestep, ihail, ice2 

  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        INTENT(INOUT) ::                                          &
                                                              th, &
                                                              qv, &
                                                              ql, &
                                                              qr, &
                                                              qi, &
                                                              qs, &
                                                              qg
!
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        INTENT(IN   ) ::                                          &
                                                             rho, &
                                                             pii, &
                                                               p, &
                                                            dz8w, &
                                                               z

  REAL, DIMENSION( ims:ime , jms:jme ),                           &
        INTENT(INOUT) ::                               rainnc,    &
                                                       rainncv,   &
                                                       snownc,    &   
                                                       snowncv,   &
                                                       sr,        &
                                                       graupelnc, &
                                                       graupelncv 

  REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN) ::       ht

  REAL, INTENT(IN   ) ::                                   dt_in, &
                                                            grav, &
                                                        rhowater, &
                                                         rhosnow

  LOGICAL, INTENT(IN), OPTIONAL :: F_QG

!  LOCAL VAR
  INTEGER ::  itaobraun, istatmin, new_ice_sat, id
  INTEGER :: i, j, k, ip, ii, ic
  INTEGER :: iskip, ih, icount, ibud
  REAL    :: hour
  REAL    , PARAMETER :: cmin=1.e-20
  REAL    :: dth, dqv, dqrest, dqall, dqall1, rhotot, a1, a2 
  LOGICAL :: flag_qg
  REAL, DIMENSION(CHUNK, kts:kte):: th1d, qv1d, ql1d, qr1d
  REAL, DIMENSION(CHUNK, kts:kte):: qi1d, qs1d, qg1d
  REAL, DIMENSION(CHUNK, kts:kte):: rho1d, pii1d, p1d

!
!c  ihail = 0    for graupel, for tropical region
!c  ihail = 1    for hail, for mid-lat region

! itaobraun: 0 for Tao's constantis, 1 for Braun's constants
!c        if ( itaobraun.eq.1 ) --> betah=0.5*beta=-.46*0.5=-0.23;   cn0=1.e-6
!c        if ( itaobraun.eq.0 ) --> betah=0.5*beta=-.6*0.5=-0.30;    cn0=1.e-8
   itaobraun = 1

!c  ice2 = 0    for 3 ice --- ice, snow and graupel/hail
!c  ice2 = 1    for 2 ice --- ice and snow only
!c  ice2 = 2    for 2 ice --- ice and graupel only, use ihail = 0 only
!c  ice2 = 3    for 0 ice --- no ice, warm only

!c  new_ice_sat = 0, 1 or 2 
    new_ice_sat = 2 

!c istatmin
    istatmin = 180

!c id = 0  without in-line staticstics
!c id = 1  with in-line staticstics
    id = 0

!c ibud = 0 no calculation of dth, dqv, dqrest and dqall
!c ibud = 1 yes
    ibud = 0

! calculte fallflux and precipiation in MKS system
   call fall_flux(dt_in, qr, qi, qs, qg, p,                   &
                      rho, z, dz8w, ht, rainnc,               &
                      rainncv, grav,itimestep,                &
                      rhowater, rhosnow,                      &
                      snownc, snowncv, sr,                    &
                      graupelnc, graupelncv,                  &
                      ihail, ice2,                            &
                      ims,ime, jms,jme, kms,kme,              & ! memory dims
                      its,ite, jts,jte, kts,kte               ) ! tile   dims
!-----------------------------------------------------------------------

!c  set up constants used internally in GCE
   call consat_s (ihail, itaobraun)

!c Negative values correction
   iskip = 1
   if (iskip.eq.0) then
      call negcor(qv,rho,dz8w,ims,ime,jms,jme,kms,kme, &
                           itimestep,1,             &
                           its,ite,jts,jte,kts,kte)
      call negcor(ql,rho,dz8w,ims,ime,jms,jme,kms,kme, &
                           itimestep,2,             &
                           its,ite,jts,jte,kts,kte)
      call negcor(qr,rho,dz8w,ims,ime,jms,jme,kms,kme, &
                           itimestep,3,             &
                           its,ite,jts,jte,kts,kte)
      call negcor(qi,rho,dz8w,ims,ime,jms,jme,kms,kme, &
                           itimestep,4,             &
                           its,ite,jts,jte,kts,kte)
      call negcor(qs,rho,dz8w,ims,ime,jms,jme,kms,kme, &
                           itimestep,5,             &
                           its,ite,jts,jte,kts,kte)
      call negcor(qg,rho,dz8w,ims,ime,jms,jme,kms,kme, &
                           itimestep,6,             &
                           its,ite,jts,jte,kts,kte)
   endif ! iskip

!$OMP PARALLEL DO &
!$OMP PRIVATE ( ic, j, ii, i, k ) &
!$OMP PRIVATE ( th1d, qv1d, ql1d, qr1d ) &
!$OMP PRIVATE ( qi1d, qs1d, qg1d ) &
!$OMP PRIVATE ( rho1d, pii1d, p1d ) &
!$OMP SCHEDULE(dynamic,1)
!      ip_loop:
      DO ip = 1,((1+(ite-its+1)/CHUNK)*CHUNK)*(jte-jts+1),CHUNK ! i-dim contains '(1+(ite-its+1)/CHUNK)' blocks of size 'CHUNK'
       ! tid = omp_get_thread_num() + 1   ! not currently used but available and useful for debugging
       j  = jts+(ip-1)/((1+(ite-its+1)/CHUNK)*CHUNK)
       IF ( j .ge. jts .and. j .le. jte ) THEN ! j: [jts, jte]
        ii = its+mod((ip-1),((1+(ite-its+1)/CHUNK)*CHUNK)) ! ii: [its, ((1+(ite-its+1)/CHUNK)*CHUNK)]

         do k = kts, kte
          DO ic=1,min(CHUNK,ite-ii+1) ! ensure max. value for 'i' is 'ite'
            i = ii+ic -1
            th1d(ic,k) = th(i,k,j)
            qv1d(ic,k) = qv(i,k,j)
            ql1d(ic,k) = ql(i,k,j)
            qr1d(ic,k) = qr(i,k,j)
            qi1d(ic,k) = qi(i,k,j)
            qs1d(ic,k) = qs(i,k,j)
            qg1d(ic,k) = qg(i,k,j)
            rho1d(ic,k) = rho(i,k,j)
            pii1d(ic,k) = pii(i,k,j)
            p1d(ic,k) = p(i,k,j)
          ENDDO
         enddo

         IF ( min(CHUNK,ite-ii+1) .gt. 0 ) THEN
!c microphysics in GCE
   call SATICEL_S( dt_in, ihail, itaobraun, ice2, istatmin,     &
                   new_ice_sat, id,                             &
                   th1d, qv1d, ql1d, qr1d,                      &
                   qi1d, qs1d, qg1d,                            &
                   rho1d, pii1d, p1d, itimestep,                & 
                   kts, kte, kms, kme, min(CHUNK,ite-ii+1) )
         ENDIF

         do k = kts, kte
          DO ic=1,min(CHUNK,ite-ii+1)
            i = ii+ic -1
            th(i,k,j) = th1d(ic,k)
            qv(i,k,j) = qv1d(ic,k)
            ql(i,k,j) = ql1d(ic,k)
            qr(i,k,j) = qr1d(ic,k)
            qi(i,k,j) = qi1d(ic,k)
            qs(i,k,j) = qs1d(ic,k)
            qg(i,k,j) = qg1d(ic,k)
          ENDDO
         enddo

         ENDIF
      ENDDO ! ip_loop


   END SUBROUTINE gsfcgce

!-----------------------------------------------------------------------
   SUBROUTINE fall_flux ( dt, qr, qi, qs, qg, p,              &
                      rho, z, dz8w, topo, rainnc,             &
                      rainncv, grav, itimestep,               &
                      rhowater, rhosnow,                      &
                      snownc, snowncv, sr,                    &
                      graupelnc, graupelncv,                  &
                      ihail, ice2,                            &
                      ims,ime, jms,jme, kms,kme,              & ! memory dims
                      its,ite, jts,jte, kts,kte               ) ! tile   dims
!-----------------------------------------------------------------------
! adopted from Jiun-Dar Chern's codes for Purdue Regional Model
! adopted by Jainn J. Shi, 6/10/2005
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER, INTENT(IN   )               :: ihail, ice2,                &
                                          ims,ime, jms,jme, kms,kme,  &
                                          its,ite, jts,jte, kts,kte 
  INTEGER, INTENT(IN   )               :: itimestep
  REAL,    DIMENSION( ims:ime , kms:kme , jms:jme ),                  &
           INTENT(INOUT)               :: qr, qi, qs, qg       
  REAL,    DIMENSION( ims:ime , jms:jme ),                            &
           INTENT(INOUT)               :: rainnc, rainncv,            &
                                          snownc, snowncv, sr,        &
                                          graupelnc, graupelncv
  REAL,    DIMENSION( ims:ime , kms:kme , jms:jme ),                  &
           INTENT(IN   )               :: rho, z, dz8w, p     
  REAL,    INTENT(IN   )               :: dt, grav, rhowater, rhosnow

  REAL,    DIMENSION( ims:ime , jms:jme ),                            &
           INTENT(IN   )               :: topo   

! temperary vars
  REAL,    DIMENSION( kts:kte )           :: sqrhoz
  REAL                                    :: tmp1, term0
  REAL                                :: pptrain, pptsnow,        &
                                         pptgraul, pptice
  REAL                                :: qrz, qiz, qsz, qgz, dzw, prez
  REAL,    DIMENSION( kts:kte )       :: rhoz
   INTEGER                    :: k, i, j
!
  REAL, DIMENSION( kts:kte )    :: vtr, vts, vtg, vti
  REAL                          :: dtb, pi, consta, constc, gambp4,    &
                                   gamdp4, gam4pt5, gam4bbar
!  Lin
   REAL    , PARAMETER ::     xnor = 8.0e6
   REAL    , PARAMETER ::     xnos = 1.6e7   ! Tao's value
   REAL    , PARAMETER ::                              &
             constb = 0.8, constd = 0.11, o6 = 1./6., cdrag = 0.6
!  Lin
  REAL    , PARAMETER ::     xnoh = 2.0e5    ! Tao's value
  REAL    , PARAMETER ::     rhohail = 917.
! Hobbs
  REAL    , PARAMETER ::     xnog = 4.0e6
  REAL    , PARAMETER ::     rhograul = 400.
  REAL    , PARAMETER ::     abar = 19.3, bbar = 0.37, p0 = 1.0e5
  REAL    , PARAMETER ::     rhoe_s = 1.29

! for terminal velocity flux
  INTEGER                       :: min_q, max_q
  REAL                          :: t_del_tv, del_tv, flux, fluxin, fluxout ,tmpqrz
  LOGICAL                       :: notlast

!-----------------------------------------------------------------------
!  This program calculates precipitation fluxes due to terminal velocities.
!-----------------------------------------------------------------------
   dtb=dt
   pi=acos(-1.)
   consta=2115.0*0.01**(1-constb)
   constc=78.63*0.01**(1-constd)

!  Gamma function
   gambp4=ggamma(constb+4.)
   gamdp4=ggamma(constd+4.)
   gam4pt5=ggamma(4.5)
   gam4bbar=ggamma(4.+bbar)

!***********************************************************************
! Calculate precipitation fluxes due to terminal velocities.
!***********************************************************************
!
!- Calculate termianl velocity (vt?)  of precipitation q?z
!- Find maximum vt? to determine the small delta t

!$OMP PARALLEL DO &
!$OMP FIRSTPRIVATE(dtb, pi, consta, constc, gamdp4, gam4pt5, gam4bbar) &
!$OMP PRIVATE(i,k,sqrhoz, tmp1, term0, pptrain, pptsnow) &
!$OMP PRIVATE(pptgraul, pptice, qrz, qiz, qsz, qgz) &
!$OMP PRIVATE(dzw, rhoz) &
!$OMP PRIVATE(vtr, vts, vtg, vti) &
!$OMP PRIVATE(min_q, max_q) &
!$OMP PRIVATE(t_del_tv, del_tv, flux, fluxin, fluxout , tmpqrz,notlast) &
!$OMP SCHEDULE(dynamic)
 j_loop:  do j = jts, jte
 i_loop:  do i = its, ite

   pptrain = 0.
   pptsnow = 0.
   pptgraul = 0.
   pptice  = 0.

   do k = kts, kte
      rhoz(k)=rho(i,k,j)
      sqrhoz(k)=sqrt(rhoe_s/rhoz(k))
   enddo !k

!
!-- rain
!
    t_del_tv=0.
    del_tv=dtb
    notlast=.true.
    DO while (notlast)
!
      min_q=kte
      max_q=kts-1
!
      do k=kts,kte-1
         qrz=qr(i,k,j)
         if (qrz .gt. 1.0e-8) then
            min_q=min0(min_q,k)
            max_q=max0(max_q,k)
            tmp1=sqrt(pi*rhowater*xnor/rhoz(k)/qrz)
            tmp1=sqrt(tmp1)
            vtr(k)=consta*gambp4*sqrhoz(k)/tmp1**constb
            vtr(k)=vtr(k)/6.
            if (k .eq. 1) then
               del_tv=amin1(del_tv,0.9*(z(i,k,j)-topo(i,j))/vtr(k))
            else
               del_tv=amin1(del_tv,0.9*(z(i,k,j)-z(i,k-1,j))/vtr(k))
            endif
         else
            vtr(k)=0.
         endif
      enddo

      if (max_q .ge. min_q) then
!
!- Check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dtb)

         t_del_tv=t_del_tv+del_tv
!
         if ( t_del_tv .ge. dtb ) then
              notlast=.false.
              del_tv=dtb+del_tv-t_del_tv
         endif

! use small delta t to calculate the qrz flux
! termi is the qrz flux pass in the grid box through the upper boundary
! termo is the qrz flux pass out the grid box through the lower boundary
!
         fluxin=0.
         do k=max_q,min_q,-1
            dzw=dz8w(i,k,j)
            qrz=qr(i,k,j)
            fluxout=rhoz(k)*vtr(k)*qrz
            flux=(fluxin-fluxout)/rhoz(k)/dzw
            qrz=qrz+del_tv*flux
            qrz=amax1(0.,qrz)
            qr(i,k,j)=qrz
            fluxin=fluxout
         enddo
         if (min_q .eq. 1) then
            pptrain=pptrain+fluxin*del_tv
         else
             dzw=dz8w(i,min_q-1,j)
             qrz=qr(i,min_q-1,j)
             qrz=qrz+del_tv*  &
                          fluxin/rhoz(min_q-1)/dzw
            qr(i,min_q-1,j)=qrz
         endif
!
      else
         notlast=.false.
      endif
    ENDDO

!
!-- snow
!
    t_del_tv=0.
    del_tv=dtb
    notlast=.true.

    DO while (notlast)
!
      min_q=kte
      max_q=kts-1
!
      do k=kts,kte-1
         qsz=qs(i,k,j)
         if (qsz .gt. 1.0e-8) then
            min_q=min0(min_q,k)
            max_q=max0(max_q,k)
            tmp1=sqrt(pi*rhosnow*xnos/rhoz(k)/qsz)
            tmp1=sqrt(tmp1)
            vts(k)=constc*gamdp4*sqrhoz(k)/tmp1**constd
            vts(k)=vts(k)/6.
            if (k .eq. 1) then
               del_tv=amin1(del_tv,0.9*(z(i,k,j)-topo(i,j))/vts(k))
            else
               del_tv=amin1(del_tv,0.9*(z(i,k,j)-z(i,k-1,j))/vts(k))
            endif
         else
            vts(k)=0.
         endif
      enddo

      if (max_q .ge. min_q) then
!
!
!- Check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dtb)

         t_del_tv=t_del_tv+del_tv

         if ( t_del_tv .ge. dtb ) then
              notlast=.false.
              del_tv=dtb+del_tv-t_del_tv
         endif

! use small delta t to calculate the qsz flux
! termi is the qsz flux pass in the grid box through the upper boundary
! termo is the qsz flux pass out the grid box through the lower boundary
!
         fluxin=0.
         do k=max_q,min_q,-1
            dzw=dz8w(i,k,j)
            qsz=qs(i,k,j)
            fluxout=rhoz(k)*vts(k)*qsz
            flux=(fluxin-fluxout)/rhoz(k)/dzw
            qsz=qsz+del_tv*flux
            qsz=amax1(0.,qsz)
            qs(i,k,j)=qsz
            fluxin=fluxout
         enddo
         if (min_q .eq. 1) then
            pptsnow=pptsnow+fluxin*del_tv
         else
            dzw=dz8w(i,min_q-1,j)
            qsz=qs(i,min_q-1,j)
            qsz=qsz+del_tv*  &
                         fluxin/rhoz(min_q-1)/dzw
            qs(i,min_q-1,j)=qsz
         endif
!
      else
         notlast=.false.
      endif

    ENDDO

!
!   ice2=0 --- with hail/graupel 
!   ice2=1 --- without hail/graupel 
!
  if (ice2.eq.0) then 
!
!-- If IHAIL=1, use hail.
!-- If IHAIL=0, use graupel.
!
!    if (ihail .eq. 1) then
!       xnog = xnoh
!       rhograul = rhohail
!    endif

    t_del_tv=0.
    del_tv=dtb
    notlast=.true.
!
    DO while (notlast)
!
      min_q=kte
      max_q=kts-1
!
      do k=kts,kte-1
        IF (ice2 .eq. 0) THEN
          qgz=qg(i,k,j)
        ELSE
          qgz=0.
        ENDIF
         if (qgz .gt. 1.0e-8) then
            if (ihail .eq. 1) then
!  for hail, based on Lin et al (1983)
              min_q=min0(min_q,k)
              max_q=max0(max_q,k)
              tmp1=sqrt(pi*rhohail*xnoh/rhoz(k)/qgz)
              tmp1=sqrt(tmp1)
              term0=sqrt(4.*grav*rhohail/3./rhoz(k)/cdrag)
              vtg(k)=gam4pt5*term0*sqrt(1./tmp1)
              vtg(k)=vtg(k)/6.
              if (k .eq. 1) then
                 del_tv=amin1(del_tv,0.9*(z(i,k,j)-topo(i,j))/vtg(k))
              else
                 del_tv=amin1(del_tv,0.9*(z(i,k,j)-z(i,k-1,j))/vtg(k))
              endif !k
            else
! added by JJS
! for graupel, based on RH (1984)
              min_q=min0(min_q,k)
              max_q=max0(max_q,k)
              tmp1=sqrt(pi*rhograul*xnog/rhoz(k)/qgz)
              tmp1=sqrt(tmp1)
              tmp1=tmp1**bbar
              tmp1=1./tmp1
              term0=abar*gam4bbar/6.
              prez = p(i,k,j)
              vtg(k)=term0*tmp1*(p0/prez)**0.4
              if (k .eq. 1) then
                 del_tv=amin1(del_tv,0.9*(z(i,k,j)-topo(i,j))/vtg(k))
              else
                 del_tv=amin1(del_tv,0.9*(z(i,k,j)-z(i,k-1,j))/vtg(k))
              endif !k
            endif !ihail
         else
            vtg(k)=0.
         endif !qgz
      enddo !k

      if (max_q .ge. min_q) then
!
!
!- Check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dtb)

         t_del_tv=t_del_tv+del_tv

         if ( t_del_tv .ge. dtb ) then
              notlast=.false.
              del_tv=dtb+del_tv-t_del_tv
         endif

! use small delta t to calculate the qgz flux
! termi is the qgz flux pass in the grid box through the upper boundary
! termo is the qgz flux pass out the grid box through the lower boundary
!
         fluxin=0.
         do k=max_q,min_q,-1
            IF (ice2 .eq. 0) THEN
              qgz=qg(i,k,j)
            ELSE
              qgz=0.
            ENDIF
            dzw=dz8w(i,k,j)
            fluxout=rhoz(k)*vtg(k)*qgz
            flux=(fluxin-fluxout)/rhoz(k)/dzw
            qgz=qgz+del_tv*flux
            qgz=amax1(0.,qgz)
            qg(i,k,j)=qgz
            fluxin=fluxout
         enddo
         if (min_q .eq. 1) then
            pptgraul=pptgraul+fluxin*del_tv
         else
            IF (ice2 .eq. 0) THEN
              qgz=qg(i,min_q-1,j)
            ELSE
              qgz=0.
            ENDIF
            dzw=dz8w(i,min_q-1,j)
            qgz=qgz+del_tv*  &
                         fluxin/rhoz(min_q-1)/dzw
            qg(i,min_q-1,j)=qgz
         endif
!
      else
         notlast=.false.
      endif
!
    ENDDO
 ENDIF !ice2
!
!-- cloud ice  (03/21/02) follow Vaughan T.J. Phillips at GFDL
!

    t_del_tv=0.
    del_tv=dtb
    notlast=.true.
!
    DO while (notlast)
!
      min_q=kte
      max_q=kts-1
!
      do k=kts,kte-1
         qiz=qi(i,k,j)
         if (qiz .gt. 1.0e-8) then
            min_q=min0(min_q,k)
            max_q=max0(max_q,k)
            vti(k)= 3.29 * (rhoz(k)* qiz)** 0.16  ! Heymsfield and Donner
            if (k .eq. 1) then
               del_tv=amin1(del_tv,0.9*(z(i,k,j)-topo(i,j))/vti(k))
            else
               del_tv=amin1(del_tv,0.9*(z(i,k,j)-z(i,k-1,j))/vti(k))
            endif
         else
            vti(k)=0.
         endif
      enddo

      if (max_q .ge. min_q) then
!
!
!- Check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dtb)

         t_del_tv=t_del_tv+del_tv

         if ( t_del_tv .ge. dtb ) then
              notlast=.false.
              del_tv=dtb+del_tv-t_del_tv
         endif

! use small delta t to calculate the qiz flux
! termi is the qiz flux pass in the grid box through the upper boundary
! termo is the qiz flux pass out the grid box through the lower boundary
!

         fluxin=0.
         do k=max_q,min_q,-1
            qiz=qi(i,k,j)
            dzw=dz8w(i,k,j)
            fluxout=rhoz(k)*vti(k)*qiz
            flux=(fluxin-fluxout)/rhoz(k)/dzw
            qiz=qiz+del_tv*flux
            qiz=amax1(0.,qiz)
            qi(i,k,j)=qiz
            fluxin=fluxout
         enddo
         if (min_q .eq. 1) then
            pptice=pptice+fluxin*del_tv
         else
            qiz=qi(i,min_q-1,j)
            dzw=dz8w(i,min_q-1,j)
            qiz=qiz+del_tv*  &
                         fluxin/rhoz(min_q-1)/dzw
            qi(i,min_q-1,j)=qiz
         endif
!
      else
         notlast=.false.
      endif
!
   ENDDO !notlast
   snowncv(i,j) = pptsnow
   snownc(i,j) = snownc(i,j) + pptsnow
   graupelncv(i,j) = pptgraul
   graupelnc(i,j) = graupelnc(i,j) + pptgraul 
   RAINNCV(i,j) = pptrain + pptsnow + pptgraul + pptice                 
   RAINNC(i,j)  = RAINNC(i,j) + pptrain + pptsnow + pptgraul + pptice
   sr(i,j) = 0.
   if (RAINNCV(i,j) .gt. 0.) sr(i,j) = (pptsnow + pptgraul + pptice) / RAINNCV(i,j) 

  ENDDO i_loop
  ENDDO j_loop

  END SUBROUTINE fall_flux

!----------------------------------------------------------------
   REAL FUNCTION ggamma(X)
!----------------------------------------------------------------
   IMPLICIT NONE
!----------------------------------------------------------------
      REAL, INTENT(IN   ) :: x
      REAL, DIMENSION(8)  :: B
      INTEGER             ::j, K1
      REAL                ::PF, G1TO2 ,TEMP

      DATA B/-.577191652,.988205891,-.897056937,.918206857,  &
             -.756704078,.482199394,-.193527818,.035868343/

      PF=1.
      TEMP=X
      DO 10 J=1,200
      IF (TEMP .LE. 2) GO TO 20
      TEMP=TEMP-1.
   10 PF=PF*TEMP
  100 FORMAT(//,5X,'module_gsfcgce: INPUT TO GAMMA FUNCTION TOO LARGE, X=',E12.5)
   20 G1TO2=1.
      TEMP=TEMP - 1.
      DO 30 K1=1,8
   30 G1TO2=G1TO2 + B(K1)*TEMP**K1
      ggamma=PF*G1TO2

      END FUNCTION ggamma

!-----------------------------------------------------------------------
!c Correction of negative values  
   SUBROUTINE negcor ( X, rho, dz8w,                         &
                      ims,ime, jms,jme, kms,kme,              & ! memory dims
                      itimestep, ics,                         &
                      its,ite, jts,jte, kts,kte               ) ! tile   dims
!-----------------------------------------------------------------------
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        INTENT(INOUT) ::                                     X   
  REAL, DIMENSION( ims:ime , kms:kme , jms:jme ),                 &
        INTENT(IN   ) ::                              rho, dz8w  
  integer, INTENT(IN   ) ::                           itimestep, ics 

!c Local variables
  REAL   ::   A0, A1, A2

  A1=0.
  A2=0.

!$OMP PARALLEL DO &
!$OMP SCHEDULE(dynamic)
  do j=jts,jte
    do k=kts,kte
        do i=its,ite
        A1=A1+max(X(i,k,j), 0.)*rho(i,k,j)*dz8w(i,k,j)
        A2=A2+max(-X(i,k,j), 0.)*rho(i,k,j)*dz8w(i,k,j)
        enddo
     enddo
  enddo

  A0=0.0

  if (A1.NE.0.0.and.A1.GT.A2) then 
     A0=(A1-A2)/A1

!$OMP PARALLEL DO &
!$OMP SCHEDULE(dynamic)
   do j=jts,jte
     do k=kts,kte
           do i=its,ite
           X(i,k,j)=A0*AMAX1(X(i,k,j), 0.0)
           enddo
        enddo
     enddo
  endif

  END SUBROUTINE negcor

 SUBROUTINE consat_s (ihail,itaobraun)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
!   Tao, W.-K., and J. Simpson, 1989: Modeling study of a tropical     c
!   squall-type convective line. J. Atmos. Sci., 46, 177-202.          c
!                                                                      c
!   Tao, W.-K., J. Simpson and M. McCumber, 1989: An ice-water         c
!   saturation adjustment. Mon. Wea. Rev., 117, 231-235.               c
!                                                                      c
!   Tao, W.-K., and J. Simpson, 1993: The Goddard Cumulus Ensemble     c
!   Model. Part I: Model description. Terrestrial, Atmospheric and     c
!   Oceanic Sciences, 4, 35-72.                                        c
!                                                                      c
!   Tao, W.-K., J. Simpson, D. Baker, S. Braun, M.-D. Chou, B.         c
!   Ferrier,D. Johnson, A. Khain, S. Lang,  B. Lynn, C.-L. Shie,       c
!   D. Starr, C.-H. Sui, Y. Wang and P. Wetzel, 2003: Microphysics,    c
!   radiation and surface processes in the Goddard Cumulus Ensemble    c
!   (GCE) model, A Special Issue on Non-hydrostatic Mesoscale          c
!   Modeling, Meteorology and Atmospheric Physics, 82, 97-137.         c
!                                                                      c
!   Lang, S., W.-K. Tao, R. Cifelli, W. Olson, J. Halverson, S.        c
!   Rutledge, and J. Simpson, 2007: Improving simulations of           c
!   convective system from TRMM LBA: Easterly and Westerly regimes.    c
!   J. Atmos. Sci., 64, 1141-1164.                                     c
!                                                                      c
!   Coded by Tao (1989-2003), modified by S. Lang (2006/07)            c
!                                                                      c
!   Implemented into WRF  by Roger Shi 2006/2007                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!        itaobraun=0   ! see Tao and Simpson (1993)
!        itaobraun=1   ! see Tao et al. (2003)

 integer :: itaobraun
 real    :: cn0

!JJS 1/3/2008  vvvvv
!JJS   the following common blocks have been moved to the top of
!JJS   module_mp_gsfcgce_driver_instat.F
!
! real,   dimension (1:31) ::  a1, a2
! data a1/.7939e-7,.7841e-6,.3369e-5,.4336e-5,.5285e-5,.3728e-5, &
!         .1852e-5,.2991e-6,.4248e-6,.7434e-6,.1812e-5,.4394e-5,.9145e-5, &
!         .1725e-4,.3348e-4,.1725e-4,.9175e-5,.4412e-5,.2252e-5,.9115e-6, &
!         .4876e-6,.3473e-6,.4758e-6,.6306e-6,.8573e-6,.7868e-6,.7192e-6, &
!         .6513e-6,.5956e-6,.5333e-6,.4834e-6/
! data a2/.4006,.4831,.5320,.5307,.5319,.5249,.4888,.3894,.4047, &
!         .4318,.4771,.5183,.5463,.5651,.5813,.5655,.5478,.5203,.4906, &
!         .4447,.4126,.3960,.4149,.4320,.4506,.4483,.4460,.4433,.4413, &
!         .4382,.4361/
!JJS 1/3/2008  ^^^^^


!     ******************************************************************
!JJS
      al = 2.5e10
      cp = 1.004e7
      rd1 = 1.e-3
      rd2 = 2.2
!JJS
      cpi=4.*atan(1.)
      cpi2=cpi*cpi
      grvt=980.
      cd1=6.e-1
      cd2=4.*grvt/(3.*cd1)
      tca=2.43e3
      dwv=.226
      dva=1.718e-4
      amw=18.016
      ars=8.314e7
      scv=2.2904487
      t0=273.16
      t00=238.16
      alv=2.5e10
      alf=3.336e9
      als=2.8336e10
      avc=alv/cp
      afc=alf/cp
      asc=als/cp
      rw=4.615e6
      cw=4.187e7
      ci=2.093e7
      c76=7.66
      c358=35.86
      c172=17.26939
      c409=4098.026
      c218=21.87456
      c580=5807.695
      c610=6.1078e3
      c149=1.496286e-5
      c879=8.794142
      c141=1.4144354e7
!***   DEFINE THE COEFFICIENTS USED IN TERMINAL VELOCITY
!***   DEFINE THE DENSITY AND SIZE DISTRIBUTION OF PRECIPITATION
!**********   HAIL OR GRAUPEL PARAMETERS   **********
      if(ihail .eq. 1) then
         roqg=.9
         ag=sqrt(cd2*roqg)
         bg=.5
         tng=.002
      else
         roqg=.4
         ag=351.2
!        AG=372.3 ! if ice913=1 6/15/02 tao's
         bg=.37
         tng=.04
      endif
!**********         SNOW PARAMETERS        **********
!ccshie 6/15/02 tao's
!      TNS=1.
!      TNS=.08 ! if ice913=1, tao's
      tns=.16 ! if ice913=0, tao's
      roqs=.1
      as=78.63
      bs=.11
!**********         RAIN PARAMETERS        **********
      aw=2115.
      bw=.8
      roqr=1.
      tnw=.08
!*****************************************************************
      bgh=.5*bg
      bsh=.5*bs
      bwh=.5*bw
      bgq=.25*bg
      bsq=.25*bs
      bwq=.25*bw
!**********GAMMA FUNCTION CALCULATIONS*************
      ga3b  = gammagce(3.+bw)
      ga4b  = gammagce(4.+bw)
      ga6b  = gammagce(6.+bw)
      ga5bh = gammagce((5.+bw)/2.)
      ga3g  = gammagce(3.+bg)
      ga4g  = gammagce(4.+bg)
      ga5gh = gammagce((5.+bg)/2.)
      ga3d  = gammagce(3.+bs)
      ga4d  = gammagce(4.+bs)
      ga5dh = gammagce((5.+bs)/2.)
!CCCCC        LIN ET AL., 1983 OR LORD ET AL., 1984   CCCCCCCCCCCCCCCCC
      ac1=aw
!JJS
      ac2=ag
      ac3=as
!JJS
      bc1=bw
      cc1=as
      dc1=bs
      zrc=(cpi*roqr*tnw)**0.25
      zsc=(cpi*roqs*tns)**0.25
      zgc=(cpi*roqg*tng)**0.25
      vrc=aw*ga4b/(6.*zrc**bw)
      vsc=as*ga4d/(6.*zsc**bs)
      vgc=ag*ga4g/(6.*zgc**bg)
!     ****************************
      rn1=9.4e-15 ! 6/15/02 tao's
      bnd1=6.e-4
      rn2=1.e-3
      bnd2=2.0e-3 ! if ice913=0 6/15/02 tao's
      rn3=.25*cpi*tns*cc1*ga3d
      esw=1.
      rn4=.25*cpi*esw*tns*cc1*ga3d
      eri=.1  ! 6/17/02 tao's ice913=0 (not 1)
      rn5=.25*cpi*eri*tnw*ac1*ga3b
      ami=1./(24.*6.e-9) ! 6/15/02 tao's
      rn6=cpi2*eri*tnw*ac1*roqr*ga6b*ami
      esr=.5 ! 6/15/02 for ice913=0 tao's
      rn7=cpi2*esr*tnw*tns*roqs
      esr=1.
      rn8=cpi2*esr*tnw*tns*roqr
      rn9=cpi2*tns*tng*roqs
      rn10=2.*cpi*tns
      rn101=.31*ga5dh*sqrt(cc1)
      rn10a=als*als/rw
!JJS
       rn10b=alv/tca
       rn10c=ars/(dwv*amw)
!JJS
      rn11=2.*cpi*tns/alf
      rn11a=cw/alf
      ami50=3.84e-6 ! 6/15/02 tao's
      ami40=3.08e-8 ! 6/15/02 tao's
      eiw=1.
      ui50=100. ! 6/15/02 tao's
      ri50=2.*5.e-3
      cmn=1.05e-15
      rn12=cpi*eiw*ui50*ri50**2

      do 10 k=1,31
         y1=1.-aa2(k)
         rn13(k)=aa1(k)*y1/(ami50**y1-ami40**y1)
         rn12a(k)=rn13(k)/ami50
         rn12b(k)=aa1(k)*ami50**aa2(k)
         rn25a(k)=aa1(k)*cmn**aa2(k)
   10 continue

      egw=1.
      rn14=.25*cpi*egw*tng*ga3g*ag
      egi=.1
      rn15=.25*cpi*egi*tng*ga3g*ag
      egi=1.
      rn15a=.25*cpi*egi*tng*ga3g*ag
      egr=1.
      rn16=cpi2*egr*tng*tnw*roqr
      rn17=2.*cpi*tng
      rn17a=.31*ga5gh*sqrt(ag)
      rn17b=cw-ci
      rn17c=cw
      apri=.66
      bpri=1.e-4
      bpri=0.5*bpri ! 6/17/02 tao's
      rn18=20.*cpi2*bpri*tnw*roqr
      rn18a=apri
      rn19=2.*cpi*tng/alf
      rn19a=.31*ga5gh*sqrt(ag)
      rn19b=cw/alf
!
       rnn191=.78
       rnn192=.31*ga5gh*sqrt(ac2/dva)
!
      rn20=2.*cpi*tng
      rn20a=als*als/rw
      rn20b=.31*ga5gh*sqrt(ag)
      bnd3=2.e-3
      rn21=1.e3*1.569e-12/0.15
      erw=1.
      rn22=.25*cpi*erw*ac1*tnw*ga3b
      rn23=2.*cpi*tnw
      rn23a=.31*ga5bh*sqrt(ac1)
      rn23b=alv*alv/rw


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc
!cc        "c0" in routine      "consat" (2d), "consatrh" (3d)
!cc        if ( itaobraun.eq.1 ) --> betah=0.5*beta=-.46*0.5=-0.23;   cn0=1.e-6
!cc        if ( itaobraun.eq.0 ) --> betah=0.5*beta=-.6*0.5=-0.30;    cn0=1.e-8

       if (itaobraun .eq. 0) then
         cn0=1.e-8
         beta=-.6
       elseif (itaobraun .eq. 1) then
         cn0=1.e-6
         beta=-.46
       endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      rn25=cn0
      rn30a=alv*als*amw/(tca*ars)
      rn30b=alv/tca
      rn30c=ars/(dwv*amw)
      rn31=1.e-17

      rn32=4.*51.545e-4
!
      rn30=2.*cpi*tng
      rnn30a=alv*alv*amw/(tca*ars)
!
      rn33=4.*tns
       rn331=.65
       rn332=.44*sqrt(ac3/dva)*ga5dh
!

    return
 END SUBROUTINE consat_s

 SUBROUTINE saticel_s (dt, ihail, itaobraun, ice2, istatmin, &
                       new_ice_sat, id, &
                       ptwrf, qvwrf, qlwrf, qrwrf, &
                       qiwrf, qswrf, qgwrf, &
                       rho_mks, pi_mks, p0_mks,itimestep, &
                       kts, kte, kms, kme, irestrict )
    IMPLICIT NONE
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
!   Tao, W.-K., and J. Simpson, 1989: Modeling study of a tropical     c
!   squall-type convective line. J. Atmos. Sci., 46, 177-202.          c
!                                                                      c
!   Tao, W.-K., J. Simpson and M. McCumber, 1989: An ice-water         c
!   saturation adjustment. Mon. Wea. Rev., 117, 231-235.               c
!                                                                      c
!   Tao, W.-K., and J. Simpson, 1993: The Goddard Cumulus Ensemble     c
!   Model. Part I: Model description. Terrestrial, Atmospheric and     c
!   Oceanic Sciences, 4, 35-72.                                        c
!                                                                      c
!   Tao, W.-K., J. Simpson, D. Baker, S. Braun, M.-D. Chou, B.         c
!   Ferrier,D. Johnson, A. Khain, S. Lang,  B. Lynn, C.-L. Shie,       c
!   D. Starr, C.-H. Sui, Y. Wang and P. Wetzel, 2003: Microphysics,    c
!   radiation and surface processes in the Goddard Cumulus Ensemble    c
!   (GCE) model, A Special Issue on Non-hydrostatic Mesoscale          c
!   Modeling, Meteorology and Atmospheric Physics, 82, 97-137.         c
!                                                                      c
!   Lang, S., W.-K. Tao, R. Cifelli, W. Olson, J. Halverson, S.        c
!   Rutledge, and J. Simpson, 2007: Improving simulations of           c
!   convective system from TRMM LBA: Easterly and Westerly regimes.    c
!   J. Atmos. Sci., 64, 1141-1164.                                     c
!                                                                      c
!   Tao, W.-K., J. J. Shi,  S. Lang, C. Peters-Lidard, A. Hou, S.      c
!   Braun, and J. Simpson, 2007: New, improved bulk-microphysical      c
!   schemes for studying precipitation processes in WRF. Part I:       c
!   Comparisons with other schemes. to appear on Mon. Wea. Rev.        C
!                                                                      c
!   Coded by Tao (1989-2003), modified by S. Lang (2006/07)            c
!                                                                      c
!   Implemented into WRF  by Roger Shi 2006/2007                       c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!      COMPUTE ICE PHASE MICROPHYSICS AND SATURATION PROCESSES
!
  integer,    parameter ::  nt=2880, nt2=2*nt 

  integer, INTENT(IN)  ::   ice2,ihail,new_ice_sat,id,istatmin
  integer, INTENT(INOUT)  ::   itaobraun
  integer, INTENT(IN)  ::   itimestep
  real, dimension (CHUNK) ::  tairccri
  integer, INTENT(IN) :: kts, kte, kms, kme, irestrict

  real     ::   cn0, dt
  real, INTENT(INOUT), dimension (CHUNK, kms:kme ) ::  ptwrf, qvwrf 
  real, INTENT(INOUT), dimension (CHUNK, kms:kme ) ::  qlwrf, qrwrf,        &
                                                  qiwrf, qswrf, qgwrf
  real, INTENT(IN), dimension (CHUNK, kms:kme) ::  rho_mks
  real, INTENT(IN), dimension (CHUNK, kms:kme) ::  pi_mks
  real, INTENT(IN), dimension (CHUNK, kms:kme) ::  p0_mks

  integer i,j,k, kp

  real, dimension (CHUNK) :: rp0, pi0, pir, pr0, r00, rr0, rrs, rrq, fv0
  real, dimension (CHUNK) :: fvs, zrr, zsr, zgr, cp409, cv409, cp580, cs580
  real, dimension (CHUNK) :: alvr, afcp, avcp, ascp, vrcf, vscf, vgcf, vgcr
  real, dimension (CHUNK) :: dwvp, r3f ,r4f ,r5f ,r6f ,r7r ,r8r ,r9r, r101f
  real, dimension (CHUNK) :: r10ar, r11rt, r12r, r14r, r14f, r15r, r15ar
  real, dimension (CHUNK) :: r15af, r16r, r17r, r17aq, r17as, r18r, r19rt
  real, dimension (CHUNK) :: r19aq, r19as, r20bq, r20bs, r22f, r23af, r23br
  real, dimension (CHUNK) :: r25rt, r31r, r32rt
  real, dimension (CHUNK) :: a1, a2, del, ee1, ee2, temp, r_nci
  real :: ami100 ,ami40 ,ami50,betah &
   ,bg3 ,bgh5 ,bs3 ,bs6 ,bsh5 ,bw3 ,bw6 ,bwh5 ,cmin ,cmin1 ,cmin2  &
   ,d2t , f2 ,f3 ,ft &
   ,qb0,r10t ,r11at &
   ,r19bt ,r20t ,r23t ,r25a ,r2ice &
   ,rt0 ,tb0

  real, dimension (CHUNK,kts:kte) ::  fv
  real, dimension (CHUNK,kts:kte) ::  dpt, dqv
  real, dimension (CHUNK,kts:kte) ::  qcl, qrn,      &
                                                qci, qcs, qcg
  real, dimension (CHUNK) ::        &
           vg,      zg,       &
           ps,      pg,       &
          prn,     psn,       &
        pwacs,   wgacr,       &
        pidep,    pint,       &
          qsi,     ssi,       &
          esi,     esw,       &
          qsw,      pr,       &
          ssw,   pihom,       &
         pidw,   pimlt,       &
        psaut,   qracs,       &
        psaci,   psacw,       &
        qsacw,   praci,       &
        pmlts,   pmltg,       &
        asss,       y1,    y2
  real, dimension (CHUNK) ::        &
        praut,   pracw,       &
         psfw,    psfi,       &
        dgacs,   dgacw,       &
        dgaci,   dgacr,       &
        pgacs,   wgacs,       &
        qgacw,   wgaci,       &
        qgacr,   pgwet,       &
        pgaut,   pracs,       &
        psacr,   qsacr,       &
         pgfr,   psmlt,       &
        pgmlt,   psdep,       &
        pgdep,   piacr,       &
           y5,     scv,       &
          tca,     dwv,       &
          egs,      y3,       &
           y4,     ddb
  real, dimension (CHUNK) ::        &
           pt,      qv,       &
           qc,      qr,       &
           qi,      qs,       &
           qg,    tair,       &
        tairc,   rtair,       &
          dep,      dd,       &
          dd1,     qvs,       &
           dm,      rq,       &
        rsub1,     col,       &
          cnd,     ern,       &
         dlt1,    dlt2,       &
         dlt3,    dlt4,       &
           zr,      vr,       &
           zs,      vs,       &
                 pssub,       &
        pgsub,     dda
  real, dimension (CHUNK,kts:kte) ::  rho
  real, dimension (CHUNK,kts:kte) ::  p0, pi, f0
  real, dimension (nt) ::    & 
          tqc,     tqr,    tqi,    tqs,    tqg
  real, dimension (CHUNK) ::     &
           y0,     ts0,   qss0
  integer, dimension (CHUNK) ::        it  
  integer, dimension (CHUNK, 4) ::    ics 
  integer :: iwarm
  real :: r2is, r2ig

!JJS  convert from mks to cgs, and move from WRF grid to GCE grid
      do k=kts,kte
!dir$ vector aligned
         DO i=1,irestrict
         rho(i,k)=rho_mks(i,k)*0.001
         p0(i,k)=p0_mks(i,k)*10.0
         pi(i,k)=pi_mks(i,k)
         dpt(i,k)=ptwrf(i,k)
         dqv(i,k)=qvwrf(i,k)
         qcl(i,k)=qlwrf(i,k)
         qrn(i,k)=qrwrf(i,k)
         qci(i,k)=qiwrf(i,k)
         qcs(i,k)=qswrf(i,k)
         qcg(i,k)=qgwrf(i,k)
         enddo !i
      enddo !k

      do k=kts,kte
!dir$ vector aligned
         DO i=1,irestrict
         fv(i,k)=sqrt(rho(i,2)/rho(i,k))
         enddo !i
      enddo !k
!JJS

!     ******   THREE CLASSES OF ICE-PHASE   (LIN ET AL, 1983)  *********
         d2t=dt
!       itaobraun=0 ! original pint and pidep & see Tao and Simpson 1993
        itaobraun=1 ! see Tao et al. (2003)
!
       if ( itaobraun.eq.0 ) then
          cn0=1.e-8
       elseif ( itaobraun.eq.1 ) then
          cn0=1.e-6
       endif
         r2ig=1.
         r2is=1.
          if (ice2 .eq. 1) then
              r2ig=0.
              r2is=1.
          endif
          if (ice2 .eq. 2) then
              r2ig=1.
              r2is=0.
          endif
!C  TAO 2007 END
     
!JJS  10/7/2008
!   ICE2=3 ! no ice, warm rain only
    iwarm = 0
    if (ice2 .eq. 3 ) iwarm = 1
      cmin=1.e-19
      cmin1=1.e-20
      cmin2=1.e-12

!JJScap $doacross local(j,i)
!dir$ vector aligned
       DO i=1,irestrict
          it(i)=1
       enddo

      f2=rd1*d2t
      f3=rd2*d2t

      ft=dt/d2t
      rt0=1./(t0-t00)
      bw3=bw+3.
      bs3=bs+3.
      bg3=bg+3.
      bsh5=2.5+bsh
      bgh5=2.5+bgh
      bwh5=2.5+bwh
      bw6=bw+6.
      bs6=bs+6.
      betah=.5*beta
      r10t=rn10*d2t
      r11at=rn11a*d2t
      r19bt=rn19b*d2t
      r20t=-rn20*d2t
      r23t=-rn23*d2t
      r25a=rn25

!     ami50 for use in PINT
       ami50=3.76e-8
       ami100=1.51e-7
       ami40=2.41e-8

!C    ******************************************************************


      do k=kts,kte
         kp=k+1
         tb0=0.
         qb0=0.
!dir$ vector aligned
         DO i=1,irestrict
         rp0(i)=3.799052e3/p0(i,k)
         pi0(i)=pi(i,k)
         pir(i)=1./(pi(i,k))
         pr0(i)=1./p0(i,k)
         r00(i)=rho(i,k)
         rr0(i)=1./rho(i,k)
         rrs(i)=sqrt(rr0(i))
         rrq(i)=sqrt(rrs(i))
         f0(i,k)=al/cp/pi(i,k)
         fv0(i)=fv(i,k)
         fvs(i)=sqrt(fv(i,k))
         zrr(i)=1.e5*zrc*rrq(i)
         zsr(i)=1.e5*zsc*rrq(i)
         zgr(i)=1.e5*zgc*rrq(i)
         cp409(i)=c409*pi0(i)
         cv409(i)=c409*avc
         cp580(i)=c580*pi0(i)
         cs580(i)=c580*asc
         alvr(i)=r00(i)*alv
         afcp(i)=afc*pir(i)
         avcp(i)=avc*pir(i)
         ascp(i)=asc*pir(i)
         vrcf(i)=vrc*fv0(i)
         vscf(i)=vsc*fv0(i)
         vgcf(i)=vgc*fv0(i)
         vgcr(i)=vgc*rrs(i)
         dwvp(i)=c879*pr0(i)
         r3f(i)=rn3*fv0(i)
         r4f(i)=rn4*fv0(i)
         r5f(i)=rn5*fv0(i)
         r6f(i)=rn6*fv0(i)
         r7r(i)=rn7*rr0(i)
         r8r(i)=rn8*rr0(i)
         r9r(i)=rn9*rr0(i)
         r101f(i)=rn101*fvs(i)
         r10ar(i)=rn10a*r00(i)
         r11rt(i)=rn11*rr0(i)*d2t
         r12r(i)=rn12*r00(i)
         r14r(i)=rn14*rrs(i)
         r14f(i)=rn14*fv0(i)
         r15r(i)=rn15*rrs(i)
         r15ar(i)=rn15a*rrs(i)
         r15af(i)=rn15a*fv0(i)
         r16r(i)=rn16*rr0(i)
         r17r(i)=rn17*rr0(i)
         r17aq(i)=rn17a*rrq(i)
         r17as(i)=rn17a*fvs(i)
         r18r(i)=rn18*rr0(i)
         r19rt(i)=rn19*rr0(i)*d2t
         r19aq(i)=rn19a*rrq(i)
         r19as(i)=rn19a*fvs(i)
         r20bq(i)=rn20b*rrq(i)
         r20bs(i)=rn20b*fvs(i)
         r22f(i)=rn22*fv0(i)
         r23af(i)=rn23a*fvs(i)
         r23br(i)=rn23b*r00(i)
         r25rt(i)=rn25*rr0(i)*d2t
         r31r(i)=rn31*rr0(i)
         r32rt(i)=rn32*d2t*rrs(i)
        pt(i)=dpt(i,k)
        qv(i)=dqv(i,k)
        qc(i)=qcl(i,k)
        qr(i)=qrn(i,k)
        qi(i)=qci(i,k)
        qs(i)=qcs(i,k)
        qg(i)=qcg(i,k)
         if (qc(i) .le.  cmin1) qc(i)=0.0
         if (qr(i) .le.  cmin1) qr(i)=0.0
         if (qi(i) .le.  cmin1) qi(i)=0.0
         if (qs(i) .le.  cmin1) qs(i)=0.0
         if (qg(i) .le.  cmin1) qg(i)=0.0
        tair(i)=(pt(i)+tb0)*pi0(i)
        tairc(i)=tair(i)-t0
         zr(i)=zrr(i)
         zs(i)=zsr(i)
         zg(i)=zgr(i)
         vr(i)=0.0
         vs(i)=0.0
         vg(i)=0.0

         enddo
!JJS 10/7/2008     vvvvv
    IF (IWARM .EQ. 1) THEN
!JJS   for calculating processes related to warm rain only
!dir$ vector aligned
         DO i=1,irestrict
                qi(i)=0.0
                qs(i)=0.0
                qg(i)=0.0
                dep(i)=0.
                pint(i)=0.
                psdep(i)=0.
                pgdep(i)=0.
                dd1(i)=0.
                pgsub(i)=0.
                psmlt(i)=0.
                pgmlt(i)=0.
                pimlt(i)=0.
                psacw(i)=0.
                piacr(i)=0.
                psfw(i)=0.
                pgfr(i)=0.
                dgacw(i)=0.
                dgacr(i)=0.
                psacr(i)=0.
                wgacr(i)=0.
                pihom(i)=0.
                pidw(i)=0.

                if (qr(i) .gt. cmin1) then
                   dd(i)=r00(i)*qr(i)
                   y1(i)=dd(i)**.25
                   zr(i)=zrc/y1(i)
                   vr(i)=max(vrcf(i)*dd(i)**bwq, 0.)
                endif
!* 21 * PRAUT   AUTOCONVERSION OF QC TO QR                        **21**
!* 22 * PRACW : ACCRETION OF QC BY QR                             **22**
                pracw(i)=0.
                praut(i)=0.0
                pracw(i)=r22f(i)*qc(i)/zr(i)**bw3
                y1(i)=qc(i)-bnd3
                if (y1(i).gt.0.0) then
                    praut(i)=r00(i)*y1(i)*y1(i)/(1.2e-4+rn21/y1(i))
                 endif

!C********   HANDLING THE NEGATIVE CLOUD WATER (QC)    ******************
                 Y1(i)=QC(i)/D2T
                 PRAUT(i)=MIN(Y1(i), PRAUT(i))
                 PRACW(i)=MIN(Y1(i), PRACW(i))
                 Y1(i)=(PRAUT(i)+PRACW(i))*D2T
               
               if (qc(i) .lt. y1(i) .and. y1(i) .ge. cmin2) then
                   y2(i)=qc(i)/(y1(i)+cmin2)
                   praut(i)=praut(i)*y2(i)
                   pracw(i)=pracw(i)*y2(i)
                   qc(i)=0.0
               else
                  qc(i)=qc(i)-y1(i)
               endif
               PR(i)=(PRAUT(i)+PRACW(i))*D2T
               QR(i)=QR(i)+PR(i)
                        
!*****   TAO ET AL (1989) SATURATION TECHNIQUE  ***********************
           
           cnd(i)=0.0
           tair(i)=(pt(i)+tb0)*pi0(i)
              y1(i)=1./(tair(i)-c358)
              qsw(i)=rp0(i)*exp(c172-c409*y1(i))
              dd(i)=cp409(i)*y1(i)*y1(i)
              dm(i)=qv(i)+qb0-qsw(i)
              cnd(i)=dm(i)/(1.+avcp(i)*dd(i)*qsw(i))
!c    ******   condensation or evaporation of qc  ******
              cnd(i)=max(-qc(i), cnd(i))
                         pt(i)=pt(i)+avcp(i)*cnd(i)
             qv(i)=qv(i)-cnd(i)
                         qc(i)=qc(i)+cnd(i)

!C     ******   EVAPORATION   ******
!* 23 * ERN : EVAPORATION OF QR (SUBSATURATION)                   **23**
            ern(i)=0.0
            if(qr(i).gt.0.0) then
               tair(i)=(pt(i)+tb0)*pi0(i)
               rtair(i)=1./(tair(i)-c358)
               qsw(i)=rp0(i)*exp(c172-c409*rtair(i))
               ssw(i)=(qv(i)+qb0)/qsw(i)-1.0
               dm(i)=qv(i)+qb0-qsw(i)
               rsub1(i)=cv409(i)*qsw(i)*rtair(i)*rtair(i)
               dd1(i)=max(-dm(i)/(1.+rsub1(i)), 0.0)
               y1(i)=.78/zr(i)**2+r23af(i)*scv(i)/zr(i)**bwh5
               y2(i)=r23br(i)/(tca(i)*tair(i)**2)+1./(dwv(i) &
                       *qsw(i))
!cccc
               ern(i)=r23t*ssw(i)*y1(i)/y2(i)
               ern(i)=min(dd1(i),qr(i),max(ern(i),0.))
               pt(i)=pt(i)-avcp(i)*ern(i)
               qv(i)=qv(i)+ern(i)
               qr(i)=qr(i)-ern(i)
            endif
         end do ! i
       ELSE       ! part of if (iwarm.eq.1) then
!JJS 10/7/2008     ^^^^^
!dir$ vector aligned
         DO i=1,irestrict
!JJS   for calculating processes related to both ice and warm rain

!     ***   COMPUTE ZR,ZS,ZG,VR,VS,VG      *****************************

            if (qr(i) .gt. cmin1) then
               dd(i)=r00(i)*qr(i)
               y1(i)=dd(i)**.25
               zr(i)=zrc/y1(i)
               vr(i)=max(vrcf(i)*dd(i)**bwq, 0.)
            endif

            if (qs(i) .gt. cmin1) then
               dd(i)=r00(i)*qs(i)
               y1(i)=dd(i)**.25
               zs(i)=zsc/y1(i)
               vs(i)=max(vscf(i)*dd(i)**bsq, 0.)
            endif

            if (qg(i) .gt. cmin1) then
               dd(i)=r00(i)*qg(i)
               y1(i)=dd(i)**.25
               zg(i)=zgc/y1(i)
               if(ihail .eq. 1) then
                  vg(i)=max(vgcr(i)*dd(i)**bgq, 0.)
               else
                  vg(i)=max(vgcf(i)*dd(i)**bgq, 0.)
               endif
            endif

            if (qr(i) .le. cmin2) vr(i)=0.0
            if (qs(i) .le. cmin2) vs(i)=0.0
            if (qg(i) .le. cmin2) vg(i)=0.0

!     ******************************************************************
!     ***   Y1 : DYNAMIC VISCOSITY OF AIR (U)
!     ***   DWV : DIFFUSIVITY OF WATER VAPOR IN AIR (PI)
!     ***   TCA : THERMAL CONDUCTIVITY OF AIR (KA)
!     ***   Y2 : KINETIC VISCOSITY (V)
            y1(i)=c149*tair(i)**1.5/(tair(i)+120.)
            dwv(i)=dwvp(i)*tair(i)**1.81
            tca(i)=c141*y1(i)
            scv(i)=1./((rr0(i)*y1(i))**.1666667*dwv(i)**.3333333)
!*  1 * PSAUT : AUTOCONVERSION OF QI TO QS                        ***1**
!*  3 * PSACI : ACCRETION OF QI TO QS                             ***3**
!*  4 * PSACW : ACCRETION OF QC BY QS (RIMING) (QSACW FOR PSMLT)  ***4**
!*  5 * PRACI : ACCRETION OF QI BY QR                             ***5**
!*  6 * PIACR : ACCRETION OF QR OR QG BY QI                       ***6**
            psaut(i)=0.0
            psaci(i)=0.0
            praci(i)=0.0
            piacr(i)=0.0
            psacw(i)=0.0
            qsacw(i)=0.0
            dd(i)=1./zs(i)**bs3
            if (tair(i).lt.t0) then
               esi(i)=exp(.025*tairc(i))
               psaut(i)=r2is*max(rn1*esi(i)*(qi(i)-bnd1) ,0.0)
               psaci(i)=r2is*r3f(i)*esi(i)*qi(i)*dd(i)
!JJS 3/30/06
!    to cut water to snow accretion by half
               psacw(i)=r2is*0.5*r4f(i)*qc(i)*dd(i)
!JJS 3/30/06
               praci(i)=r2is*r5f(i)*qi(i)/zr(i)**bw3
               piacr(i)=r2is*r6f(i)*qi(i)*(zr(i)**(-bw6))
            else
               qsacw(i)=r2is*r4f(i)*qc(i)*dd(i)
            endif

!* 21 * PRAUT   AUTOCONVERSION OF QC TO QR                        **21**
!* 22 * PRACW : ACCRETION OF QC BY QR                             **22**

            pracw(i)=r22f(i)*qc(i)/zr(i)**bw3
            praut(i)=0.0
            y1(i)=qc(i)-bnd3
            if (y1(i).gt.0.0) then
               praut(i)=r00(i)*y1(i)*y1(i)/(1.2e-4+rn21/y1(i))
            endif

!* 12 * PSFW : BERGERON PROCESSES FOR QS (KOENING, 1971)          **12**
!* 13 * PSFI : BERGERON PROCESSES FOR QS                          **13**
            psfw(i)=0.0
            psfi(i)=0.0
            pidep(i)=0.0

            if(tair(i).lt.t0.and.qi(i).gt.cmin) then
               y1(i)=max( min(tairc(i), -1.), -31.)
               it(i)=int(abs(y1(i)))
               y1(i)=rn12a(it(i))
               y2(i)=rn12b(it(i))
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          psfw(i)=r2is* &
                    max(d2t*y1(i)*(y2(i)+r12r(i)*qc(i))*qi(i),0.0)
               rtair(i)=1./(tair(i)-c76)
               y2(i)=exp(c218-c580*rtair(i))
               qsi(i)=rp0(i)*y2(i)
               esi(i)=c610*y2(i)
               ssi(i)=(qv(i)+qb0)/qsi(i)-1.
               r_nci(i)=min(1.e-6*exp(-.46*tairc(i)),1.)
               dm(i)=max( (qv(i)+qb0-qsi(i)), 0.)
               rsub1(i)=cs580(i)*qsi(i)*rtair(i)*rtair(i)
               y3(i)=1./tair(i)
          dd(i)=y3(i)*(rn30a*y3(i)-rn30b)+rn30c*tair(i)/esi(i)
               y1(i)=206.18*ssi(i)/dd(i)
               pidep(i)=y1(i)*sqrt(r_nci(i)*qi(i)/r00(i))
               dep(i)=dm(i)/(1.+rsub1(i))/d2t
               if(dm(i).gt.cmin2) then
                  a2(i)=1.
                if(pidep(i).gt.dep(i).and.pidep(i).gt.cmin2) then
                     a2(i)=dep(i)/pidep(i)
                     pidep(i)=dep(i)
                endif
                  psfi(i)=r2is*a2(i)*.5*qi(i)*y1(i)/(sqrt(ami100) &
                          -sqrt(ami40))
                  elseif(dm(i).lt.-cmin2) then
!
!        SUBLIMATION TERMS USED ONLY WHEN SATURATION ADJUSTMENT FOR ICE
!        IS TURNED OFF
!
                  pidep(i)=0.
                  psfi(i)=0.
               else
                  pidep(i)=0.
                  psfi(i)=0.
               endif
            endif
!TTT***** QG=QG+MIN(PGDRY,PGWET)
!*  9 * PGACS : ACCRETION OF QS BY QG (DGACS,WGACS: DRY AND WET)  ***9**
!* 14 * DGACW : ACCRETION OF QC BY QG (QGACW FOR PGMLT)           **14**
!* 16 * DGACR : ACCRETION OF QR TO QG (QGACR FOR PGMLT)           **16**
            if(qc(i)+qr(i).lt.1.e-4) then
               ee1(i)=.01
              else
                 ee1(i)=1.
              endif
            ee2(i)=0.09
            egs(i)=ee1(i)*exp(ee2(i)*tairc(i))
            if (tair(i).ge.t0) egs(i)=1.0
            y1(i)=abs(vg(i)-vs(i))
            y2(i)=zs(i)*zg(i)
            y3(i)=5./y2(i)
            y4(i)=.08*y3(i)*y3(i)
            y5(i)=.05*y3(i)*y4(i)
            dd(i)=y1(i)*(y3(i)/zs(i)**5+y4(i)/zs(i)**3 &
                    +y5(i)/zs(i))
            pgacs(i)=r2ig*r2is*r9r(i)*egs(i)*dd(i)
!JJS 1/3/06 from Steve and Chunglin
            if (ihail.eq.1) then
               dgacs(i)=pgacs(i)
            else
               dgacs(i)=0.
            endif
!JJS 1/3/06 from Steve and Chunglin
            wgacs(i)=r2ig*r2is*r9r(i)*dd(i)
            y1(i)=1./zg(i)**bg3

            if(ihail .eq. 1) then
               dgacw(i)=r2ig*max(r14r(i)*qc(i)*y1(i), 0.0)
            else
               dgacw(i)=r2ig*max(r14f(i)*qc(i)*y1(i), 0.0)
            endif
         end do ! i
!dir$ vector aligned
         DO i=1,irestrict
            qgacw(i)=dgacw(i)
            y1(i)=abs(vg(i)-vr(i))
            y2(i)=zr(i)*zg(i)
            y3(i)=5./y2(i)
            y4(i)=.08*y3(i)*y3(i)
            y5(i)=.05*y3(i)*y4(i)
            dd(i)=r16r(i)*y1(i)*(y3(i)/zr(i)**5+y4(i)/zr(i)**3 &
                    +y5(i)/zr(i))
            dgacr(i)=r2ig*max(dd(i), 0.0)
            qgacr(i)=dgacr(i)

            if (tair(i).ge.t0) then
               dgacs(i)=0.0
               wgacs(i)=0.0
               dgacw(i)=0.0
               dgacr(i)=0.0
            else
               pgacs(i)=0.0
               qgacw(i)=0.0
               qgacr(i)=0.0
            endif
!*******PGDRY : DGACW+DGACI+DGACR+DGACS                           ******
!* 15 * DGACI : ACCRETION OF QI BY QG (WGACI FOR WET GROWTH)      **15**
!* 17 * PGWET : WET GROWTH OF QG                                  **17**
            dgaci(i)=0.0
            wgaci(i)=0.0
            pgwet(i)=0.0
         end do ! i
!dir$ vector aligned
         DO i=1,irestrict
            if (tair(i).lt.t0) then
               y1(i)=qi(i)/zg(i)**bg3
               if (ihail.eq.1) then
                  dgaci(i)=r2ig*r15r(i)*y1(i)
                  wgaci(i)=r2ig*r15ar(i)*y1(i)
               else
                   dgaci(i)=0.
                  wgaci(i)=r2ig*r15af(i)*y1(i)
               endif
!
               if (tairc(i).ge.-50.) then
!                if (alf+rn17c*tairc(i) .eq. 0.) then
!                   write(91,*) itimestep, i,j,k, alf, rn17c, tairc(i)
!                endif
                y1(i)=1./(alf+rn17c*tairc(i))
                if (ihail.eq.1) then
                   y3(i)=.78/zg(i)**2+r17aq(i)*scv(i)/zg(i)**bgh5
                else
                   y3(i)=.78/zg(i)**2+r17as(i)*scv(i)/zg(i)**bgh5
                endif
                y4(i)=alvr(i)*dwv(i)*(rp0(i)-(qv(i)+qb0)) &
                        -tca(i)*tairc(i)
                dd(i)=y1(i)*(r17r(i)*y4(i)*y3(i) &
                       +(wgaci(i)+wgacs(i))*(alf+rn17b*tairc(i)))
                pgwet(i)=r2ig*max(dd(i), 0.0)
               endif
            endif
         end do ! i
!dir$ vector aligned
         DO i=1,irestrict
!********   HANDLING THE NEGATIVE CLOUD WATER (QC)    ******************
!********   HANDLING THE NEGATIVE CLOUD ICE (QI)      ******************
            y1(i)=qc(i)/d2t
            psacw(i)=min(y1(i), psacw(i))
            praut(i)=min(y1(i), praut(i))
            pracw(i)=min(y1(i), pracw(i))
            psfw(i)= min(y1(i), psfw(i))
            dgacw(i)=min(y1(i), dgacw(i))
            qsacw(i)=min(y1(i), qsacw(i))
            qgacw(i)=min(y1(i), qgacw(i))

            y1(i)=(psacw(i)+praut(i)+pracw(i)+psfw(i) &
                    +dgacw(i)+qsacw(i)+qgacw(i))*d2t
            qc(i)=qc(i)-y1(i)

            if (qc(i) .lt. 0.0) then
               a1(i)=1.
               if (y1(i) .ne. 0.0) a1(i)=qc(i)/y1(i)+1.
               psacw(i)=psacw(i)*a1(i)
               praut(i)=praut(i)*a1(i)
               pracw(i)=pracw(i)*a1(i)
               psfw(i)=psfw(i)*a1(i)
               dgacw(i)=dgacw(i)*a1(i)
               qsacw(i)=qsacw(i)*a1(i)
               qgacw(i)=qgacw(i)*a1(i)
               qc(i)=0.0
            endif
!c
!
!******** SHED PROCESS (WGACR=PGWET-DGACW-WGACI-WGACS)
!c
         end do
!dir$ vector aligned
         DO i=1,irestrict
            wgacr(i)=pgwet(i)-dgacw(i)-wgaci(i)-wgacs(i)
            y2(i)=dgacw(i)+dgaci(i)+dgacr(i)+dgacs(i)
            if (pgwet(i).ge.y2(i)) then
               wgacr(i)=0.0
               wgaci(i)=0.0
               wgacs(i)=0.0
            else
               dgacr(i)=0.0
               dgaci(i)=0.0
               dgacs(i)=0.0
            endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
            y1(i)=qi(i)/d2t
            psaut(i)=min(y1(i), psaut(i))
            psaci(i)=min(y1(i), psaci(i))
            praci(i)=min(y1(i), praci(i))
            psfi(i)= min(y1(i), psfi(i))
            dgaci(i)=min(y1(i), dgaci(i))
            wgaci(i)=min(y1(i), wgaci(i))
!
            y2(i)=(psaut(i)+psaci(i)+praci(i)+psfi(i) &
                   +dgaci(i)+wgaci(i))*d2t
            qi(i)=qi(i)-y2(i)+pidep(i)*d2t

            if (qi(i).lt.0.0) then
               a2(i)=1.
               if (y2(i) .ne. 0.0) a2(i)=qi(i)/y2(i)+1.
               psaut(i)=psaut(i)*a2(i)
               psaci(i)=psaci(i)*a2(i)
               praci(i)=praci(i)*a2(i)
               psfi(i)=psfi(i)*a2(i)
               dgaci(i)=dgaci(i)*a2(i)
               wgaci(i)=wgaci(i)*a2(i)
               qi(i)=0.0
            endif
!
            dlt3(i)=0.0
            dlt2(i)=0.0
            if (tair(i).lt.t0) then
               if (qr(i).lt.1.e-4) then
                  dlt3(i)=1.0
                  dlt2(i)=1.0
               endif
               if (qs(i).ge.1.e-4) then
                  dlt2(i)=0.0
               endif
            endif

            if (ice2 .eq. 1) then
                  dlt3(i)=1.0
                  dlt2(i)=1.0
            endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         end do
!dir$ vector aligned
         DO i=1,irestrict
            pr(i)=(qsacw(i)+praut(i)+pracw(i)+qgacw(i))*d2t
            ps(i)=(psaut(i)+psaci(i)+psacw(i)+psfw(i) &
                    +psfi(i)+dlt3(i)*praci(i))*d2t
            pg(i)=((1.-dlt3(i))*praci(i)+dgaci(i)+wgaci(i) &
                    +dgacw(i))*d2t

!*  7 * PRACS : ACCRETION OF QS BY QR                             ***7**
!*  8 * PSACR : ACCRETION OF QR BY QS (QSACR FOR PSMLT)           ***8**
            y1(i)=abs(vr(i)-vs(i))
            y2(i)=zr(i)*zs(i)
            y3(i)=5./y2(i)
            y4(i)=.08*y3(i)*y3(i)
            y5(i)=.05*y3(i)*y4(i)
            pracs(i)=r2ig*r2is*r7r(i)*y1(i)*(y3(i)/zs(i)**5 &
                      +y4(i)/zs(i)**3+y5(i)/zs(i))
            psacr(i)=r2is*r8r(i)*y1(i)*(y3(i)/zr(i)**5 &
                      +y4(i)/zr(i)**3+y5(i)/zr(i))
            qsacr(i)=psacr(i)

            if (tair(i).ge.t0) then
               pracs(i)=0.0
               psacr(i)=0.0
            else
               qsacr(i)=0.0
            endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!*  2 * PGAUT : AUTOCONVERSION OF QS TO QG                        ***2**
!* 18 * PGFR : FREEZING OF QR TO QG                               **18**

            pgaut(i)=0.0
            pgfr(i)=0.0

            if (tair(i) .lt. t0) then
               y2(i)=exp(rn18a*(t0-tair(i)))
!JJS              PGFR(i)=r2ig*max(R18R(i)*(Y2(i)-1.)/ZR(i)**7., 0.0)
!        modify to prevent underflow on some computers (JD)
               temp(i) = 1./zr(i)
               temp(i) = temp(i)*temp(i)*temp(i)*temp(i)*temp(i)*temp(i)*temp(i)
               pgfr(i)=r2ig*max(r18r(i)*(y2(i)-1.)* &
                                    temp(i), 0.0)
            endif

!********   HANDLING THE NEGATIVE RAIN WATER (QR)    *******************
!********   HANDLING THE NEGATIVE SNOW (QS)          *******************

            y1(i)=qr(i)/d2t
            y2(i)=-qg(i)/d2t
            piacr(i)=min(y1(i), piacr(i))
            dgacr(i)=min(y1(i), dgacr(i))
            wgacr(i)=min(y1(i), wgacr(i))
            wgacr(i)=max(y2(i), wgacr(i))
            psacr(i)=min(y1(i), psacr(i))
            pgfr(i)= min(y1(i), pgfr(i))
            del(i)=0.
            if(wgacr(i) .lt. 0.) del(i)=1.
            y1(i)=(piacr(i)+dgacr(i)+(1.-del(i))*wgacr(i) &
                    +psacr(i)+pgfr(i))*d2t
            qr(i)=qr(i)+pr(i)-y1(i)-del(i)*wgacr(i)*d2t
            if (qr(i) .lt. 0.0) then
               a1(i)=1.
               if(y1(i) .ne. 0.) a1(i)=qr(i)/y1(i)+1.
               piacr(i)=piacr(i)*a1(i)
               dgacr(i)=dgacr(i)*a1(i)
               if (wgacr(i).gt.0.) wgacr(i)=wgacr(i)*a1(i)
               pgfr(i)=pgfr(i)*a1(i)
               psacr(i)=psacr(i)*a1(i)
               qr(i)=0.0
            endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            prn(i)=d2t*((1.-dlt3(i))*piacr(i)+dgacr(i) &
                     +wgacr(i)+(1.-dlt2(i))*psacr(i)+pgfr(i))
            ps(i)=ps(i)+d2t*(dlt3(i)*piacr(i) &
                    +dlt2(i)*psacr(i))
            pracs(i)=(1.-dlt2(i))*pracs(i)
            y1(i)=qs(i)/d2t
            pgacs(i)=min(y1(i), pgacs(i))
            dgacs(i)=min(y1(i), dgacs(i))
            wgacs(i)=min(y1(i), wgacs(i))
            pgaut(i)=min(y1(i), pgaut(i))
            pracs(i)=min(y1(i), pracs(i))
            psn(i)=d2t*(pgacs(i)+dgacs(i)+wgacs(i) &
                     +pgaut(i)+pracs(i))
            qs(i)=qs(i)+ps(i)-psn(i)

            if (qs(i).lt.0.0) then
               a2(i)=1.
               if (psn(i) .ne. 0.0) a2(i)=qs(i)/psn(i)+1.
               pgacs(i)=pgacs(i)*a2(i)
               dgacs(i)=dgacs(i)*a2(i)
               wgacs(i)=wgacs(i)*a2(i)
               pgaut(i)=pgaut(i)*a2(i)
               pracs(i)=pracs(i)*a2(i)
               psn(i)=psn(i)*a2(i)
               qs(i)=0.0
            endif
            y2(i)=d2t*(psacw(i)+psfw(i)+dgacw(i)+piacr(i) &
                    +dgacr(i)+wgacr(i)+psacr(i)+pgfr(i))
            pt(i)=pt(i)+afcp(i)*y2(i)
            qg(i)=qg(i)+pg(i)+prn(i)+psn(i)

!* 11 * PSMLT : MELTING OF QS                                     **11**
!* 19 * PGMLT : MELTING OF QG TO QR                               **19**
            psmlt(i)=0.0
            pgmlt(i)=0.0
            tair(i)=(pt(i)+tb0)*pi0(i)

            if (tair(i).ge.t0) then
               tairc(i)=tair(i)-t0
               y1(i)=tca(i)*tairc(i)-alvr(i)*dwv(i) &
                               *(rp0(i)-(qv(i)+qb0))
               y2(i)=.78/zs(i)**2+r101f(i)*scv(i)/zs(i)**bsh5
               dd(i)=r11rt(i)*y1(i)*y2(i)+r11at*tairc(i) &
                       *(qsacw(i)+qsacr(i))
               psmlt(i)=r2is*max(0.0, min(dd(i), qs(i)))

               if(ihail.eq.1) then
                  y3(i)=.78/zg(i)**2+r19aq(i)*scv(i)/zg(i)**bgh5
               else
                  y3(i)=.78/zg(i)**2+r19as(i)*scv(i)/zg(i)**bgh5
               endif

               dd1(i)=r19rt(i)*y1(i)*y3(i)+r19bt*tairc(i) &
                        *(qgacw(i)+qgacr(i))
               pgmlt(i)=r2ig*max(0.0, min(dd1(i), qg(i)))
               pt(i)=pt(i)-afcp(i)*(psmlt(i)+pgmlt(i))
               qr(i)=qr(i)+psmlt(i)+pgmlt(i)
               qs(i)=qs(i)-psmlt(i)
               qg(i)=qg(i)-pgmlt(i)
            endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!* 24 * PIHOM : HOMOGENEOUS FREEZING OF QC TO QI (T < T00)        **24**
!* 25 * PIDW : DEPOSITION GROWTH OF QC TO QI ( T0 < T <= T00)     **25**
!* 26 * PIMLT : MELTING OF QI TO QC (T >= T0)                     **26**

            if (qc(i).le.cmin1) qc(i)=0.0
            if (qi(i).le.cmin1) qi(i)=0.0
            tair(i)=(pt(i)+tb0)*pi0(i)

            if(tair(i).le.t00) then
               pihom(i)=qc(i)
            else
               pihom(i)=0.0
            endif
            if(tair(i).ge.t0) then
               pimlt(i)=qi(i)
            else
               pimlt(i)=0.0
            endif
            pidw(i)=0.0

            if (tair(i).lt.t0 .and. tair(i).gt.t00) then
               tairc(i)=tair(i)-t0
               y1(i)=max( min(tairc(i), -1.), -31.)
               it(i)=int(abs(y1(i)))
               y2(i)=aa1(it(i))
               y3(i)=aa2(it(i))
               y4(i)=exp(abs(beta*tairc(i)))
               y5(i)=(r00(i)*qi(i)/(r25a*y4(i)))**y3(i)
               pidw(i)=min(r25rt(i)*y2(i)*y4(i)*y5(i), qc(i))
            endif
         end do
!dir$ vector aligned
         DO i=1,irestrict
            y1(i)=pihom(i)-pimlt(i)+pidw(i)
            pt(i)=pt(i)+afcp(i)*y1(i)+ascp(i)*(pidep(i))*d2t
            qv(i)=qv(i)-(pidep(i))*d2t
            qc(i)=qc(i)-y1(i)
            qi(i)=qi(i)+y1(i)

!* 31 * PINT  : INITIATION OF QI                                  **31**
!* 32 * PIDEP : DEPOSITION OF QI                                  **32**
!
!     CALCULATION OF PINT USES DIFFERENT VALUES OF THE INTERCEPT AND SLOPE FOR
!     THE FLETCHER EQUATION. ALSO, ONLY INITIATE MORE ICE IF THE NEW NUMBER
!     CONCENTRATION EXCEEDS THAT ALREADY PRESENT.
!* 31 * pint  : initiation of qi                                  **31**
!* 32 * pidep : deposition of qi                                  **32**
           pint(i)=0.0
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if ( itaobraun.eq.1 ) then
            tair(i)=(pt(i)+tb0)*pi0(i)
            if (tair(i) .lt. t0) then
              if (qi(i) .le. cmin2) qi(i)=0.
               tairc(i)=tair(i)-t0
               rtair(i)=1./(tair(i)-c76)
               y2(i)=exp(c218-c580*rtair(i))
              qsi(i)=rp0(i)*y2(i)
               esi(i)=c610*y2(i)
              ssi(i)=(qv(i)+qb0)/qsi(i)-1.
                        ami50=3.76e-8

!ccif ( itaobraun.eq.1 ) --> betah=0.5*beta=-.46*0.5=-0.23;   cn0=1.e-6
!ccif ( itaobraun.eq.0 ) --> betah=0.5*beta=-.6*0.5=-0.30;    cn0=1.e-8

             y1(i)=1./tair(i)

!cc insert a restriction on ice collection that ice collection
!cc should be stopped at -30 c (with cn0=1.e-6, beta=-.46)

             tairccri(i)=tairc(i)          ! in degree c
             if(tairccri(i).le.-30.) tairccri(i)=-30.

             y2(i)=exp(betah*tairccri(i))
             y3(i)=sqrt(qi(i))
             dd(i)=y1(i)*(rn10a*y1(i)-rn10b)+rn10c*tair(i) &
                                                /esi(i)
          pidep(i)=max(r32rt(i)*ssi(i)*y2(i)*y3(i)/dd(i), 0.e0)

           r_nci(i)=min(cn0*exp(beta*tairc(i)),1.)

           dd(i)=max(1.e-9*r_nci(i)/r00(i)-qi(i)*1.e-9/ami50, 0.)
                dm(i)=max( (qv(i)+qb0-qsi(i)), 0.0)
                rsub1(i)=cs580(i)*qsi(i)*rtair(i)*rtair(i)
              dep(i)=dm(i)/(1.+rsub1(i))
              pint(i)=max(min(dd(i), dm(i)), 0.)

              pint(i)=min(pint(i)+pidep(i), dep(i))

               if (pint(i) .le. cmin2) pint(i)=0.
              pt(i)=pt(i)+ascp(i)*pint(i)
              qv(i)=qv(i)-pint(i)
              qi(i)=qi(i)+pint(i)
           endif
        endif  ! if ( itaobraun.eq.1 )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if ( itaobraun.eq.0 ) then
             tair(i)=(pt(i)+tb0)*pi0(i)
             if (tair(i) .lt. t0) then
               if (qi(i) .le. cmin1) qi(i)=0.
               tairc(i)=tair(i)-t0
               dd(i)=r31r(i)*exp(beta*tairc(i))
               rtair(i)=1./(tair(i)-c76)
               y2(i)=exp(c218-c580*rtair(i))
               qsi(i)=rp0(i)*y2(i)
               esi(i)=c610*y2(i)
               ssi(i)=(qv(i)+qb0)/qsi(i)-1.
               dm(i)=max( (qv(i)+qb0-qsi(i)), 0.)
               rsub1(i)=cs580(i)*qsi(i)*rtair(i)*rtair(i)
               dep(i)=dm(i)/(1.+rsub1(i))
              pint(i)=max(min(dd(i), dm(i)), 0.)
               y1(i)=1./tair(i)
               y2(i)=exp(betah*tairc(i))
               y3(i)=sqrt(qi(i))
               dd(i)=y1(i)*(rn10a*y1(i)-rn10b) &
                     +rn10c*tair(i)/esi(i)
             pidep(i)=max(r32rt(i)*ssi(i)*y2(i)*y3(i)/dd(i), 0.)
              pint(i)=pint(i)+pidep(i)
              pint(i)=min(pint(i),dep(i))
             pt(i)=pt(i)+ascp(i)*pint(i)
             qv(i)=qv(i)-pint(i)
             qi(i)=qi(i)+pint(i)
            endif
        endif  ! if ( itaobraun.eq.0 )
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!*****   TAO ET AL (1989) SATURATION TECHNIQUE  ***********************

         if (new_ice_sat .eq. 0) then
               tair(i)=(pt(i)+tb0)*pi0(i)
               cnd(i)=rt0*(tair(i)-t00)
               dep(i)=rt0*(t0-tair(i))
               y1(i)=1./(tair(i)-c358)
               y2(i)=1./(tair(i)-c76)
               qsw(i)=rp0(i)*exp(c172-c409*y1(i))
               qsi(i)=rp0(i)*exp(c218-c580*y2(i))
               dd(i)=cp409(i)*y1(i)*y1(i)
               dd1(i)=cp580(i)*y2(i)*y2(i)
               if (qc(i).le.cmin) qc(i)=cmin
               if (qi(i).le.cmin) qi(i)=cmin
               if (tair(i).ge.t0) then
                  dep(i)=0.0
                  cnd(i)=1.
                  qi(i)=0.0
               endif

               if (tair(i).lt.t00) then
                  cnd(i)=0.0
                  dep(i)=1.
                  qc(i)=0.0
               endif

               y5(i)=avcp(i)*cnd(i)+ascp(i)*dep(i)
               y1(i)=qc(i)*qsw(i)/(qc(i)+qi(i))
               y2(i)=qi(i)*qsi(i)/(qc(i)+qi(i))
               y4(i)=dd(i)*y1(i)+dd1(i)*y2(i)
               qvs(i)=y1(i)+y2(i)
               rsub1(i)=(qv(i)+qb0-qvs(i))/(1.+y4(i)*y5(i))
               cnd(i)=cnd(i)*rsub1(i)
               dep(i)=dep(i)*rsub1(i)
               if (qc(i).le.cmin) qc(i)=0.
               if (qi(i).le.cmin) qi(i)=0.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c    ******   condensation or evaporation of qc  ******

               cnd(i)=max(-qc(i),cnd(i))

!c    ******   deposition or sublimation of qi    ******

               dep(i)=max(-qi(i),dep(i))

               pt(i)=pt(i)+avcp(i)*cnd(i)+ascp(i)*dep(i)
               qv(i)=qv(i)-cnd(i)-dep(i)
               qc(i)=qc(i)+cnd(i)
               qi(i)=qi(i)+dep(i)
         endif

         if (new_ice_sat .eq. 1) then
               tair(i)=(pt(i)+tb0)*pi0(i)
               cnd(i)=rt0*(tair(i)-t00)
               dep(i)=rt0*(t0-tair(i))
               y1(i)=1./(tair(i)-c358)
               y2(i)=1./(tair(i)-c76)
               qsw(i)=rp0(i)*exp(c172-c409*y1(i))
               qsi(i)=rp0(i)*exp(c218-c580*y2(i))
               dd(i)=cp409(i)*y1(i)*y1(i)
               dd1(i)=cp580(i)*y2(i)*y2(i)
               y5(i)=avcp(i)*cnd(i)+ascp(i)*dep(i)
               y1(i)=rt0*(tair(i)-t00)*qsw(i)
               y2(i)=rt0*(t0-tair(i))*qsi(i)
               if (tair(i).ge.t0) then
                  dep(i)=0.0
                  cnd(i)=1.
                  y2(i)=0.
                  y1(i)=qsw(i)
               endif
               if (tair(i).lt.t00) then
                  cnd(i)=0.0
                  dep(i)=1.
                  y2(i)=qsi(i)
                  y1(i)=0.
               endif
               y4(i)=dd(i)*y1(i)+dd1(i)*y2(i)
               qvs(i)=y1(i)+y2(i)
               rsub1(i)=(qv(i)+qb0-qvs(i))/(1.+y4(i)*y5(i))
               cnd(i)=cnd(i)*rsub1(i)
               dep(i)=dep(i)*rsub1(i)

!C    ******   CONDENSATION OR EVAPORATION OF QC  ******
               cnd(i)=max(-qc(i),cnd(i))

!C    ******   DEPOSITION OR SUBLIMATION OF QI    ******
               dep(i)=max(-qi(i),dep(i))
               pt(i)=pt(i)+avcp(i)*cnd(i)+ascp(i)*dep(i)
               qv(i)=qv(i)-cnd(i)-dep(i)
               qc(i)=qc(i)+cnd(i)
               qi(i)=qi(i)+dep(i)
         endif

!c
!
          if (new_ice_sat .eq. 2) then
          dep(i)=0.0
          cnd(i)=0.0
          tair(i)=(pt(i)+tb0)*pi0(i)
          if (tair(i) .ge. 253.16) then
              y1(i)=1./(tair(i)-c358)
              qsw(i)=rp0(i)*exp(c172-c409*y1(i))
              dd(i)=cp409(i)*y1(i)*y1(i)
              dm(i)=qv(i)+qb0-qsw(i)
              cnd(i)=dm(i)/(1.+avcp(i)*dd(i)*qsw(i))
!c    ******   condensation or evaporation of qc  ******
              cnd(i)=max(-qc(i), cnd(i))
             pt(i)=pt(i)+avcp(i)*cnd(i)
             qv(i)=qv(i)-cnd(i)
             qc(i)=qc(i)+cnd(i)
         endif
          if (tair(i) .le. 258.16) then
           y2(i)=1./(tair(i)-c76)
           qsi(i)=rp0(i)*exp(c218-c580*y2(i))
          dd1(i)=cp580(i)*y2(i)*y2(i)
         dep(i)=(qv(i)+qb0-qsi(i))/(1.+ascp(i)*dd1(i)*qsi(i))
!c    ******   deposition or sublimation of qi    ******
             dep(i)=max(-qi(i),dep(i))
             pt(i)=pt(i)+ascp(i)*dep(i)
             qv(i)=qv(i)-dep(i)
             qi(i)=qi(i)+dep(i)
         endif
      endif

!c
!
!* 10 * PSDEP : DEPOSITION OR SUBLIMATION OF QS                   **10**
!* 20 * PGSUB : SUBLIMATION OF QG                                 **20**
            psdep(i)=0.0
            pgdep(i)=0.0
            pssub(i)=0.0
            pgsub(i)=0.0
            tair(i)=(pt(i)+tb0)*pi0(i)

            if(tair(i).lt.t0) then
               if(qs(i).lt.cmin1) qs(i)=0.0
               if(qg(i).lt.cmin1) qg(i)=0.0
               rtair(i)=1./(tair(i)-c76)
               qsi(i)=rp0(i)*exp(c218-c580*rtair(i))
               ssi(i)=(qv(i)+qb0)/qsi(i)-1.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
               y1(i)=r10ar(i)/(tca(i)*tair(i)**2)+1./(dwv(i) &
                      *qsi(i))
               y2(i)=.78/zs(i)**2+r101f(i)*scv(i)/zs(i)**bsh5
               psdep(i)=r10t*ssi(i)*y2(i)/y1(i)
               pssub(i)=psdep(i)
               psdep(i)=r2is*max(psdep(i), 0.)
               pssub(i)=r2is*max(-qs(i), min(pssub(i), 0.))

               if(ihail.eq.1) then
                  y2(i)=.78/zg(i)**2+r20bq(i)*scv(i)/zg(i)**bgh5
               else
                  y2(i)=.78/zg(i)**2+r20bs(i)*scv(i)/zg(i)**bgh5
               endif

               pgsub(i)=r2ig*r20t*ssi(i)*y2(i)/y1(i)
               dm(i)=qv(i)+qb0-qsi(i)
               rsub1(i)=cs580(i)*qsi(i)*rtair(i)*rtair(i)

!     ********   DEPOSITION OR SUBLIMATION OF QS  **********************

               y1(i)=dm(i)/(1.+rsub1(i))
               psdep(i)=r2is*min(psdep(i),max(y1(i),0.))
               y2(i)=min(y1(i),0.)
               pssub(i)=r2is*max(pssub(i),y2(i))

!     ********   SUBLIMATION OF QG   ***********************************

               dd(i)=max((-y2(i)-qs(i)), 0.)
              pgsub(i)=r2ig*min(dd(i), qg(i), max(pgsub(i),0.))

               if(qc(i)+qi(i).gt.1.e-5) then
                  dlt1(i)=1.
               else
                  dlt1(i)=0.
               endif

               psdep(i)=dlt1(i)*psdep(i)
               pssub(i)=(1.-dlt1(i))*pssub(i)
               pgsub(i)=(1.-dlt1(i))*pgsub(i)

               pt(i)=pt(i)+ascp(i)*(psdep(i)+pssub(i)-pgsub(i))
               qv(i)=qv(i)+pgsub(i)-pssub(i)-psdep(i)
               qs(i)=qs(i)+psdep(i)+pssub(i)
               qg(i)=qg(i)-pgsub(i)
            endif

!* 23 * ERN : EVAPORATION OF QR (SUBSATURATION)                   **23**

            ern(i)=0.0

            if(qr(i).gt.0.0) then
               tair(i)=(pt(i)+tb0)*pi0(i)
               rtair(i)=1./(tair(i)-c358)
               qsw(i)=rp0(i)*exp(c172-c409*rtair(i))
               ssw(i)=(qv(i)+qb0)/qsw(i)-1.0
               dm(i)=qv(i)+qb0-qsw(i)
               rsub1(i)=cv409(i)*qsw(i)*rtair(i)*rtair(i)
               dd1(i)=max(-dm(i)/(1.+rsub1(i)), 0.0)
               y1(i)=.78/zr(i)**2+r23af(i)*scv(i)/zr(i)**bwh5
               y2(i)=r23br(i)/(tca(i)*tair(i)**2)+1./(dwv(i) &
                       *qsw(i))
!cccc
               ern(i)=r23t*ssw(i)*y1(i)/y2(i)
               ern(i)=min(dd1(i),qr(i),max(ern(i),0.))
               pt(i)=pt(i)-avcp(i)*ern(i)
               qv(i)=qv(i)+ern(i)
               qr(i)=qr(i)-ern(i)
            endif

         enddo ! i
!JJS 10/7/2008     vvvvv
    ENDIF    ! part of if (iwarm.eq.1) then
!JJS 10/7/2008     ^^^^^
!dir$ vector aligned
         DO i=1,irestrict
            if (qc(i) .le. cmin1) qc(i)=0.
            if (qr(i) .le. cmin1) qr(i)=0.
            if (qi(i) .le. cmin1) qi(i)=0.
            if (qs(i) .le. cmin1) qs(i)=0.
            if (qg(i) .le. cmin1) qg(i)=0.
            dpt(i,k)=pt(i)
            dqv(i,k)=qv(i)
            qcl(i,k)=qc(i)
            qrn(i,k)=qr(i)
            qci(i,k)=qi(i)
            qcs(i,k)=qs(i)
            qcg(i,k)=qg(i)
            dd(i)=max(-cnd(i), 0.)
            cnd(i)=max(cnd(i), 0.)
            dd1(i)=max(-dep(i), 0.)+pidep(i)*d2t
            dep(i)=max(dep(i), 0.)
 end do ! i
 end do ! k

!JJS  ****************************************************************
!JJS  convert from GCE grid back to WRF grid
      do k=kts,kte
!dir$ vector aligned
         DO i=1,irestrict
         ptwrf(i,k) = dpt(i,k)
         qvwrf(i,k) = dqv(i,k)
         qlwrf(i,k) = qcl(i,k)
         qrwrf(i,k) = qrn(i,k)
         qiwrf(i,k) = qci(i,k)
         qswrf(i,k) = qcs(i,k)
         qgwrf(i,k) = qcg(i,k)
         enddo !i
      enddo !k

!     ****************************************************************

      return
 END SUBROUTINE saticel_s

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!JJS      real function GAMMLN (xx)
  real function gammagce (xx)
!**********************************************************************
  real*8 cof(6),stp,half,one,fpf,x,tmp,ser
  data cof,stp /  76.18009173,-86.50532033,24.01409822, &
     -1.231739516,.120858003e-2,-.536382e-5, 2.50662827465 /
  data half,one,fpf / .5, 1., 5.5 /
!
      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*log(tmp)-tmp
      ser=one
      do  j=1,6
         x=x+one
        ser=ser+cof(j)/x
      enddo !j
      gammln=tmp+log(stp*ser)
!JJS
      gammagce=exp(gammln)
!JJS
      return
 END FUNCTION gammagce

END MODULE  module_mp_gsfcgce

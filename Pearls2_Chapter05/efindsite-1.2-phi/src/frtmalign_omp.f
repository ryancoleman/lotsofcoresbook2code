c===============================================================================
c         ___________.__            .____________.__  __          
c     ____\_   _____/|__| ____    __| _/   _____/|__|/  |_  ____  
c   _/ __ \|    __)  |  |/    \  / __ |\_____  \ |  \   __\/ __ \ 
c   \  ___/|     \   |  |   |  \/ /_/ |/        \|  ||  | \  ___/ 
c    \___  >___  /   |__|___|  /\____ /_______  /|__||__|  \___  >
c        \/    \/            \/      \/       \/               \/ 
c
c                                                  
c   eFindSite - ligand-binding site prediction from meta-threading
c
c   Computational Systems Biology Group
c   Department of Biological Sciences
c   Center for Computation & Technology
c   Louisiana State University
c   407 Choppin Hall, Baton Rouge, LA 70803, USA
c
c   http://www.brylinski.org
c
c   Report bugs to michal@brylinski.org wfeinstein@lsu.edu
c
c   Copyright 2013 Michal Brylinski, Wei Feinstein
c
c   This file is part of eFindSite.
c
c   eFindSite is free software: you can redistribute it and/or modify
c   it under the terms of the GNU General Public License as published by
c   the Free Software Foundation, either version 3 of the License, or
c   (at your option) any later version.
c
c   eFindSite is distributed in the hope that it will be useful,
c   but WITHOUT ANY WARRANTY; without even the implied warranty of
c   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
c   GNU General Public License for more details.
c
c   You should have received a copy of the GNU General Public License
c   along with eFindSite. If not, see <http://www.gnu.org/licenses/>.
c
c===============================================================================
c
c  Fr-TM-align converted into a subroutine

      subroutine frtmalign(tm_sco1,tm_sco2,tm_sco3,tm_sco4,tm_len1a,
     &                     tm_len1b,tm_seq1,tm_seq2,tm_align1,
     &                     tm_res1a,tm_res1b,tm_xyz1a,tm_xyz1b,tm_tt1,
     &                     tm_uu1,tm_lnorm)
      
      parameter (maxres=1000)           ! no. of residues
      parameter (maxlen=2*maxres)       ! for alignment
      parameter (npr=2)                 ! no. of proteins to align
      
      double precision tm_sco2,tm_sco3,tm_sco4,
     &                 tm_xyz1a(0:2,0:maxres-1),
     &                 tm_xyz1b(0:2,0:maxres-1),
     &                 tm_tt1(0:2),tm_uu1(0:2,0:2)
      
      integer tm_sco1,tm_len1a,tm_len1b,tm_lnorm,
     &        tm_align1(0:maxres-1),tm_res1a(0:maxres-1),
     &        tm_res1b(0:maxres-1)
      
      character*1 tm_seq1(0:maxres-1),tm_seq2(0:maxres-1)
      
      character*3 resa(maxres,npr),resn(maxres,0:1)
      character*1 seqa(maxres,npr),seqn(maxres,0:1)
      
      dimension invmap0(maxres),ires(maxres,0:1)
      dimension xtm1(maxres),ytm1(maxres),ztm1(maxres)
      dimension xtm2(maxres),ytm2(maxres),ztm2(maxres)
      dimension m1(maxres),m2(maxres)
      dimension xyz(3,maxres,npr),length(npr)
      
      common /coord/ xa(3,maxres,0:1)
      common /length/ nseq1,nseq2
      common /pdbinfo/ ires1(maxres,npr),resa,seqa
      common /secstr/ isec(maxres),jsec(maxres)     !secondary structure
      common /n1n2/ n1(maxres),n2(maxres)
      common /dinfo/ d8
      common /stru/xt(maxres),yt(maxres),zt(maxres),
     &             xb(maxres),yb(maxres),zb(maxres)
      
!$omp threadprivate(/coord/)
!$omp threadprivate(/length/)
!$omp threadprivate(/pdbinfo/)
!$omp threadprivate(/secstr/)
!$omp threadprivate(/n1n2/)
!$omp threadprivate(/dinfo/)
!$omp threadprivate(/stru/)
      
      double precision r_1(3,maxres),r_2(3,maxres),w(maxres)
      double precision u(3,3),t(3),rms !armsd is real
      data w /maxres*1.0/
      
      irun_type=0  ! TYPE of Algorithm used for Fragment alignment
                   ! 0 is slow version
                   ! 1 is fast version
      
******************************************************************
****    main program starts
******************************************************************
      
      length(1)=tm_len1a
      length(2)=tm_len1b
      
      do i=0,tm_len1a-1
       seqa(i+1,1)=tm_seq1(i)
       ires1(i+1,1)=tm_res1a(i)
       do j=1,3
        xyz(j,i+1,1)=real(tm_xyz1a(j-1,i))
       enddo
      enddo
      
      do i=0,tm_len1b-1
       seqa(i+1,2)=tm_seq2(i)
       ires1(i+1,2)=tm_res1b(i)
       do j=1,3
        xyz(j,i+1,2)=real(tm_xyz1b(j-1,i))
       enddo
      enddo
      
cc initialize the coordinates xa() for superposition

      do 1166 ipdb=1,1     ! main cycle
       do i=0,1
        ik=i+1
        do j=1,length(ik)
           do k=1,3
              xa(k,j,i)=xyz(k,j,ik)   
           enddo
              seqn(j,i)=seqa(j,ik)
              resn(j,i)=resa(j,ik)
              ires(j,i)=ires1(j,ik)
        enddo
       enddo
       nseq1=length(1)
       nseq2=length(2)
       
cc fixing d0 for search ONLY with d0 fixed for small protein
      d0_min=0.5
      aminlen=min(nseq1,nseq2)
      d8=1.5*aminlen**0.3+3.5      !remove pairs with dis>d8 during search 
      if(aminlen.gt.15) then
        d0=1.24*(aminlen-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
        d0=d0_min
      endif
      nseq=max(nseq1,nseq2)
       do i=1,nseq
        n1(i)=i
        n2(i)=i
       enddo
      d0_search=d0
        
      if (d0_search.gt.8.0)d0_search=8.0
      if (d0_search.lt.3.0)d0_search=4.5

      call super_align (d0,d0_search,invmap0,irun_type)

cc resuperose to find residues d < d8
      if (d0_search.lt.4.5)d0_search=4.5
      n_al=0
      do j=1,nseq2
         if(invmap0(j).gt.0)then
            i=invmap0(j)
            n_al=n_al+1
            xtm1(n_al)=xa(1,i,0)
            ytm1(n_al)=xa(2,i,0)
            ztm1(n_al)=xa(3,i,0)
            xtm2(n_al)=xa(1,j,1)
            ytm2(n_al)=xa(2,j,1)
            ztm2(n_al)=xa(3,j,1)
            m1(n_al)=i          !for recording residue order
            m2(n_al)=j
         endif
      enddo
      d0_input=d0
      isearch=2 
      call TMsearch (d0_input,d0_search,n_al,xtm1,ytm1,ztm1,n1,n_al,
     &     xtm2,ytm2,ztm2,n2,TM,Rcomm,Lcomm,isearch) !TM-score with dis<d8 only

cccc for Output TM-score set d0 
      d0_min=0.5                !for output
      anseq=tm_lnorm               !length for defining final TMscore
      if(anseq.gt.15)then
         d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
         d0=d0_min
      endif
      if(d0.lt.d0_min)d0=d0_min
      d0_search=d0
      if (d0_search.gt.8.0)d0_search=8.0
      if (d0_search.lt.4.5)d0_search=4.5

cccc remove dis>d8 in normal TM-score calculation for final report
      j=0
      n_eq=0
      do i=1,n_al
         dis2=sqrt((xtm1(i)-xtm2(i))**2+(ytm1(i)-ytm2(i))**2+
     &        (ztm1(i)-ztm2(i))**2)
         if(dis2.le.d8)then
            j=j+1
            xtm1(j)=xtm1(i)
            ytm1(j)=ytm1(i)
            ztm1(j)=ztm1(i)
            xtm2(j)=xtm2(i)
            ytm2(j)=ytm2(i)
            ztm2(j)=ztm2(i)
            m1(j)=m1(i)
            m2(j)=m2(i)
            if(seqn(m1(i),0).eq.seqn(m2(i),1) )then
               n_eq=n_eq+1
            endif
         endif
      enddo
      seq_id=float(n_eq)/(n_al+0.00000001)

      n8_al=j
      d0_input=d0
      isearch=3
      call TMsearch(d0_input,d0_search,n8_al,xtm1,ytm1,ztm1,n1,n8_al,
     &     xtm2,ytm2,ztm2,n2,TM8,Rcomm,Lcomm,isearch) !normal TMscore
      rmsd8_al=Rcomm
      TM8=TM8/anseq       !TM-score after cutoff
      
 1166 continue

cccc Rotation matrix

       L=0
        do i=1,n8_al
          k=m1(i)
          L=L+1
          r_1(1,L)=xa(1,k,0)
          r_1(2,L)=xa(2,k,0)
          r_1(3,L)=xa(3,k,0)
          r_2(1,L)=xtm1(i)
          r_2(2,L)=ytm1(i)
          r_2(3,L)=ztm1(i)
        enddo
        if(L.gt.3)then
           call u3b(w,r_1,r_2,L,1,rms,u,t,ier) !u rotate r_1 to r_2
        endif

c Get results back
      
      tm_sco1=n8_al
      tm_sco2=rmsd8_al
      tm_sco3=TM8
      tm_sco4=seq_id
      
      do i=1,3
       tm_tt1(i-1)=t(i)
       do j=1,3
        tm_uu1(j-1,i-1)=u(i,j)
       enddo
      enddo
      
      do i=1,tm_len1b
       tm_align1(i-1)=-1
      enddo
      
      do i=1,n8_al
       tm_align1(m2(i)-1)=m1(i)-1
      enddo
      
      end


******************************************************************
****            Subroutines
******************************************************************
cccc making initial superposition

      subroutine super_align(dx,dxs,invmap0,ir_type)
      parameter (maxres=1000)
      dimension invmap0(maxres)

      common /coord/ xa(3,maxres,0:1)
      common /length/ nseq1,nseq2
      common /secstr/ isec(maxres),jsec(maxres)     !secondary structure
       
!$omp threadprivate(/coord/)
!$omp threadprivate(/length/)
!$omp threadprivate(/secstr/)
       
       call assignssp   ! secondary structure assignment

       call fragscan(dx,dxs,invmap0,ir_type,atm)

      return
      end

cccc for getting the best alignment
      subroutine getbest(aTM,invmapi,aTMmax,invmapr)
      parameter (maxres=1000)
      dimension invmapi(maxres),invmapr(maxres)
      common /length/ nseq1,nseq2

!$omp threadprivate(/length/)

        if (aTM.gt.aTMmax) then
            aTMmax=aTM
                do j=1,nseq2
                    invmapr(j)=invmapi(j)
                enddo
        endif
      return
      end

cccc making secondary structure assignment
      subroutine assignssp
      parameter (maxres=1000)
      common /coord/ xa(3,maxres,0:1)
      common /length/ nseq1,nseq2
      common /secstr/ isec(maxres),jsec(maxres)     !secondary structure

!$omp threadprivate(/coord/)
!$omp threadprivate(/length/)
!$omp threadprivate(/secstr/)

       do i=1,nseq1
          isec(i)=1
          j1=i-2
          j2=i-1
          j3=i
          j4=i+1
          j5=i+2
          if(j1.ge.1.and.j5.le.nseq1)then
             dis13=diszy(0,j1,j3)
             dis14=diszy(0,j1,j4)
             dis15=diszy(0,j1,j5)
             dis24=diszy(0,j2,j4)
             dis25=diszy(0,j2,j5)
             dis35=diszy(0,j3,j5)
             isec(i)=make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
          endif
       enddo
       do i=1,nseq2
          jsec(i)=1
          j1=i-2
          j2=i-1
          j3=i
          j4=i+1
          j5=i+2
          if(j1.ge.1.and.j5.le.nseq2)then
             dis13=diszy(1,j1,j3)
             dis14=diszy(1,j1,j4)
             dis15=diszy(1,j1,j5)
             dis24=diszy(1,j2,j4)
             dis25=diszy(1,j2,j5)
             dis35=diszy(1,j3,j5)
             jsec(i)=make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
          endif
       enddo
       call smooth               !smooth the assignment

      return
      end

cccc make secondary str
      function make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
      make_sec=1
      delta=2.1
      if(abs(dis15-6.37).lt.delta)then
         if(abs(dis14-5.18).lt.delta)then
            if(abs(dis25-5.18).lt.delta)then
               if(abs(dis13-5.45).lt.delta)then
                  if(abs(dis24-5.45).lt.delta)then
                     if(abs(dis35-5.45).lt.delta)then
                        make_sec=2 !helix
                        return
                     endif
                  endif
               endif
            endif
         endif
      endif
      delta=1.42
      if(abs(dis15-13).lt.delta)then
         if(abs(dis14-10.4).lt.delta)then
            if(abs(dis25-10.4).lt.delta)then
               if(abs(dis13-6.1).lt.delta)then
                  if(abs(dis24-6.1).lt.delta)then
                     if(abs(dis35-6.1).lt.delta)then
                        make_sec=4 !strand
                        return
                     endif
                  endif
               endif
            endif
         endif
      endif
      if(dis15.lt.8)then
         make_sec=3
      endif

      return
      end

cccc smooth the secondary structure assignment
      subroutine smooth
      parameter (maxres=1000)
      common /length/ nseq1,nseq2
      common /secstr/ isec(maxres),jsec(maxres)     !secondary structure

!$omp threadprivate(/length/)
!$omp threadprivate(/secstr/)

***   smooth single -------------->
***   --x-- => -----
       do i=3,nseq1
          if(isec(i).eq.2.or.isec(i).eq.4)then
             j=isec(i)
             if(isec(i-2).ne.j)then
                if(isec(i-1).ne.j)then
                   if(isec(i+1).ne.j)then
                      if(isec(i+1).ne.j)then
                         isec(i)=1
                      endif
                   endif
                endif
             endif
          endif
       enddo
       do i=3,nseq2
          if(jsec(i).eq.2.or.jsec(i).eq.4)then
             j=jsec(i)
             if(jsec(i-2).ne.j)then
                if(jsec(i-1).ne.j)then
                   if(jsec(i+1).ne.j)then
                      if(jsec(i+1).ne.j)then
                         jsec(i)=1
                      endif
                   endif
                endif
             endif
          endif
       enddo

***   smooth double -------------->
***   --xx-- => ------
       do i=1,nseq1-5
          if(isec(i).ne.2)then
            if(isec(i+1).ne.2)then
                if(isec(i+2).eq.2)then
                    if(isec(i+3).eq.2)then
                        if(isec(i+4).ne.2)then
                            if(isec(i+5).ne.2)then
                                isec(i+2)=1
                                isec(i+3)=1
                            endif
                        endif
                    endif
                endif
            endif
          endif
 
          if(isec(i).ne.4)then
            if(isec(i+1).ne.4)then
                if(isec(i+2).eq.4)then
                    if(isec(i+3).eq.4)then
                        if(isec(i+4).ne.4)then
                            if(isec(i+5).ne.4)then
                                isec(i+2)=1
                                isec(i+3)=1
                            endif
                        endif
                    endif
                endif
            endif
          endif
       enddo
       do i=1,nseq2-5
          if(jsec(i).ne.2)then
            if(jsec(i+1).ne.2)then
                if(jsec(i+2).eq.2)then
                    if(jsec(i+3).eq.2)then
                        if(jsec(i+4).ne.2)then
                            if(jsec(i+5).ne.2)then
                                jsec(i+2)=1
                                jsec(i+3)=1
                            endif
                        endif
                    endif
                endif
            endif
          endif
 
          if(jsec(i).ne.4)then
            if(jsec(i+1).ne.4)then
                if(jsec(i+2).eq.4)then
                    if(jsec(i+3).eq.4)then
                        if(jsec(i+4).ne.4)then
                            if(jsec(i+5).ne.4)then
                                jsec(i+2)=1
                                jsec(i+3)=1
                            endif
                        endif
                    endif
                endif
            endif
          endif
       enddo

***   connect -------------->
***   x-x => xxx
       do i=1,nseq1-2
          if(isec(i).eq.2)then
            if(isec(i+1).ne.2)then
                if(isec(i+2).eq.2)then
                    isec(i+1)=2
                endif
            endif
          endif
 
          if(isec(i).eq.4)then
            if(isec(i+1).ne.4)then
                if(isec(i+2).eq.4)then
                    isec(i+1)=4
                endif
            endif
          endif
       enddo
       do i=1,nseq2-2
          if(jsec(i).eq.2)then
            if(jsec(i+1).ne.2)then
                if(jsec(i+2).eq.2)then
                    jsec(i+1)=2
                endif
            endif
          endif
 
          if(jsec(i).eq.4)then
            if(jsec(i+1).ne.4)then
                if(jsec(i+2).eq.4)then
                    jsec(i+1)=4
                endif
            endif
          endif
       enddo
 
      return
      end

      
cccc distance calculation

      function diszy(i,i1,i2)
       parameter (maxres=1000)
       common /coord/ xa(3,maxres,0:1)

!$omp threadprivate(/coord/)

        diszy=sqrt((xa(1,i1,i)-xa(1,i2,i))**2
     &     +(xa(2,i1,i)-xa(2,i2,i))**2
     &     +(xa(3,i1,i)-xa(3,i2,i))**2)
       return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc  Subroutine to calculate GL-score and secondary structure identity
cccc  for a given fragment length
cccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_nGL (ist1,istp1,ist2,istp2,GL,ssa,ssan)  
      parameter (maxres=1000)           ! no. of residues  
      dimension xo1(maxres),yo1(maxres),zo1(maxres)
      dimension xo2(maxres),yo2(maxres),zo2(maxres)
      dimension dis2(maxres)
      common /coord/ xa(3,maxres,0:1)
      common /length/ nseq1,nseq2
      common /secstr/ isec(maxres),jsec(maxres)

!$omp threadprivate(/coord/)
!$omp threadprivate(/length/)
!$omp threadprivate(/secstr/)

cccc  RMSD:
      double precision r_1(3,maxres),r_2(3,maxres),w(maxres)
      double precision u(3,3),t(3),rms !armsd is real
      data w /maxres*1.0/

       da=0.5
       d02=(da)**2
       GL=0.0
       ssa=0.0

       idiff1=(istp1-ist1)+1
       idiff2=(istp2-ist2)+1
       istp2_mod=istp2
    
       if (idiff1.lt.idiff2) then   ! if first fragment is less than 
                                    ! second one make second one equal to first fragment
        istp2_mod=ist2+idiff1-1
       endif
        
       i=ist1
       n_al=0
       iss=0
       issn=0
       do j=ist2,istp2_mod
         n_al=n_al+1
         r_1(1,n_al)=xa(1,i,0)
         r_1(2,n_al)=xa(2,i,0)
         r_1(3,n_al)=xa(3,i,0)
         r_2(1,n_al)=xa(1,j,1)
         r_2(2,n_al)=xa(2,j,1)
         r_2(3,n_al)=xa(3,j,1)
         xo1(n_al)=xa(1,i,0)
         yo1(n_al)=xa(2,i,0)
         zo1(n_al)=xa(3,i,0)
         xo2(n_al)=xa(1,j,1)
         yo2(n_al)=xa(2,j,1)
         zo2(n_al)=xa(3,j,1)
         if (isec(i).eq.jsec(j).and.isec(i).gt.1)iss=iss+1
         if (isec(i).eq.jsec(j))issn=issn+1
         i=i+1
       enddo

       ssa=float(iss)/float(n_al)
       ssan=float(issn)/float(n_al)

       call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier)  !u rotate r_1 to r_2
       GL1=0.0
       do i=1,n_al
          dis2(i)=0.0
          xx=real(t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i))
          yy=real(t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i))
          zz=real(t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i))
          dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
          GL1=GL1+(1/(1+dis2(i)/(d02)))
       enddo
       GL1=GL1/float(n_al)

       dxs=4.5
       d002=(dxs**2)

   20  continue
       j=0
       do i=1,n_al
         if (dis2(i).le.d002) then
            j=j+1
            r_1(1,j)=xo1(i)
            r_1(2,j)=yo1(i)
            r_1(3,j)=zo1(i)
            r_2(1,j)=xo2(i)
            r_2(2,j)=yo2(i)
            r_2(3,j)=zo2(i)
        endif
       enddo
       if (j.le.3) then
        d002=d002+0.5
       goto 20
      endif
       L=j

       call u3b (w,r_1,r_2,L,1,rms,u,t,ier)  !u rotate r_1 to r_2

       GL2=0.0
       do i=1,n_al
          dis2(i)=0.0
          xx=real(t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i))
          yy=real(t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i))
          zz=real(t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i))
          dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
          GL2=GL2+(1/(1+dis2(i)/(d02)))
       enddo
       GL2=GL2/float(n_al)

   21  continue
       j=0
       do i=1,n_al
         if (dis2(i).le.d002) then
            j=j+1
            r_1(1,j)=xo1(i)
            r_1(2,j)=yo1(i)
            r_1(3,j)=zo1(i)
            r_2(1,j)=xo2(i)
            r_2(2,j)=yo2(i)
            r_2(3,j)=zo2(i)
        endif
       enddo
       if (j.le.3) then
        d002=d002+0.5
       goto 21
      endif
       L=j

       call u3b (w,r_1,r_2,L,1,rms,u,t,ier)  !u rotate r_1 to r_2

       GL3=0.0
       do i=1,n_al
          dis2(i)=0.0
          xx=real(t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i))
          yy=real(t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i))
          zz=real(t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i))
          dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
          GL3=GL3+(1/(1+dis2(i)/(d02)))
       enddo
       GL3=GL3/float(n_al)
       
        GL=GL1
        if (GL2.gt.GL) GL=GL2
        if (GL3.gt.GL) GL=GL3
        
       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc
cccc Subroutine for aligning the fragments based on GL-score matrix or
cccc SSA+GL-score matrix.
cccc
cccc A variation Waterman-Eggert algorithm is applied to generate
cccc suboptimal alignments after getting the best optimal alignments.
cccc The step of generating suboptimal alignment:
cccc 1. Recompute the F matrix with not allowing match for those position,
cccc    which are aligned in the previous alignment.
cccc 2. Backtrace the F matrix generating the alignment
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       subroutine fragdp (nfr1,nfr2,smtx,gap,ifrali,ntot)
       parameter (maxfr=1200)
       dimension smtx(maxfr,maxfr),ifrali(100,maxfr)
       dimension idir(0:maxfr,0:maxfr),lcmap(maxfr)
       dimension val(0:maxfr,0:maxfr)
       dimension lf(0:maxfr,0:maxfr)

       ! initialize the F matrix
         val(0,0)=0.0
         do i=1,nfr1
          idir(i,0)=0
          val(i,0)=0.0
         enddo

         do j=1,nfr2
          idir(0,j)=0
          val(0,j)=0.0
          lcmap(j)=-1
         enddo
        
       ! compute F matrix
            do j=1,nfr2
                do i=1,nfr1
                 d=val(i-1,j-1)+smtx(i,j)
                 h=val(i-1,j)+gap
                 v=val(i,j-1)+gap 
                 lf(i,j)=0
             
                 if (d.ge.h.and.d.ge.v) then
                     val(i,j)=d
                     idir(i,j)=1
                 else
                     if (h.ge.v) then
                         val(i,j)=h
                         idir(i,j)=2
                     else
                         val(i,j)=v
                         idir(i,j)=3
                     endif
                 endif
                enddo
            enddo
        
        amax=0.0
        call traceback (nfr1,nfr2,val,idir,lf,lcmap,1,amax)
        ntot=ntot+1
        do ii=1,nfr2
            ifrali(ntot,ii)=lcmap(ii)
            lcmap(ii)=-1
        enddo
        
        call recomputefmatrix (nfr1,nfr2,smtx,lf,val,gap,idir)
        call traceback (nfr1,nfr2,val,idir,lf,lcmap,2,amax)

        ntot=ntot+1
        do ii=1,nfr2
            ifrali(ntot,ii)=lcmap(ii)
            lcmap(ii)=-1
        enddo

        call recomputefmatrix (nfr1,nfr2,smtx,lf,val,gap,idir)
        call traceback (nfr1,nfr2,val,idir,lf,lcmap,3,amax)

        ntot=ntot+1
        do ii=1,nfr2
            ifrali(ntot,ii)=lcmap(ii)
            lcmap(ii)=-1
        enddo

        call recomputefmatrix (nfr1,nfr2,smtx,lf,val,gap,idir)
        call traceback (nfr1,nfr2,val,idir,lf,lcmap,4,amax)

        ntot=ntot+1
        do ii=1,nfr2
            ifrali(ntot,ii)=lcmap(ii)
            lcmap(ii)=-1
        enddo

       return
       end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc Routine to filter the redundant fragment alignments
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine filter (ifrali,ntot,nfr2,nret,jali) 
        parameter (maxfr=1200)
        dimension ifrali (100,maxfr),jali(100,maxfr),istate(100)

        do ii=1,ntot
            istate(ii)=1 ! all are unique alignment
        enddo
        
         do ii=1,ntot-1
            do jj=ii+1,ntot
             if (istate(jj).eq.1) then
               nlc=0
               do kk=1,nfr2
                if (ifrali(ii,kk).eq.ifrali(jj,kk)) nlc=nlc+1
               enddo
               if (nlc.eq.nfr2) istate(jj)=0
             endif
            enddo
         enddo

         nret=0
         do ii=1,ntot
            if (istate(ii).eq.1) then
                nret=nret+1
                do jj=1,nfr2
                    jali(nret,jj)=ifrali(ii,jj)
                enddo
            endif
         enddo
        
        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc Subroutine for calculating the F matrix and returning the traceback matrix
cccc 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       subroutine recomputefmatrix (nfr1,nfr2,smtx,lf,val,gap,idir)
       parameter (maxfr=1200)
       dimension smtx(maxfr,maxfr)
       dimension idir(0:maxfr,0:maxfr)
       dimension val(0:maxfr,0:maxfr)
       dimension lf(0:maxfr,0:maxfr)

       ! initialize F matrix
         val(0,0)=0.0
         do i=1,nfr1
            idir(i,0)=0
            val(i,0)=0.0
         enddo

         do j=1,nfr2
            idir(0,j)=0
            val(0,j)=0.0
         enddo

      ! fill F matrix using a modified recursion relation

        do j=1,nfr2
             do i=1,nfr1
              if (lf(i,j).eq.0) then
                 d=val(i-1,j-1)+smtx(i,j)
              else
                 d=0.0
              endif
             
              h=val(i-1,j)+gap
              v=val(i,j-1)+gap 
             
               if (d.ge.h.and.d.ge.v) then
                 val(i,j)=d
                 idir(i,j)=1
               else
                 if (h.ge.v) then
                     val(i,j)=h
                     idir(i,j)=2
                 else
                     val(i,j)=v
                     idir(i,j)=3
                 endif
               endif
             enddo
        enddo

       return
       end

cccc Subroutine to traceback the filled F matrix 
cccc 
       subroutine traceback (nfr1,nfr2,val,idir,lf,lcmap,iround,amax)
       parameter (maxfr=1200)
       dimension idir(0:maxfr,0:maxfr),val(0:maxfr,0:maxfr)
       dimension lf(0:maxfr,0:maxfr),lcmap(maxfr)
       
        ! Traceback GLOBAL optimal/suboptimal alignment with no penalty for end gaps
        ! find the highest score on the last row and column and start
        ! alignment from that position

        i=nfr1
        j=nfr2
        amax2=0.0
        if (iround.eq.1) then
            amax=40000.00
        endif

        do ii=1,nfr1
            if (val(ii,nfr2).gt.amax2.and.val(ii,nfr2).lt.amax) then
                amax2=val(ii,nfr2)
                i=ii
            endif
        enddo

        do jj=1,nfr2
            if (val(nfr1,jj).gt.amax2.and.val(nfr1,jj).lt.amax) then
                amax2=val(nfr1,jj)
                j=jj
                i=nfr1
            endif
        enddo

        if (iround.eq.1) then
            amax=amax2
        endif
    
        do while ( (i.gt.0) .and. (j.gt.0) )
            if ( idir(i,j) .eq. 1 ) then
                lcmap(j)=i
                lf(i,j)=1
                i=i-1
                j=j-1
            elseif (idir(i,j) .eq. 2) then
                i=i-1
            elseif (idir(i,j) .eq. 3) then
                j=j-1
            endif
        enddo

       return
       end

************************************************************************
**
** The protein is divided into set of non-overlapping fragments of a
** given length and we would comapre one fragment against all others
** using cheap version of TM-score (GL-score). This score is finally 
** used for doing a repeated match DP to get many fragment alignment.
** The aligned sets of fragment is joined to make the initial seed 
** alignment. The initial alignment is expanded with a routine of
** extension using DP and TM-score. Now, we are using a new DP routine.
**
************************************************************************

       subroutine fragscan (dx,dxs,mapr,irt,TMmax)
       parameter (maxres=1000)
       parameter (maxfr=1200)   ! max. no. of fragments
       parameter (mfrali=100)   ! max. no. of aligned fragments

       dimension mapr(maxres)
       dimension dpst(maxfr,maxfr),sspt(maxfr,maxfr),iali(100,maxfr)
       dimension smtx(maxfr,maxfr),ifrali(100,maxfr),gap(4),
     &           sspta(maxfr,maxfr)
       
       common /length/ nseq1,nseq2
       common /secstr/ isec(maxres),jsec(maxres)
       
!$omp threadprivate(/length/)
!$omp threadprivate(/secstr/)
       
       fr=float(min(nseq1,nseq2))/float(max(nseq1,nseq2))
       minseq=min(nseq1,nseq2)
       
       if (minseq.lt.100) then
        ipp=8
       else
        ipp=12
       endif

       lenfr=ipp
       nfr1=nseq1/lenfr ! no. of fragments for seq1
       nfr2=nseq2/lenfr ! no. of fragments for seq2
       glmax=0.0
       ntot=0
       slmax=0.0
       slmaxn=0.0
    
       do ii=1,nfr2
        istart2=(ii-1)*lenfr+1
        iend2=istart2+lenfr-1
        if (ii.eq.nfr2) iend2=nseq2

            do jj=1,nfr1
             istart1=(jj-1)*lenfr+1
             iend1=istart1+lenfr-1
             if (jj.eq.nfr1) iend1=nseq1

             call get_nGL(istart1,iend1,istart2,iend2,GL,ssa,ssna)
              
              if (GL.gt.glmax) then
                glmax=GL
              endif
               if (ii.ne.nfr2.and.jj.ne.nfr1) then
                if (ssa.gt.slmax) slmax=ssa
                if (ssna.gt.slmaxn) slmaxn=ssna
               endif

              dpst(jj,ii)=GL
              sspt(jj,ii)=ssa
              sspta(jj,ii)=ssna
              ntot=ntot+1
            enddo
       enddo

       gap(1)=-0.6
       gap(2)=-0.1
       gap(3)=0.0
        
       tmmax=0.0
       ntot=0 
        do ii=1,nfr1
            do jj=1,nfr2
                smtx(ii,jj)=dpst(ii,jj)/glmax
            enddo
        enddo

       do ipp=1,3
        gap_open=gap(ipp)
        call fragdp (nfr1,nfr2,smtx,gap_open,ifrali,ntot) ! returns ntot 
       enddo

       if (slmax .lt. 0.3) then
        do ii=1,nfr1
         do jj=1,nfr2
            sspt(ii,jj)=sspta(ii,jj)
         enddo
        enddo
        slmax=slmaxn
       endif
       
        do ii=1,nfr1
            do jj=1,nfr2
                smtx(ii,jj)=0.5*(sspt(ii,jj)/slmax)+
     &                      0.5*(dpst(ii,jj)/glmax)
            enddo
        enddo

       do ipp=1,3
        gap_open=gap(ipp)
        call fragdp (nfr1,nfr2,smtx,gap_open,ifrali,ntot) ! returns ntot 
       enddo

       call filter (ifrali,ntot,nfr2,nali,iali)
       call calbesttm (dx,dxs,nfr1,nfr2,lenfr,nali,iali,irt,mapr,tmmax)

       end

************************************************************************
*****          Get the alignment from the fragments and do heuristic
*****          iteration
************************************************************************

       subroutine calbesttm (dx,dxs,nfr1,nfr2,nlen,ntot,ifrali,
     &                       irt,mapr,tmmax)
       parameter (maxres=1000)
       parameter (maxfr=1200)
       dimension ifrali(100,maxfr),mapr(maxres),invmap(maxres)
       dimension lcmap(maxfr),invmapr(maxres)

       common /length/ nseq1,nseq2
       
!$omp threadprivate(/length/)
       
       tmmax=0.0

        do iali=1,ntot
            do ii=1,nfr2
                lcmap(ii)=ifrali(iali,ii)
            enddo

            call fillinvmap (invmap)

            do ii=1,nfr2
                ifr1=lcmap(ii)
                ifr2=ii
               if (ifr1.gt.0) then
                istart1=(ifr1-1)*nlen+1
                istart2=(ifr2-1)*nlen+1
                iend1=istart1+nlen-1
                iend2=istart2+nlen-1
                if (ifr1.eq.nfr1) iend1=nseq1
                if (ifr2.eq.nfr2) iend2=nseq2
                if ((iend1-istart1).lt.(iend2-istart2)) then
                    iend2=(iend1-istart1)+istart2
                endif
                 itmp=istart1
                 do jj=istart2,iend2
                    invmap(jj)=itmp
                    itmp=itmp+1
                 enddo
               endif
            enddo
            
            tminp=tmmax
            
            if (irt .eq. 0) then
                call caltmsc (dx,dxs,invmap,tminp,invmapr,b1tmmax)
            elseif (irt.eq.1) then
                call caltmsc_fast (dx,dxs,invmap,tminp,invmapr,b1tmmax)
            endif

            if (b1tmmax.gt.tmmax) then
                do ii=1,nseq2
                    mapr(ii)=invmapr(ii)
                enddo
                tmmax=b1tmmax
            endif

        enddo ! for total no. of alignments

       return
       end

*********************************************************
**      Routine for slow calculation of best alignment
**      Scans many alignments!!
*********************************************************

       subroutine caltmsc (dx,dxs,map,tmch,mapr,tmscm)
       parameter (maxres=1000)
       dimension map(maxres),mapr(maxres),score(maxres,maxres)
       dimension scorea(maxres,maxres),gap(2),invmap_i(maxres)
       dimension invmapr(maxres),scoret(maxres,maxres)
    
       common /length/ nseq1,nseq2
       common /secstr/ isec(maxres),jsec(maxres)
       
!$omp threadprivate(/length/)
!$omp threadprivate(/secstr/)
       
       gap(1)=-1.0
       gap(2)=0.0
       
       tmscm=0.0
       call get_initial3 (dx,map,scoret,scorea)
       tminp_lc=tmch
       istart=1
       iend=2

       if (tmch .gt. 0.50) then
            iend=1
       endif

       do ii=istart,iend
         gap_open=gap(ii)
         call dp (scorea,gap_open,invmap_i,dpout)
         call get_score (dx,dxs,invmap_i,score,TM)
         tTM=TM
         if (tTM .gt. tminp_lc) then
            call make_iter(dx,dxs,score,12,invmap_i,tTM,invmapr)
            call getbest (tTM,invmapr,tmscm,mapr)
            if (tmscm .gt. tminp_lc) tminp_lc=tmscm
         else
            call make_iter(dx,dxs,score,7,invmap_i,tTM,invmapr)
            call getbest (tTM,invmapr,tmscm,mapr)
            if (tmscm .gt. tminp_lc) tminp_lc=tmscm
         endif
       enddo

       if (tmscm .gt. tminp_lc) tminp_lc=tmscm

       if (tminp_lc .gt. 0.50) then
            istart=2
       endif

       do ii=istart,iend
         gap_open=gap(ii)
         call dp (scoret,gap_open,invmap_i,dpout)
         call get_score (dx,dxs,invmap_i,score,TM)
         tTM=TM
         if (tTM .gt. tmscm) then
            call make_iter(dx,dxs,score,12,invmap_i,tTM,invmapr)
            call getbest (tTM,invmapr,tmscm,mapr)
         else
            call make_iter(dx,dxs,score,7,invmap_i,tTM,invmapr)
            call getbest (tTM,invmapr,tmscm,mapr)
         endif
       enddo

       if (tmscm .gt. tminp_lc) tminp_lc=tmscm

       if (tminp_lc .lt. 0.50) then
       
       call get_score (dx,dxs,map,scorea,TM)
       do ii=2,2
         gap_open=gap(ii)
         call dp (scorea,gap_open,invmap_i,dpout)
         call get_score (dx,dxs,invmap_i,score,TM)
         tTM=TM
         if (tTM .gt. tmscm) then
            call make_iter(dx,dxs,score,12,invmap_i,tTM,invmapr)
            call getbest (tTM,invmapr,tmscm,mapr)
         else
            call make_iter(dx,dxs,score,6,invmap_i,tTM,invmapr)
            call getbest (tTM,invmapr,tmscm,mapr)
         endif
       enddo
       
       endif

       return
       end

*********************************************************
**** General routine to fill the starting alignment !!
*********************************************************
       subroutine fillinvmap (invmapr)
       parameter (maxres=1000)
       dimension invmapr(maxres)
       common /length/ nseq1,nseq2

!$omp threadprivate(/length/)

        do jj=1,nseq2
            invmapr(jj)=-1
        enddo

       return
       end

*********************************************************
**      Routine for fast calculation of best alignment
**      Scans less number of alignments!!
*********************************************************

       subroutine caltmsc_fast (dx,dxs,map,tmch,mapr,tmscm)
       parameter (maxres=1000)
       dimension map(maxres),mapr(maxres),score(maxres,maxres)
       dimension scorea(maxres,maxres),gap(2),invmap_i(maxres)
       dimension invmapr(maxres),scoret(maxres,maxres)

       common /length/ nseq1,nseq2
       common /secstr/ isec(maxres),jsec(maxres)

!$omp threadprivate(/length/)
!$omp threadprivate(/secstr/)

       gap(1)=-1.0
       gap(2)=0.0

       tmscm=0.0
       call get_initial3 (dx,map,scoret,scorea)
       tminp_lc=tmch

       do ii=1,1
         gap_open=gap(ii)
         call dp (scorea,gap_open,invmap_i,dpout)
         call get_score (dx,dxs,invmap_i,score,TM)
         tTM=TM
         if (tTM .gt. tminp_lc) then
            call make_iter(dx,dxs,score,14,invmap_i,tTM,invmapr)
            call getbest (tTM,invmapr,tmscm,mapr)
            if (tmscm .gt. tminp_lc) tminp_lc=tmscm
         else
            call make_iter(dx,dxs,score,7,invmap_i,tTM,invmapr)
            call getbest (tTM,invmapr,tmscm,mapr)
            if (tmscm .gt. tminp_lc) tminp_lc=tmscm
         endif
       enddo

       if (tminp_lc .lt. tmscm ) tminp_lc=tmscm
                    
       do ii=2,2
         gap_open=gap(ii)
         call dp (scoret,gap_open,invmap_i,dpout)
         call get_score (dx,dxs,invmap_i,score,TM)
         tTM=TM
         if (tTM .gt. tmscm) then
            call make_iter(dx,dxs,score,14,invmap_i,tTM,invmapr)
            call getbest (tTM,invmapr,tmscm,mapr)
         else
            call make_iter(dx,dxs,score,7,invmap_i,tTM,invmapr)
            call getbest (tTM,invmapr,tmscm,mapr)
         endif
       enddo

       if (tminp_lc .lt. tmscm ) tminp_lc=tmscm

       if (tminp_lc .lt. 0.50) then
        call get_score (dx,dxs,map,scorea,TM)
            do ii=2,2
                gap_open=gap(ii)
                call dp (scorea,gap_open,invmap_i,dpout)
                call get_score (dx,dxs,invmap_i,score,TM)
                tTM=TM
                if (tTM .gt. tmscm) then
                    call make_iter(dx,dxs,score,14,invmap_i,
     &                             tTM,invmapr)
                    call getbest (tTM,invmapr,tmscm,mapr)
                else
                    call make_iter(dx,dxs,score,7,invmap_i,tTM,invmapr)
                    call getbest (tTM,invmapr,tmscm,mapr)
                endif
            enddo
       endif

       return
       end
 
******************************************************************


cccc Subroutines for Fr-TM-align 
cccc get_initial3,dp,make_iter,get_score1,get_score

*************************************************************
** 	DP routine (score matrix, gap open, alignment)
**  returns invmap_i
** Fix: In the Iteration routine if the number of residues aligned is
**      less than 5 residues the method aborts checking for any 
**      improvement in the alignment.
*************************************************************
      subroutine dp (score,gap_open,invmap_i,dpout)
      
      parameter (maxres=1000)             ! no. of residues  
      dimension score(maxres,maxres),invmap_i(maxres)
      logical*1 dir
      dimension dir(0:maxres,0:maxres),val(0:maxres,0:maxres)
      real h,v

      common /length/ nseq1,nseq2

!$omp threadprivate(/length/)

cccc initialize matrix
      do ii=1,nseq2
        invmap_i(ii)=-1
      enddo

       val(0,0)=0.0
       do i=1,nseq1
         dir(i,0)=.false.
         val(i,0)=0.0
       enddo
       do j=1,nseq2
         dir(0,j)=.false.
         val(0,j)=0.0
         invmap_i(j)=-1
       enddo

cccc fill the matrix and path

       do j=1,nseq2
        do i=1,nseq1
           d=val(i-1,j-1)+score(i,j)
           h=val(i-1,j)
           v=val(i,j-1)

            if(dir(i-1,j))h=h+gap_open
            if(dir(i,j-1))v=v+gap_open
 
                if((d.ge.h).and.(d.ge.v)) then
                    dir(i,j)=.true.
                    val(i,j)=d
                else
                    dir(i,j)=.false.
                        if(v.ge.h)then
                            val(i,j)=v
                        else
                            val(i,j)=h
                        end if
                endif
        enddo
       enddo
 
cccc   extract the alignment:
       i=nseq1
       j=nseq2

       amax=0.0
       icho=0
       minseq=nseq1
       if (minseq .lt. nseq2) then
         icho=1
       endif
      
       if (icho .eq. 1) then
        do jj=nseq2,1,-1
           if (val(nseq1,jj).gt.amax) then
             amax=val(nseq1,jj)
             j=jj
           endif
      enddo
       else
        do ii=nseq1,1,-1
          if (val(ii,nseq2).gt.amax) then
             amax=val(ii,nseq2)
             i=ii
          endif
        enddo
       endif

       do while((i.gt.0).and.(j.gt.0))
         if(dir(i,j))then
           invmap_i(j)=i
           i=i-1
           j=j-1
         else
           h=val(i-1,j)
           v=val(i,j-1)

           if(dir(i-1,j))h=h+gap_open
           if(dir(i,j-1))v=v+gap_open

           if(v.ge.h) then
             j=j-1
           else
             i=i-1
           endif
         endif
       enddo
 
       j=nseq2
       do while(invmap_i(j).lt.1)
          j=j-1
       enddo
       dpout=val(invmap_i(j),j)
 
      return
      end 


**************************************************************
**  a kind of mixing of two previous alignments.
**  Takes the alignment and retuns the alignment!!
**************************************************************

      subroutine get_initial3 (d0,invmapi,score,scorer)
      parameter (maxres=1000)             ! no. of residues  
      dimension invmap_i(maxres),invmapi(maxres)
      dimension score(maxres,maxres),scorer(maxres,maxres)

      common /coord/ xa(3,maxres,0:1)
      common /secstr/ isec(maxres),jsec(maxres)
      common /length/ nseq1,nseq2

!$omp threadprivate(/coord/)
!$omp threadprivate(/secstr/)
!$omp threadprivate(/length/)

cccc  score matrix 
      do j=1,nseq2
         invmap_i(j)=invmapi(j)
      enddo

      call get_score1 (d0,invmap_i,score)          !get score(i,j) using RMSD martix

      do i=1,nseq1
         do j=1,nseq2
            if(isec(i).eq.jsec(j))then
               scorer(i,j)=0.5+score(i,j)
            else
               scorer(i,j)=score(i,j)
            endif
         enddo
      enddo

      return
      end

******************************************************************
**  with invmap_i(i) calculate score(i,j) using RMSD rotation
******************************************************************
      subroutine get_score1 (dx,invmapi,score)
      parameter (maxres=1000)             ! no. of residues  
      dimension invmapi(maxres)
      dimension score(maxres,maxres)

      common /coord/ xa(3,maxres,0:1)
      common /length/ nseq1,nseq2

!$omp threadprivate(/coord/)
!$omp threadprivate(/length/)

cccc  RMSD:
      double precision r_1(3,maxres),r_2(3,maxres),w(maxres)
      double precision u(3,3),t(3),rms !armsd is real
      data w /maxres*1.0/

cccc  calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=1,nseq2
         i=invmapi(j)            !j aligned to i
         if(i.gt.0)then
            n_al=n_al+1
            r_1(1,n_al)=xa(1,i,0) ! for rotation matrix
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
            r_2(1,n_al)=xa(1,j,1)
            r_2(2,n_al)=xa(2,j,1)
            r_2(3,n_al)=xa(3,j,1)
         endif
      enddo

cccc  calculate score matrix score(i,j)------------------>
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2
      
      d0_min=0.5
      d01=dx+1.5
      if(d01.lt.d0_min)d01=d0_min
      d02=d01**2

      do i=1,nseq1
         xx=real(t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+
     &           u(1,3)*xa(3,i,0))
         yy=real(t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+
     &           u(2,3)*xa(3,i,0))
         zz=real(t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+
     &           u(3,3)*xa(3,i,0))
         do j=1,nseq2
            dd=(xx-xa(1,j,1))**2+(yy-xa(2,j,1))**2+
     &         (zz-xa(3,j,1))**2
            score(i,j)=1/(1+(dd/d02))
         enddo
      enddo

      return
      end

******************************************************************
**  with invmap_i(i) calculate TM-score and martix score(i,j) for rotation
**  isearch takes 3 values for 3 different kind of TM-score calculation
**  isearch=1, searching with long jumps of 40 residues in TMsearch routine
******************************************************************

      subroutine get_score (d0,d0_s,invmap_i,score,TM)
      parameter (maxres=1000)             ! no. of residues  
      dimension invmap_i(maxres)
      dimension score(maxres,maxres)
      dimension xtm1(maxres),ytm1(maxres),ztm1(maxres)
      dimension xtm2(maxres),ytm2(maxres),ztm2(maxres)

      common /coord/ xa(3,maxres,0:1)
      common /n1n2/ n1(maxres),n2(maxres)
      common /length/ nseq1,nseq2

!$omp threadprivate(/coord/)
!$omp threadprivate(/n1n2/)
!$omp threadprivate(/length/)

cccc  RMSD:
      double precision r_1(3,maxres),r_2(3,maxres),w(maxres)
      double precision u(3,3),t(3),rms !armsd is real
      data w /maxres*1.0/

cccc  calculate RMSD between aligned structures and rotate the structures -->
      anseq=min(nseq1,nseq2)    ! for normalization of TM-score
      d02=d0**2

      n_al=0
      do j=1,nseq2
         i=invmap_i(j)            !j aligned to i
         if(i.gt.0)then
            n_al=n_al+1

            xtm1(n_al)=xa(1,i,0) !for TM-score
            ytm1(n_al)=xa(2,i,0)
            ztm1(n_al)=xa(3,i,0)
            xtm2(n_al)=xa(1,j,1)
            ytm2(n_al)=xa(2,j,1)
            ztm2(n_al)=xa(3,j,1)

ccc   for rotation matrix:
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
         endif
      enddo
      d0_input=d0
      isearch=1     ! VERY IMP. determines way TM-score is calculate and reported !

      call TMsearch (d0_input,d0_s,n_al,xtm1,ytm1,ztm1,n1,
     &     n_al,xtm2,ytm2,ztm2,n2,TM,Rcomm,Lcomm,isearch) !simplified search engine 
        
       TM=(TM/anseq)

cccc output score matrix score(i,j)
cccc protein 1 or r_1 is need to rotated since in previous the rotation is for aligned subset!!
        do i=1,n_al
            r_2(1,i)=xtm1(i)
            r_2(2,i)=ytm1(i)
            r_2(3,i)=ztm1(i)
        enddo

        call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2

          do i=1,nseq1
            xx=real(t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+
     &              u(1,3)*xa(3,i,0))
            yy=real(t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+
     &              u(2,3)*xa(3,i,0))
            zz=real(t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+
     &              u(3,3)*xa(3,i,0))
            do j=1,nseq2
                dd=(xx-xa(1,j,1))**2+(yy-xa(2,j,1))**2+
     &             (zz-xa(3,j,1))**2
                score(i,j)=1.0/(1+dd/d02)
            enddo   
          enddo

       return
       end

******************************************************************
** Make iteration for a given alignment with two gap opening's
**
******************************************************************

       subroutine make_iter (d0,d0_s,score,niter,invmapi,aTM,invmapr)
       parameter (maxres=1000)
       dimension gap(10),score(maxres,maxres),invmapi(maxres)
       dimension invmap_i(maxres),invmap(maxres),invmapr(maxres)

       common /length/ nseq1,nseq2

!$omp threadprivate(/length/)

       ngap=2
       gap(1)=-0.6
       gap(2)=0.0

cccc copy the alignment
       do j=1,nseq2
        invmap(j)=invmapi(j)
       enddo

        do 111 i_gap=1,ngap
        gap_open=gap(i_gap)

            do 222 id=1,niter
                call dp(score,gap_open,invmap_i,dpout)  ! returns invmap_i
                ilc_count=0
                do ik=1,nseq2
                    if (invmap_i(ik).gt.0)ilc_count=ilc_count+1
                enddo
                TM_old=0
                if (ilc_count.le.5) goto 33
                call get_score (d0,d0_s,invmap_i,score,TM)   ! refills the score matrix using invmap_i
                    if (TM.gt.aTM) then
                        aTM=TM
                        do j=1,nseq2
                            invmap(j)=invmap_i(j)
                        enddo
                    endif
                    if (id.gt.1) then
                        diff=abs(TM-TM_old)
                        if (diff.lt.0.000001)goto 33
                    endif
                TM_old=TM   
  222       continue
   33     continue
  111   continue

        do j=1,nseq2
          invmapr(j)=invmap(j)
        enddo

       return
       end

cccc Subroutines for TM-align
cccc TMsearch, cal_tmscore

      subroutine TMsearch(dx,dxs,L1,x1,y1,z1,n1,L2,x2,y2,z2,n2,
     &     TM,Rcomm,Lcomm,isearch)

      parameter (maxres=1000)
      dimension x1(maxres),y1(maxres),z1(maxres),n1(maxres)
      dimension x2(maxres),y2(maxres),z2(maxres),n2(maxres)
      dimension xa(maxres),ya(maxres),za(maxres)
      dimension k_ali(maxres),k_ali0(maxres)
      dimension L_ini(100)
      double precision score,score_max
      dimension iL0(maxres),i_ali(maxres)

      common /stru/xt(maxres),yt(maxres),zt(maxres),
     &             xb(maxres),yb(maxres),zb(maxres) !for cal_tmscore routine
      common /align/ iA(maxres),iB(maxres)  !shared in next subroutine
      common /length/ nseq1,nseq2
      
!$omp threadprivate(/stru/)
!$omp threadprivate(/align/)
!$omp threadprivate(/length/)
      
      logical nstat

ccc   RMSD:
      double precision r_1(3,maxres),r_2(3,maxres),w(maxres)
      double precision u(3,3),t(3),rms !armsd is real
      data w /maxres*1.0/
      
      ka0=0
      ijump=0
      xdum1=n1(1)
      xdum1=n2(1)
      

cccc  convert input data. Special case here of L1.eq.L2
      nstat=.false. ! used in cal_tmscore routine for decision 
                    ! for including residues for TM-score calculation
      nseqA=L1
      nseqB=L2
      do i=1,nseqA
         xa(i)=x1(i)
         ya(i)=y1(i)
         za(i)=z1(i)
         xb(i)=x2(i)
         yb(i)=y2(i)
         zb(i)=z2(i)
         iA(i)=i
         iB(i)=i
      enddo
      n_ali=L1                  !number of aligned residues
      Lcomm=L1

cccc d0 parameters 
      d0_min=0.5
      d0=dx
      if(d0.lt.d0_min)d0=d0_min
       d0_search=dxs

cccc  iterative parameters ----->
      n_it=20                   !maximum number of iterations
      n_init_max=6              !maximum number of L_init
      n_init=0                  ! number of fragments L/2**(n)
      L_ini_min=4

      if(n_ali.lt.4)L_ini_min=n_ali

      do i=1,n_init_max-1
         n_init=n_init+1
         L_ini(n_init)=n_ali/2**(n_init-1)
         if(L_ini(n_init).le.L_ini_min)then
            L_ini(n_init)=L_ini_min
            goto 402
         endif
      enddo
      n_init=n_init+1
      L_ini(n_init)=L_ini_min
 402  continue

cccc  find the maximum score starting from local structures superposition

      score_max=-1              !TM-score
      if (isearch.eq.1) then
        ijump=40
        if (min(nseq1,nseq2).le.80)ijump=30
      endif
      if (isearch.ge.2) ijump=1
      if (isearch.eq.1.or.isearch.eq.2) nstat=.true.

      do 333 i_init=1,n_init
        L_init=L_ini(i_init)    ! Length to be used for scanning L
        iL_max=n_ali-L_init+1   ! maximum number of scans (Nali-L+1), shifting by one residue
c      if (isearch.gt.2) write(*,*)iL_max,score_max,d0_search,d0

        if (isearch.eq.1) then
         k=0
            do i=1,iL_max,ijump     ! this is the for the real quick result
                k=k+1               ! when isearch.eq.1
                iL0(k)=i            ! storing the starting residue for scan
            enddo

            if(iL0(k).lt.iL_max)then
                k=k+1
                iL0(k)=iL_max   ! fixing the last run, quite arbitarary!!
            endif
            n_shift=k           ! number of shifts
        else
           n_shift=iL_max
        endif    

              do 300 i_shift=1,n_shift
                    if (isearch.eq.1) then
                     iL=iL0(i_shift)
                    else
                     iL=i_shift
                    endif

                    LL=0
                    ka=0
                      do i=1,L_init
                          k=iL+i-1          ![1,n_ali] common aligned
                          r_1(1,i)=xa(iA(k))
                          r_1(2,i)=ya(iA(k))
                          r_1(3,i)=za(iA(k))
                          r_2(1,i)=xb(iB(k))
                          r_2(2,i)=yb(iB(k))
                          r_2(3,i)=zb(iB(k))
                          LL=LL+1
                          ka=ka+1
                          k_ali(ka)=k
                      enddo

                    call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2

                    if(i_init.eq.1)then  !global superposition
                     armsd=real(dsqrt(rms/LL))
                     Rcomm=armsd
                    endif

                    do j=1,nseqA
                      xt(j)=real(t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+
     &                           u(1,3)*za(j))
                      yt(j)=real(t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+
     &                           u(2,3)*za(j))
                      zt(j)=real(t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+
     &                           u(3,3)*za(j))
                    enddo

                    d=d0_search-1                                 ! reduce the search rms first
                    call cal_tmscore(d,d0,n_ali,nstat,score,i_ali,
     &                               n_cut) !init, get scores, n_cut+i_ali(i) for iteration
                    if(score_max.lt.score)then
                        score_max=score
                        ka0=ka
                        do i=1,ka0
                            k_ali0(i)=k_ali(i)  ! alignment stored
                        enddo
                    endif
cccc  iteration for extending
                d=d0_search+1
                do 301 it=1,n_it
                    LL=0
                    ka=0
                    do i=1,n_cut
                        m=i_ali(i)          ![1,n_ali]
                        r_1(1,i)=xa(iA(m))
                        r_1(2,i)=ya(iA(m))
                        r_1(3,i)=za(iA(m))
                        r_2(1,i)=xb(iB(m))
                        r_2(2,i)=yb(iB(m))
                        r_2(3,i)=zb(iB(m))
                        ka=ka+1
                        k_ali(ka)=m
                        LL=LL+1
                    enddo
                    call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2

                    do j=1,nseqA
                     xt(j)=real(t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+
     &                          u(1,3)*za(j))
                     yt(j)=real(t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+
     &                          u(2,3)*za(j))
                     zt(j)=real(t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+
     &                          u(3,3)*za(j))
                    enddo

                    call cal_tmscore(d,d0,n_ali,nstat,score,i_ali,
     &                               n_cut) !get scores, n_cut+i_ali(i) for iteration

                    if(score_max.lt.score)then
                        score_max=score
                        ka0=ka
                        do i=1,ka
                            k_ali0(i)=k_ali(i)
                        enddo
                    endif

                    if(it.eq.n_it)goto 302

                    if(n_cut.eq.ka)then
                        neq=0
                        do i=1,n_cut
                            if(i_ali(i).eq.k_ali(i))neq=neq+1
                        enddo
                        if(n_cut.eq.neq)goto 302
                    endif
 301            continue             !for iteration
 302           continue
 300          continue                !for shift
 333  continue                  !for initial length, L_ali/M

cccc  return the final rotation with rotated (xtm1 or x1)
      LL=0
      do i=1,ka0
         m=k_ali0(i)            !record of the best alignment
         r_1(1,i)=xa(iA(m))
         r_1(2,i)=ya(iA(m))
         r_1(3,i)=za(iA(m))
         r_2(1,i)=xb(iB(m))
         r_2(2,i)=yb(iB(m))
         r_2(3,i)=zb(iB(m))
         LL=LL+1
      enddo

      call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2

      do j=1,nseqA
         x1(j)=real(t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j))
         y1(j)=real(t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j))
         z1(j)=real(t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j))
      enddo
      
      TM=real(score_max)

      return
      end

cccc rotuine to calculate TM-score given the superposed coordinates
cccc IMP: superposed coordinates are passed in common
cccc d is for selecting the residues with certain cutoff and d0 is for
cccc TM-score calculation. nstat is type of calculation d or d8

      subroutine cal_tmscore (d,d0,n_ali,nstat,score,i_ali,n_cut)
      parameter (maxres=1000)
      double precision score
      dimension i_ali(maxres)

      common/stru/xt(maxres),yt(maxres),zt(maxres),
     &            xb(maxres),yb(maxres),zb(maxres)
      common /align/ iA(maxres),iB(maxres)
      common /length/ nseq1,nseq2
      common /dinfo/ d8
      
!$omp threadprivate(/stru/)
!$omp threadprivate(/align/)
!$omp threadprivate(/length/)
!$omp threadprivate(/dinfo/)
      
      logical nstat

      d_cp=d        ! for preventing superposition of 2 residues and less
      nseqB=n_ali   ! normalized for TMscore

  116 continue
        n_cut=0                   !number of residue-pairs dis<d, for iteration
        score_sum=0               !TMscore
        do k=1,n_ali
            i=iA(k)                ![1,nseqA] reoder number of structureA
            j=iB(k)                ![1,nseqB]
            dis=sqrt((xt(i)-xb(j))**2+(yt(i)-yb(j))**2+
     &               (zt(i)-zb(j))**2)
                if(dis.lt.d_cp)then
                    n_cut=n_cut+1
                    i_ali(n_cut)=k      ![1,n_ali], mark the residue-pairs in dis<d
                endif

                if(nstat) then
                    if(dis.le.d8)then
                        score_sum=score_sum+1.0/(1+(dis/d0)**2)
                    endif
                else
                        score_sum=score_sum+1.0/(1+(dis/d0)**2)
                endif
        enddo

        if(n_cut.gt.0.and.n_cut.lt.3) then
            d_cp=d_cp+0.5
            goto 116
        endif
        score=score_sum !TM-score

      return
      end

cccccccccccccccc Calculate sum of (r_d-r_m)^2 cccccccccccccccccccccccccc
c  w    - w(m) is weight for atom pair  c m           (given)
c  x    - x(i,m) are coordinates of atom c m in set x       (given)
c  y    - y(i,m) are coordinates of atom c m in set y       (given)
c  n    - n is number of atom pairs                         (given)
c  mode  - 0:calculate rms only                             (given)
c          1:calculate rms,u,t                              (takes longer)
c  rms   - sum of w*(ux+t-y)**2 over all atom pairs         (result)
c  u    - u(i,j) is   rotation  matrix for best superposition  (result)
c  t    - t(i)   is translation vector for best superposition  (result)
c  ier  - 0: a unique optimal superposition has been determined(result)
c       -1: superposition is not unique but optimal
c       -2: no result obtained because of negative weights w
c           or all weights equal to zero.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine u3b(w, x, y, n, mode, rms, u, t, ier)
      integer ip(9), ip2312(4), i, j, k, l, m1, m, ier, n, mode
      double precision w(n), x(3, n), y(3, n), u(3, 3), t(3), rms, sigma
      double precision r(3, 3), xc(3), yc(3), wc, a(3, 3), b(3, 3), e0, 
     &e(3), e1, e2, e3, d, spur, det, cof, h, g, cth, sth, sqrth, p, tol
     &, rr(6), rr1, rr2, rr3, rr4, rr5, rr6, ss(6), ss1, ss2, ss3, ss4, 
     &ss5, ss6, zero, one, two, three, sqrt3
      equivalence (rr(1), rr1), (rr(2), rr2), (rr(3), rr3), (rr(4), rr4)
     &, (rr(5), rr5), (rr(6), rr6), (ss(1), ss1), (ss(2), ss2), (ss(3), 
     &ss3), (ss(4), ss4), (ss(5), ss5), (ss(6), ss6), (e(1), e1), (e(2)
     &, e2), (e(3), e3)
      data sqrt3 / 1.73205080756888d+00 /
      data tol / 1.0d-2 /
      data zero / 0.0d+00 /
      data one / 1.0d+00 /
      data two / 2.0d+00 /
      data three / 3.0d+00 /
      data ip / 1, 2, 4, 2, 3, 5, 4, 5, 6 /
      data ip2312 / 2, 3, 1, 2 /
c 156 "rms.for"
      wc = zero
      rms = 0.0
      e0 = zero
      do 1 i = 1, 3
      xc(i) = zero
      yc(i) = zero
      t(i) = 0.0
      do 1 j = 1, 3
      d = zero
      if (i .eq. j) d = one
      u(i,j) = d
      a(i,j) = d
    1 r(i,j) = zero
      ier = -1
c**** DETERMINE CENTROIDS OF BOTH VECTOR SETS X AND Y
c 170 "rms.for"
      if (n .lt. 1) return 
c 172 "rms.for"
      ier = -2
      do 2 m = 1, n
      if (w(m) .lt. 0.0) return 
      wc = wc + w(m)
      do 2 i = 1, 3
      xc(i) = xc(i) + (w(m) * x(i,m))
    2 yc(i) = yc(i) + (w(m) * y(i,m))
      if (wc .le. zero) return 
      do 3 i = 1, 3
      xc(i) = xc(i) / wc
c**** DETERMINE CORRELATION MATRIX R BETWEEN VECTOR SETS Y AND X
c 182 "rms.for"
    3 yc(i) = yc(i) / wc
c 184 "rms.for"
      do 4 m = 1, n
      do 4 i = 1, 3
      e0 = e0 + (w(m) * (((x(i,m) - xc(i)) ** 2) + ((y(i,m) - yc(i)) ** 
     &2)))
c 187 "rms.for"
      d = w(m) * (y(i,m) - yc(i))
      do 4 j = 1, 3
c**** CALCULATE DETERMINANT OF R(I,J)
c 189 "rms.for"
    4 r(i,j) = r(i,j) + (d * (x(j,m) - xc(j)))
c 191 "rms.for"
      det = ((r(1,1) * ((r(2,2) * r(3,3)) - (r(2,3) * r(3,2)))) - (r(1,2
     &) * ((r(2,1) * r(3,3)) - (r(2,3) * r(3,1))))) + (r(1,3) * ((r(2,1)
     & * r(3,2)) - (r(2,2) * r(3,1))))
c**** FORM UPPER TRIANGLE OF TRANSPOSED(R)*R
c 194 "rms.for"
      sigma = det
c 196 "rms.for"
      m = 0
      do 5 j = 1, 3
      do 5 i = 1, j
      m = m + 1
c***************** EIGENVALUES *****************************************
c**** FORM CHARACTERISTIC CUBIC  X**3-3*SPUR*X**2+3*COF*X-DET=0
c 200 "rms.for"
    5 rr(m) = ((r(1,i) * r(1,j)) + (r(2,i) * r(2,j))) + (r(3,i) * r(3,j)
     &)
c 203 "rms.for"
      spur = ((rr1 + rr3) + rr6) / three
      cof = ((((((rr3 * rr6) - (rr5 * rr5)) + (rr1 * rr6)) - (rr4 * rr4)
     &) + (rr1 * rr3)) - (rr2 * rr2)) / three
c 205 "rms.for"
      det = det * det
      do 6 i = 1, 3
    6 e(i) = spur
c**** REDUCE CUBIC TO STANDARD FORM Y**3-3HY+2G=0 BY PUTTING X=Y+SPUR
c 208 "rms.for"
      if (spur .le. zero) goto 40
c 210 "rms.for"
      d = spur * spur
      h = d - cof
c**** SOLVE CUBIC. ROOTS ARE E1,E2,E3 IN DECREASING ORDER
c 212 "rms.for"
      g = (((spur * cof) - det) / two) - (spur * h)
c 214 "rms.for"
      if (h .le. zero) goto 8
      sqrth = dsqrt(h)
      d = ((h * h) * h) - (g * g)
      if (d .lt. zero) d = zero
      d = datan2(dsqrt(d),- g) / three
      cth = sqrth * dcos(d)
      sth = (sqrth * sqrt3) * dsin(d)
      e1 = (spur + cth) + cth
      e2 = (spur - cth) + sth
      e3 = (spur - cth) - sth
c.....HANDLE SPECIAL CASE OF 3 IDENTICAL ROOTS
c 224 "rms.for"
      if (mode) 10, 50, 10
c**************** EIGENVECTORS *****************************************
c 226 "rms.for"
    8 if (mode) 30, 50, 30
c 228 "rms.for"
   10 do 15 l = 1, 3, 2
      d = e(l)
      ss1 = ((d - rr3) * (d - rr6)) - (rr5 * rr5)
      ss2 = ((d - rr6) * rr2) + (rr4 * rr5)
      ss3 = ((d - rr1) * (d - rr6)) - (rr4 * rr4)
      ss4 = ((d - rr3) * rr4) + (rr2 * rr5)
      ss5 = ((d - rr1) * rr5) + (rr2 * rr4)
      ss6 = ((d - rr1) * (d - rr3)) - (rr2 * rr2)
      j = 1
      if (dabs(ss1) .ge. dabs(ss3)) goto 12
      j = 2
      if (dabs(ss3) .ge. dabs(ss6)) goto 13
   11 j = 3
      goto 13
   12 if (dabs(ss1) .lt. dabs(ss6)) goto 11
   13 d = zero
      j = 3 * (j - 1)
      do 14 i = 1, 3
      k = ip(i + j)
      a(i,l) = ss(k)
   14 d = d + (ss(k) * ss(k))
      if (d .gt. zero) d = one / dsqrt(d)
      do 15 i = 1, 3
   15 a(i,l) = a(i,l) * d
      d = ((a(1,1) * a(1,3)) + (a(2,1) * a(2,3))) + (a(3,1) * a(3,3))
      m1 = 3
      m = 1
      if ((e1 - e2) .gt. (e2 - e3)) goto 16
      m1 = 1
      m = 3
   16 p = zero
      do 17 i = 1, 3
      a(i,m1) = a(i,m1) - (d * a(i,m))
   17 p = p + (a(i,m1) ** 2)
      if (p .le. tol) goto 19
      p = one / dsqrt(p)
      do 18 i = 1, 3
   18 a(i,m1) = a(i,m1) * p
      goto 21
   19 p = one
      do 20 i = 1, 3
      if (p .lt. dabs(a(i,m))) goto 20
      p = dabs(a(i,m))
      j = i
   20 continue
      k = ip2312(j)
      l = ip2312(j + 1)
      p = dsqrt((a(k,m) ** 2) + (a(l,m) ** 2))
      if (p .le. tol) goto 40
      a(j,m1) = zero
      a(k,m1) = - (a(l,m) / p)
      a(l,m1) = a(k,m) / p
   21 a(1,2) = (a(2,3) * a(3,1)) - (a(2,1) * a(3,3))
      a(2,2) = (a(3,3) * a(1,1)) - (a(3,1) * a(1,3))
c****************** ROTATION MATRIX ************************************
c 282 "rms.for"
      a(3,2) = (a(1,3) * a(2,1)) - (a(1,1) * a(2,3))
c 284 "rms.for"
   30 do 32 l = 1, 2
      d = zero
      do 31 i = 1, 3
      b(i,l) = ((r(i,1) * a(1,l)) + (r(i,2) * a(2,l))) + (r(i,3) * a(3,l
     &))
c 288 "rms.for"
   31 d = d + (b(i,l) ** 2)
      if (d .gt. zero) d = one / dsqrt(d)
      do 32 i = 1, 3
   32 b(i,l) = b(i,l) * d
      d = ((b(1,1) * b(1,2)) + (b(2,1) * b(2,2))) + (b(3,1) * b(3,2))
      p = zero
      do 33 i = 1, 3
      b(i,2) = b(i,2) - (d * b(i,1))
   33 p = p + (b(i,2) ** 2)
      if (p .le. tol) goto 35
      p = one / dsqrt(p)
      do 34 i = 1, 3
   34 b(i,2) = b(i,2) * p
      goto 37
   35 p = one
      do 36 i = 1, 3
      if (p .lt. dabs(b(i,1))) goto 36
      p = dabs(b(i,1))
      j = i
   36 continue
      k = ip2312(j)
      l = ip2312(j + 1)
      p = dsqrt((b(k,1) ** 2) + (b(l,1) ** 2))
      if (p .le. tol) goto 40
      b(j,2) = zero
      b(k,2) = - (b(l,1) / p)
      b(l,2) = b(k,1) / p
   37 b(1,3) = (b(2,1) * b(3,2)) - (b(2,2) * b(3,1))
      b(2,3) = (b(3,1) * b(1,2)) - (b(3,2) * b(1,1))
      b(3,3) = (b(1,1) * b(2,2)) - (b(1,2) * b(2,1))
      do 39 i = 1, 3
      do 39 j = 1, 3
c****************** TRANSLATION VECTOR *********************************
c 320 "rms.for"
   39 u(i,j) = ((b(i,1) * a(j,1)) + (b(i,2) * a(j,2))) + (b(i,3) * a(j,3
     &))
   40 do 41 i = 1, 3
c********************** RMS ERROR **************************************
c 323 "rms.for"
   41 t(i) = ((yc(i) - (u(i,1) * xc(1))) - (u(i,2) * xc(2))) - (u(i,3)
     & * xc(3))
   50 do 51 i = 1, 3
      if (e(i) .lt. zero) e(i) = zero
   51 e(i) = dsqrt(e(i))
      ier = 0
      if (e2 .le. (e1 * 1.0d-05)) ier = -1
      d = e3
      if (sigma .ge. 0.0) goto 52
      d = - d
      if ((e2 - e3) .le. (e1 * 1.0d-05)) ier = -1
   52 d = (d + e2) + e1
      rms = (e0 - d) - d
      if (rms .lt. 0.0) rms = 0.0
      return 
c.....END U3B...........................................................
c----------------------------------------------------------
c                       THE END
c----------------------------------------------------------
c 338 "rms.for"
      end

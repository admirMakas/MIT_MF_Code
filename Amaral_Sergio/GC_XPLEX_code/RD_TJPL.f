C $Id: RD_TJPL.f,v 1.1.1.1 2009/06/09 21:51:51 daven Exp $
      SUBROUTINE RD_TJPL(NJ1,NAMFIL)
C-----------------------------------------------------------------------
c  Read in wavelength bins, solar fluxes, Rayleigh parameters, temperature-
c  dependent cross sections and Rayleigh/aerosol scattering phase functions
c  with temperature dependences. Current data originates from JPL'97
C-----------------------------------------------------------------------
c
c     NAMFIL   Name of spectral data file (jv_spec.dat)
c     NJ1      Channel number for reading data file
c     NJVAL    Number of species to calculate J-values for
c     NWWW     Number of wavelength bins, from NW1:NW2
c     WBIN     Boundaries of wavelength bins
c     WL       Centres of wavelength bins - 'effective wavelength'
c     FL       Solar flux incident on top of atmosphere (cm-2.s-1)
c     QRAYL    Rayleigh parameters (effective cross-section) (cm2)
c     QBC      Black Carbon absorption extinct. (specific cross-sect.) (m2/g)
c     QO2      O2 cross-sections
c     QO3      O3 cross-sections
c     Q1D      O3 => O(1D) quantum yield
c     TQQ      Temperature for supplied cross sections
c     QQQ      Supplied cross sections in each wavelength bin (cm2)
c     NAA      Number of categories for scattering phase functions
c     QAA      Aerosol scattering phase functions
c     NK       Number of wavelengths at which functions supplied (set as 4)
c     WAA      Wavelengths for the NK supplied phase functions
c     PAA      Phase function: first 8 terms of expansion
c     RAA      Effective radius associated with aerosol type
c     SSA      Single scattering albedo
c
c     npdep    Number of pressure dependencies
c     zpdep    Pressure dependencies by wavelength bin
c     jpdep    Index of cross sections requiring pressure dependence
c     lpdep    Label for pressure dependence
c
c  NOTES:
c  (1 ) Updated to include new pressure-dependancy function for GLYX and MGLY. 
c        (tmf, 1/7/09)
c  (2 ) Added a pressure-dependancy function selector 'pdepf'. 
c        (tmf, ccc, 1/7/09)
C-----------------------------------------------------------------------
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "cmn_fj.h"
#     include "jv_cmn.h"

      integer i, j, k, iw, nk, nqqq, nwww, nj1
      character*7  lpdep(7)
      character*11 NAMFIL
      !real*8 D_WBIN(NW+1),D_WL(NW),D_FL(NW),D_QRAYL(NW),D_QBC(NW),
      !& D_TQQ(3,NS),D_QO2(NW,3),D_QO3(NW,3),D_Q1D(NW,3),D_QQQ(NW,2,NS-3),
      !& d_zpdep(NW,7),d_mglypdep(NW,3),D_WAA(4,NP),D_QAA(4,NP),
      !& D_RAA(4,NP),D_SSA(4,NP),D_PAA(8,4,NP)
      do J=1,NS
        do K=1,3
          D_TQQ(K,J) = 0.d0
          TQQ(K,J) = 0.d0
        enddo
      enddo
C-------------spectral data---------------------------------------------
      open(NJ1, FILE=NAMFIL)
      read(NJ1,'(A)') TITLE0
      write(6,'(1X,A)') TITLE0
      read(NJ1,'(10X,14I5)') NJVAL,NWWW,NW1,NW2
      if(NJVAL.gt.NS) then
        write(6,300) NJVAL,NS
        stop
      endif
C------------NQQQ = no. additional J-values from X-sects (O2,O3P,O3D+NQQQ)
C- NQQQ is changed to NJVAL-1 because there are 2 dummy species at the end
C used for acetone pressure dependency only. (ccc, 4/20/09)
C- prior to 4/20/09
C      NQQQ = NJVAL-3
      NQQQ = NJVAL-1
      read(NJ1,102) (D_WBIN(IW),IW=1,NWWW)
      read(NJ1,102) (D_WBIN(IW+1),IW=1,NWWW)
      read(NJ1,102) (D_WL(IW),IW=1,NWWW)
      read(NJ1,102) (D_FL(IW),IW=1,NWWW)
      read(NJ1,102) (D_QRAYL(IW),IW=1,NWWW)
      read(NJ1,102) (D_QBC(IW),IW=1,NWWW)   !  From Loiusse et al. [JGR, 1996]
      WBIN(:) = (D_WBIN(:))
      WL(:) = (D_WL(:))
      FL(:) = (D_FL(:))
      QRAYL(:) = (D_QRAYL(:))
      QBC(:) = (D_QBC(:))
c
C---Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
      do K=1,3
        read(NJ1,103) TITLEJ(K,1),D_TQQ(K,1), (D_QO2(IW,K),IW=1,NWWW)
        TQQ(K,1) = (D_TQQ(K,1))
        QO2(:,K) = (D_QO2(:,K))
      enddo
      do K=1,3
        read(NJ1,103) TITLEJ(K,2),D_TQQ(K,2), (D_QO3(IW,K),IW=1,NWWW)
        TQQ(K,2) = (D_TQQ(K,2))
        QO3(:,K) = (D_QO3(:,K))
      enddo
      do K=1,3
        read(NJ1,103) TITLEJ(K,3),D_TQQ(K,3), (D_Q1D(IW,K),IW=1,NWWW)
        TQQ(K,3) = (D_TQQ(K,3))
        Q1D(:,K) = (D_Q1D(:,K))
      enddo
      do K=1,3
        write(6,200) titlej(1,k),(tqq(i,k),i=1,3)
      enddo
c
C---Read remaining species:  X-sections at 2 T's
      do J=1,NQQQ
      read(NJ1,103) TITLEJ(1,J+3),D_TQQ(1,J+3),(D_QQQ(IW,1,J),IW=1,NWWW)
      read(NJ1,103) TITLEJ(2,J+3),D_TQQ(2,J+3),(D_QQQ(IW,2,J),IW=1,NWWW)
      TQQ(1,J+3) = (D_TQQ(1,J+3))
      QQQ(:,1,J) = (D_QQQ(:,1,J))
      TQQ(2,J+3) = (D_TQQ(2,J+3))
      QQQ(:,2,J) = (D_QQQ(:,2,J))
        write(6,200) titlej(1,j+3),(tqq(i,j+3),i=1,2)
      enddo
      read(NJ1,'(A)') TITLE0
c
c---Pressure dependencies
      read(NJ1,104) npdep
      do k=1,npdep
         read(NJ1,105) lpdep(k), pdepf(k), (d_zpdep(iw,k),iw=1,nwww)
         zpdep(:,k) = (d_zpdep(:,k))
         write(6,201)  lpdep(k), pdepf(k), (zpdep(iw,k),iw=1,nwww)

         !--------------------------------------
         ! Special treatment for MGLY pressure dependency
         ! (tmf, 11/16/06)
         !--------------------------------------
         if ( pdepf(k) .eq. 4 ) then           
            ! pass zpdep to mglypdep
            mglypdep(:,1) = zpdep(:,k)
            read(NJ1,105) lpdep(k),pdepf(k),(d_mglypdep(iw,2),iw=1,nwww)
            read(NJ1,105) lpdep(k),pdepf(k),(d_mglypdep(iw,3),iw=1,nwww)
            mglypdep(:,2) = (d_mglypdep(:,2))
            mglypdep(:,3) = (d_mglypdep(:,3))
         endif
      enddo
      read(NJ1,'(A)') TITLE0

c
c---Zero index arrays
      do j=1,jppj
        jind(j)=0
      enddo
      do j=1,NJVAL
        jpdep(j)=0
      enddo
c
C---Set mapping index
      do j=1,NJVAL
        do k=1,jppj
          if (jlabel(k).eq.titlej(1,j)) jind(k)=j
        enddo
        do k=1,npdep
          if (lpdep(k).eq.titlej(1,j)) jpdep(j)=k
        enddo
      enddo
      do k=1,jppj
        if(jfacta(k).eq.0.d0)
     &             write(6,*) 'Not using photolysis reaction ',k
        if(jind(k).eq.0) then
          if(jfacta(k).eq.0.d0) then
            jind(k)=1
          else
            write(6,*) 'Which J-rate for photolysis reaction ',k,' ?'
            stop
          endif
        endif
      enddo
c
C---Read aerosol phase functions:
      read(NJ1,'(A10,I5,/)') TITLE0,NAA
      NK=4        ! Fix number of wavelengths at 4
      do j=1,NAA
        read(NJ1,110) TITLEA(j)
        do k=1,NK
          read(NJ1,*) D_WAA(k,j),D_QAA(k,j),D_RAA(k,j),D_SSA(k,j),
     &                                      (D_PAA(i,k,j),i=1,8)
          WAA(k,j) = (D_WAA(k,j))
          QAA(k,j) = (D_QAA(k,j))
          RAA(k,j) = (D_RAA(k,j))
          SSA(k,j) = (D_SSA(k,j))
          PAA(:,k,j) = (D_PAA(:,k,j))
        enddo
      enddo
c
      write(6,*) 'Aerosol phase functions & wavelengths'
      do J=1,NAA
        write(6,'(1x,A8,I2,A,18F8.1)')
     $                   TITLEA(J),J,'  wavel=',(WAA(K,J),K=1,NK)
        write(6,'(9x,I2,A,18F8.4)') J,'  Qext =',(QAA(K,J),K=1,NK)
      enddo
C--------
C Modify reading and writing formats 105 & 201 for pressure dependancy 
c (ccc, 1/7/09)
      !print*,'WBIN',WBIN
      !print*,'WL',WL
      !print*,'FL',FL
      !print*,'QRAYL',QRAYL
      !print*,'QBC',QBC
      !print*,'TQQ',TQQ
      !print*,'QO2',QO2
      !print*,'QO3',QO3
      !print*,'Q1D',Q1D
      !print*,'QQQ',QQQ
      !print*,'zpdep',zpdep
      !print*,'mglypdep',mglypdep
      !print*,'WAA',WAA
      !print*,'QAA',QAA
      !print*,'RAA',RAA
      !print*,'SSA',SSA
      !print*,'PAA',PAA
  101 FORMAT(8E10.3)
  102 FORMAT(10X,7E10.3)
  103 FORMAT(A7,F3.0,7E10.3)
c 103 FORMAT(A7,F3.0,7E10.3/(10X,7E10.3))
  104 FORMAT(13x,i2)
  105 FORMAT(A7,2x,I1,7E10.3)
  110 format(3x,a20)
  200 format(1x,' x-sect:',a10,3(3x,2f6.2))
  201 format(1x,' pr.dep:',a10,1x,I1,7(1pE10.3))
  300 format(' Number of x-sections supplied to Fast-J: ',i3,/,
     &       ' Maximum number allowed (NS) only set to: ',i3,
     &       ' - increase in jv_cmn.h')
      close(NJ1)
      return
      end

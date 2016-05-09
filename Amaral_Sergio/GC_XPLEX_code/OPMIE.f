C $Id: OPMIE.f,v 1.1.1.1 2009/06/09 21:51:54 daven Exp $
      SUBROUTINE OPMIE(KW,WAVEL,XQO2,XQO3,FMEAN)
C-----------------------------------------------------------------------
C  NEW Mie code for J's, only uses 8-term expansion, 4-Gauss pts
C  Currently allow up to NP aerosol phase functions (at all altitudes) to
C  be associated with optical depth AER(1:NC) = aerosol opt.depth @ 1000 nm
C  
C  Pick Mie-wavelength with phase function and Qext:
C
C  01 RAYLE = Rayleigh phase
C  02 ISOTR = isotropic
C  03 ABSRB = fully absorbing 'soot', wavelength indep.
C  04 S_Bkg = backgrnd stratospheric sulfate (n=1.46,log-norm:r=.09um/sigma=.6)
C  05 S_Vol = volcanic stratospheric sulfate (n=1.46,log-norm:r=.08um/sigma=.8)
C  06 W_H01 = water haze (H1/Deirm.) (n=1.335, gamma:  r-mode=0.1um /alpha=2)
C  07 W_H04 = water haze (H1/Deirm.) (n=1.335, gamma:  r-mode=0.4um /alpha=2)
C  08 W_C02 = water cloud (C1/Deirm.) (n=1.335, gamma:  r-mode=2.0um /alpha=6)
C  09 W_C04 = water cloud (C1/Deirm.) (n=1.335, gamma:  r-mode=4.0um /alpha=6)
C  10 W_C08 = water cloud (C1/Deirm.) (n=1.335, gamma:  r-mode=8.0um /alpha=6)
C  11 W_C13 = water cloud (C1/Deirm.) (n=1.335, gamma:  r-mode=13.3um /alpha=6)
C  12 W_L06 = water cloud (Lacis) (n=1.335, r-mode=5.5um / alpha=11/3)
C  13 Ice-H = hexagonal ice cloud (Mishchenko)
C  14 Ice-I = irregular ice cloud (Mishchenko)
C
C  Choice of aerosol index MIEDX is made in SET_AER; optical depths are
C  apportioned to the AER array in SET_PROF
C
C-----------------------------------------------------------------------
C  FUNCTION RAYLAY(WAVE)---RAYLEIGH CROSS-SECTION for wave > 170 nm
C       WSQI = 1.E6/(WAVE*WAVE)
C       REFRM1 = 1.0E-6*(64.328+29498.1/(146.-WSQI)+255.4/(41.-WSQI))
C       RAYLAY = 5.40E-21*(REFRM1*WSQI)**2
C-----------------------------------------------------------------------
c
c     DTAUX    Local optical depth of each CTM level
c     PIRAY    Contribution of Rayleigh scattering to extinction
c     PIAER    Contribution of Aerosol scattering to extinction
c     TTAU     Optical depth of air vertically above each point (to top of atm)
c     FTAU     Attenuation of solar beam
c     POMEGA   Scattering phase function
c     FMEAN    Mean actinic flux at desired levels
c
C-----------------------------------------------------------------------
      USE MYTYPE
      USE COMPLEXIFY
      USE ERROR_MOD
   
      IMPLICIT NONE

#     include "cmn_fj.h"
#     include "jv_cmn.h"
#     include "jv_mie.h"

      integer jndlev(lpar),jaddlv(nc),jaddto(nc+1)
      integer KW,km,i,j,k,l,ix,j1
      TYPE (XPLEX) QXMIE(MX),XLAER(MX),SSALB(MX)
      TYPE (XPLEX) xlo2,xlo3,xlray,xltau,zk,taudn,tauup,zk2
      TYPE (XPLEX) WAVEL,XQO2(NB),XQO3(NB),FMEAN(lpar),
     & POMEGAJ(2*M__,NC+1)
      TYPE(XPLEX) DTAUX(NB),PIRAY(NB),PIAER(MX,NB),TTAU(NC+1),FTAU(NC+1)
      TYPE (XPLEX) ftaulog,dttau,dpomega(2*M__)
      TYPE (XPLEX) ftaulog2,dttau2,dpomega2(2*M__)

      ! For KLUDGE to fix the # of added levels (phs, 7/1/08)
      INTEGER :: loc(1)
      
ci
C---Pick nearest Mie wavelength, no interpolation--------------
                              KM=1
      if( WAVEL .gt. 355.d0 ) KM=2
      if( WAVEL .gt. 500.d0 ) KM=3
C     if( WAVEL .gt. 800.d0 ) KM=4  !drop the 1000 nm wavelength
c
C---For Mie code scale extinction at 1000 nm to wavelength WAVEL (QXMIE)
      do I=1,MX
        QXMIE(I) = QAA(KM,MIEDX(I))/QAA(4,MIEDX(I))
        SSALB(I) = SSA(KM,MIEDX(I))
        !if (isnan(QXMIE(I)) .or. isnan(SSALB(I))) then
        !   print*,'QXMIE(I),SSALB(I)',QXMIE(I),SSALB(I)
        !   CALL GEOS_CHEM_STOP
        !endif
      enddo
c
C---Reinitialize arrays
      do j=1,nc+1
        ttau(j)=0.d0
        ftau(j)=0.d0
      enddo
c
C---Set up total optical depth over each CTM level, DTAUX
      J1 = NLBATM
      do J=J1,NB
        XLO3=DO3(J)*XQO3(J)
        XLO2=DM(J)*XQO2(J)*0.20948d0
        XLRAY=DM(J)*QRAYL(KW)
c  Zero absorption for testing purposes
c        call NOABS(XLO3,XLO2,XLRAY,AER(1,j),RFLECT)
        do I=1,MX
          XLAER(I)=AER(I,J)*QXMIE(I)
        enddo
c  Total optical depth from all elements
        DTAUX(J)=XLO3+XLO2+XLRAY
      
        do I=1,MX
          DTAUX(J)=DTAUX(J)+XLAER(I)
          
        enddo
c  Fractional extinction for Rayleigh scattering and each aerosol type
        PIRAY(J)=XLRAY/DTAUX(J)
      
        do I=1,MX
          PIAER(I,J)=SSALB(I)*XLAER(I)/DTAUX(J)
        enddo
      enddo
c
C---Define the scattering phase fn. with mix of Rayleigh(1) & Mie(MIEDX)
C   No. of quadrature pts fixed at 4 (M__), expansion of phase fn @ 8
      N = M__
      MFIT = 2*M__
      do j=j1,NB
        do i=1,MFIT
          pomegaj(i,j) = PIRAY(J)*PAA(I,KM,1)
               
          do k=1,MX
            pomegaj(i,j) = pomegaj(i,j) + PIAER(K,J)*PAA(I,KM,MIEDX(K))
          enddo
        enddo
      enddo
c
C---Calculate attenuated incident beam EXP(-TTAU/U0) and flux on surface
      do J=J1,NB
        if(AMF(J,J).gt.0.0D0) then
          XLTAU=0.0D0
          do I=1,NB
            XLTAU=XLTAU + DTAUX(I)*AMF(I,J)
          enddo
          if(XLTAU.gt.450.d0) then   ! for compilers with no underflow trapping
            FTAU(j)=0.d0
          else
            FTAU(J)=EXP(-XLTAU)
!            if (isnan(FTAU(J))) then
!                print*,'FTAU(J),-XLTAU,XLTAU,EXP(-XLTAU)',
!     &       FTAU(J),-XLTAU,XLTAU,EXP(-XLTAU)
!            endif    
          endif
        else
          FTAU(J)=0.0D0
        endif
      enddo
      if(U0.gt.0.D0) then
        ZFLUX = U0*FTAU(J1)*RFLECT/(1.d0+RFLECT)
      else
        ZFLUX = 0.d0
      endif
c
C------------------------------------------------------------------------
c  Take optical properties on CTM layers and convert to a photolysis
c  level grid corresponding to layer centres and boundaries. This is
c  required so that J-values can be calculated for the centre of CTM
c  layers; the index of these layers is kept in the jndlev array.
C------------------------------------------------------------------------
c
c  Set lower boundary and levels to calculate J-values at 
      J1=2*J1-1
      do j=1,lpar
        jndlev(j)=2*j
      enddo
c
c  Calculate column optical depths above each level, TTAU
      TTAU(NC+1)=0.0D0
      do J=NC,J1,-1
        I=(J+1)/2
        TTAU(J)=TTAU(J+1) + 0.5d0*DTAUX(I)
        jaddlv(j)=int(0.5d0*DTAUX(I)/dtaumax)
c  Subdivide cloud-top levels if required
! NOTE: Don't add more than DTAUSUB-1 (=9) sublevels (phs)
        if(jadsub(j).gt.0) then
          jadsub(j)=min(jaddlv(j)+1,int(dtausub))*(int(dsubdiv)-1) 
          jaddlv(j)=jaddlv(j)+jadsub(j)
        endif
      enddo
c
c  Calculate attenuated beam, FTAU, level boundaries then level centres
      FTAU(NC+1)=1.0D0
      do J=NC-1,J1,-2
        I=(J+1)/2
        FTAU(J)=FTAU(I)
      enddo
      do J=NC,J1,-2
        !if (abs(dimag(FTAU(J+1)))<1d-150) THEN
        !FTAU(J+1) = xplx(dble(FTAU(J+1)),0d0)
        !endif
        !if (abs(dimag(FTAU(J-1)))<1d-150) THEN
        !FTAU(J-1) = xplx(dble(FTAU(J-1)),0d0)
        !endif
        FTAU(J)=sqrt(FTAU(J+1)*FTAU(J-1))
!        if (isnan(FTAU(J))) then
!            print*,'FTAU(J) opmie 198',FTAU(J)
!        endif
      enddo
c
c  Calculate scattering properties, level centres then level boundaries
c  using an inverse interpolation to give correctly-weighted values
      do j=NC,J1,-2
        do i=1,MFIT
          pomegaj(i,j) = pomegaj(i,j/2)
        enddo
      enddo
      do j=J1+2,nc,2
        taudn = ttau(j-1)-ttau(j)
        tauup = ttau(j)-ttau(j+1)
        do i=1,MFIT
          pomegaj(i,j) = (pomegaj(i,j-1)*taudn + 
     $                    pomegaj(i,j+1)*tauup) / (taudn+tauup)
        enddo
      enddo
c  Define lower and upper boundaries
      do i=1,MFIT
        pomegaj(i,J1)   = pomegaj(i,J1+1)
        pomegaj(i,nc+1) = pomegaj(i,nc)
      enddo
c
C------------------------------------------------------------------------
c  Calculate cumulative total and define levels we want J-values at.
c  Sum upwards for levels, and then downwards for Mie code readjustments.
c
c     jaddlv(i)   Number of new levels to add between (i) and (i+1)
c     jaddto(i)   Total number of new levels to add to and above level (i)
c     jndlev(j)   Level needed for J-value for CTM layer (j)
c
C------------------------------------------------------------------------
c
c  Reinitialize level arrays
      do j=1,nc+1
        jaddto(j)=0
      enddo
c
      jaddto(J1)=jaddlv(J1)
      do j=J1+1,nc
        jaddto(j)=jaddto(j-1)+jaddlv(j)
      enddo

!==============================================================================
! KLUDGE TO LIMIT THE NUMBER OF ADDED LEVELS (phs, 7/1/08)
!
! PART 1: We need to replace the .gt. with .ge in this IF test
!
      if((jaddto(nc)+nc).GE.nl) then
         write(6,1500)  jaddto(nc)+nc, 'NL',NL
!
! PART 2: We just trim the largest JADDLV until the condition is satisfied 
!         instead of simply stopping.  Remove the STOP statement.
!
         !-------------------
         ! Prior to 7/1/08:
         !stop
         !-------------------

         ! trim
         do while( (SUM( jaddlv(J1:nc) ) + NC) >= NL )
            loc=maxloc(jaddlv)
            jaddlv(loc(1))=jaddlv(loc(1))-1
         enddo

         ! then refill JADDTO
         jaddto(J1)=jaddlv(J1)
         do j=J1+1,nc
            jaddto(j)=jaddto(j-1)+jaddlv(j)
         enddo
         
!        ! Debug: double check
!        write(6,*) jaddto(nc)+nc
!        if((jaddto(nc)+nc).gt.nl) 
!     &      write(6,*)'OPMIE kludge: trap not working'
!==============================================================================
      endif

c     write(6,1300) jndlev
c     write(6,1300) jaddto
      do i=1,lpar
         jndlev(i)=jndlev(i)+jaddto(jndlev(i)-1)
      enddo

      ! this is just a transposition of the jaddto vector (phs)
      jaddto(nc)=jaddlv(nc)
      do j=nc-1,J1,-1
         jaddto(j)=jaddto(j+1)+jaddlv(j)
      enddo
c     write(6,1300) jndlev
c     write(6,1300) jaddto
c
C---------------------SET UP FOR MIE CODE-------------------------------
c
c  Transpose the ascending TTAU grid to a descending ZTAU grid.
c  Double the resolution - TTAU points become the odd points on the
c  ZTAU grid, even points needed for asymm phase fn soln, contain 'h'.
c  Odd point added at top of grid for unattenuated beam   (Z='inf')
c  
c        Surface:   TTAU(1)   now use ZTAU(2*NC+1)
c        Top:       TTAU(NC)  now use ZTAU(3)
c        Infinity:            now use ZTAU(1)
c
c  Mie scattering code only used from surface to level NC
C------------------------------------------------------------------------
C
c  Initialise all Fast-J optical property arrays
      do k=1,N__
        do i=1,MFIT
          pomega(i,k) = 0.d0
        enddo
        ztau(k) = 0.d0
        fz(k)   = 0.d0
      enddo
c
c  Ascend through atmosphere transposing grid and adding extra points
      do j=J1,nc+1
        k = 2*(nc+1-j)+2*jaddto(j)+1
        ztau(k)= ttau(j)
        fz(k)  = ftau(j)
!        if (isnan(fz(k))) then
!           print*,'fz(k) opmie 309',fz(k),ftau(j)
!           CALL GEOS_CHEM_STOP
!        endif
        do i=1,MFIT
          pomega(i,k) = pomegaj(i,j)
        
        enddo
      enddo
c
c  Check profiles if desired
c      ND = 2*(NC+jaddto(J1)-J1)  + 3
c      if(kw.eq.1) call CH_PROF
c
C------------------------------------------------------------------------
c    Insert new levels, working downwards from the top of the atmosphere
c  to the surface (down in 'j', up in 'k'). This allows ztau and pomega
c  to be incremented linearly (in a +ve sense), and the flux fz to be
c  attenuated top-down (avoiding problems where lower level fluxes are
c  zero).
c
c    zk        fractional increment in level
c    dttau     change in ttau per increment    (linear, positive)
c    dpomega   change in pomega per increment  (linear)
c    ftaulog   change in ftau per increment    (exponential, normally < 1)
c
C------------------------------------------------------------------------
c
      do j=nc,J1,-1
          zk = 0.5d0/(1.d0+xplx(jaddlv(j)-jadsub(j)))
          dttau = (ttau(j)-ttau(j+1))*zk
          do i=1,MFIT
            dpomega(i) = (pomegaj(i,j)-pomegaj(i,j+1))*zk
          enddo
c  Filter attenuation factor - set minimum at 1.0d-05
          if(ftau(j+1).eq.0.d0) then
            ftaulog=0.d0
          else
            ftaulog = ftau(j)/ftau(j+1)
            if(ftaulog.lt.1.d-150) then
              ftaulog=1.0d-05
            else
              ftaulog=exp(log(ftaulog)*zk)
            endif
          endif
          k = 2*(nc-j+jaddto(j)-jaddlv(j))+1   !  k at level j+1
          l = 0
c  Additional subdivision of first level if required
          if(jadsub(j).ne.0) then
            l=jadsub(j)/nint(dsubdiv-1)
            zk2=1.d0/dsubdiv
            dttau2=dttau*zk2
            ftaulog2=ftaulog**zk2
            do i=1,MFIT
              dpomega2(i)=dpomega(i)*zk2
            enddo
            do ix=1,2*(jadsub(j)+l)
              ztau(k+1) = ztau(k) + dttau2 
        !if (abs(dimag(fz(k)))<1d-150) THEN
        !fz(k) = xplx(dble(fz(k)),0d0)
        !endif
        !if (abs(dimag(FTAU(J-1)))<1d-150) THEN
        !FTAU(J-1) = xplx(dble(FTAU(J-1)),0d0)
        !endif
              
              fz(k+1) = fz(k)*ftaulog2
              do i=1,MFIT
                pomega(i,k+1) = pomega(i,k) + dpomega2(i)
              enddo
              k = k+1
            enddo
          endif
          l = 2*(jaddlv(j)-jadsub(j)-l)+1
c
c  Add values at all intermediate levels
          do ix=1,l
            ztau(k+1) = ztau(k) + dttau
            fz(k+1) = fz(k)*ftaulog
            do i=1,MFIT
              pomega(i,k+1) = pomega(i,k) + dpomega(i)
            enddo
            k = k+1
          enddo
c
c  Alternate method to attenuate fluxes, fz, using 2nd-order finite
c  difference scheme - just need to comment in section below
c          ix = 2*(jaddlv(j)-jadsub(j))+1
c          if(l.le.0) then
c            l=k-ix-1
c          else
c            l=k-ix
c          endif
c          call efold(ftau(j+1),ftau(j),ix+1,fz(l))
c          if(jadsub(j).ne.0) then
c            k = 2*(nc-j+jaddto(j)-jaddlv(j))+1 !  k at level j+1
c            ix=2*(jadsub(j)+(jadsub(j)/nint(dsubdiv-1)))
c            call efold(ftau(j+1),fz(k+ix),ix,fz(k))
c          endif
c
      enddo
c
C---Update total number of levels and check doesn't exceed N__
      ND = 2*(NC+jaddto(J1)-J1)  + 3

!==============================================================================
! KLUDGE TO LIMIT THE NUMBER OF ADDED LEVELS (phs, 7/1/08)
!
! PART 3: Test to make sure that we haven't added more levels than the
!         dimension of the common block (i.e. ND <= N__).
!         
!         NOTE: this test should always be passed now that .ge. is 
!         used instead of .gt. in PART 1.
!
      if(nd.gt.N__) then
         write(6,1500) ND, 'N__',N__
         stop
      endif
!==============================================================================
c
C---Add boundary/ground layer to ensure no negative J's caused by
C---too large a TTAU-step in the 2nd-order lower b.c.
      ZTAU(ND+1) = ZTAU(ND)*1.000005d0
      ZTAU(ND+2) = ZTAU(ND)*1.000010d0
      zk=max(abs(U0),0.01d0)
      zk=exp(-ZTAU(ND)*5.d-6/zk)
      FZ(ND+1) = FZ(ND)*zk
      FZ(ND+2) = FZ(ND+1)*zk
!      if (isnan(FZ(ND))) then
!         print*,'FZ(ND) opmie 448',FZ(ND:ND+2),zk
!         CALL GEOS_CHEM_STOP
!      endif
      do I=1,MFIT
        POMEGA(I,ND+1)   = POMEGA(I,ND)
        POMEGA(I,ND+2)   = POMEGA(I,ND)
      enddo
      ND = ND+2
c
      ZU0 = U0
      ZREFL = RFLECT
c
C-----------------------------------------
      CALL MIESCT
C-----------------------------------------
c  Accumulate attenuation for selected levels
      l=2*(NC+jaddto(J1))+3
      do j=1,lpar
        k=l-(2*jndlev(j))
        if(k.gt.ND-2) then
          FMEAN(j) = 0.d0
        else
          FMEAN(j) = FJ(k)
        endif
      enddo
c
      return
 1000 format(1x,i3,3(2x,1pe10.4),1x,i3)
 1300 format(1x,50(i3))
 1500 format(' Too many levels in photolysis code: need ',i5,' but ',a,
     $       ' dimensioned as ',i5)
      END

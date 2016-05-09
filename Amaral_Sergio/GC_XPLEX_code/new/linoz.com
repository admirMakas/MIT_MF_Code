C $Id: linoz.com,v 1.1.1.1 2009/06/09 21:51:54 daven Exp $
C $Log: linoz.com,v $
C Revision 1.1.1.1  2009/06/09 21:51:54  daven
C Initial import
C
C Revision 2.23  2000/05/24 23:09:33  pjc
C Changed criteria for using Linoz: now must have [Ox]>150ppb AND Level>=9.
C
C Revision 2.10  2000/03/23 20:39:04  pjc
C Initial version created out of McLinden's original files.
C

C common block for linoz. Created by Philip Cameron-Smith, 00/1/14.

      INTEGER   nfields_linoz,nlevels_linoz,nlat_linoz,nmonths_linoz
      PARAMETER(nfields_linoz=7)  ! Number of linoz fields.
      PARAMETER(nlevels_linoz=25) ! Number of levels in linoz fields.
      PARAMETER(nlat_linoz=18)    ! Number of latitudes in linoz fields.
      PARAMETER(nmonths_linoz=12) !Number of months in linoz fields.

      TYPE (XPLEX) TPARM(nlevels_linoz,nlat_linoz,nmonths_linoz,nfields_linoz)
      TYPE (XPLEX) TLSTT(JJPAR,LLPAR,nfields_linoz)
      COMMON/linoz_fields/TPARM,TLSTT

      TYPE (XPLEX) linoz_min_alt      !Minimum altitude covered by linoz data.
      PARAMETER(linoz_min_alt=10) ! units=[km]
      INTEGER linoz_min_lev     ! Minimum GCM level linoz can cover.
      COMMON/linoz_levels/linoz_min_lev

C*PJC* Need to define the minimum Level at which Linoz can be used.
C      NB: Linoz data goes down to ~277mbar, so any part of a layer below 
C          this has effectively no Linoz chemistry.
      INTEGER Linoz_min_L
      PARAMETER(Linoz_min_L=9)
C*PJC* Define ozone tropopause, below which Linoz not used.
      TYPE (XPLEX)  Linoz_min_Ox
      PARAMETER(Linoz_min_Ox=150D-9) ! VMR, so 150D-9 = 150 ppb.


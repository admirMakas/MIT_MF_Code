! $Id: jv_index.f,v 1.1.1.1 2009/06/09 21:51:52 daven Exp $
      SUBROUTINE JV_INDEX
!
!******************************************************************************
!  Subroutine JV_INDEX computes the mapping between the CTM indices
!  (from "chem.dat") for J-values to the FAST-J indices (from "ratj.d")
!  for J-values.  (bmy, 10/5/98, 10/16/06)
!
!  NOTES:
!  (1 ) Assumes the ordering of a species with several branches in 
!        "ratj.d" is the same as in "chem.dat".
!  (2 ) Updated comments, cosmetic changes (bmy, 11/15/01)
!  (3 ) NAMESPEC is now NAMEGAS for SMVGEAR II.   We don't need to reference 
!        CMN anymore. Now loop from NCS = 1..NCSGAS (bdf, bmy, 4/8/03)
!  (4 ) Now reset NCS to NCSURBAN after loop (dbm, bmy, 10/16/06)
!******************************************************************************
!
      USE MYTYPE
      USE COMPLEXIFY
      IMPLICIT NONE

#     include "cmn_fj.h"  ! CMN_SIZE  
#     include "comode.h"  ! SMVGEAR II arrays

      ! Local variables
      INTEGER            :: I, IFNC, IBRCH, N, NK
      CHARACTER (LEN=4)  :: SPECNAME

      !=================================================================
      ! JV_INDEX begins here!
      !=================================================================

      ! Zero the RINDEX array
      RINDEX(:) = 0

      ! Loop over photolysis rxns (urban chemistry only)
      DO NCS = 1, NCSGAS
      DO I   = 1, NPHOT

         !==============================================================
         ! I        = Index of photo rxns    from "globchem.dat"
         ! NK       = Absolute rxn number (adds offset to I)
         ! SPECNAME = Name of species I,     from "globchem.dat"
         ! IBRCH    = Branch # of species I, from "globchem.dat"
         !==============================================================
         NK       = NRATES(NCS) + I
         SPECNAME = NAMEGAS(IRM(1,NK,NCS)) 
         IFNC     = DEFPRAT(NK,NCS) + 0.01d0
         IBRCH    = 10d0*( DEFPRAT(NK,NCS) - IFNC ) + 0.5d0

         !==============================================================
         ! N      = Index of photolysis reactions as listed in "ratj.d"
         ! RNAMES = Name of species N,            as listed in "ratj.d" 
         ! BRANCH = Branch number of species N,   as listed in "ratj.d" 
         !  
         ! If the species names and branch numbers from both "chem.dat" 
         ! and "ratj.d" match, then store N (the "ratj.d" index) in the 
         ! Ith element of RINDEX.
         !  
         ! Thus, when looping over I (the chem.dat" indices), as is 
         ! done in FJFUNC.F, RINDEX(I) will access the correct J-value 
         ! according to the ordering in "ratj.d".
         !==============================================================
         DO N = 1, JPPJ
            IF ( SPECNAME == RNAMES(N) .and. IBRCH == BRANCH(N) ) THEN
               RINDEX(I) = N

               WRITE ( 6, 100 ) I,         SPECNAME,  IBRCH, 
     &                          RINDEX(I), RNAMES(N), BRANCH(N)
 100           FORMAT('Harvard #: ', i3, 1x, a4, ' Branch: ', i2, 
     &                ' --->  Fast-J #: ', i3, 1x, a4, ' Branch: ',i2 )
               EXIT
            ENDIF
         ENDDO
      ENDDO  
      ENDDO

      ! Reset NCS to NCSURBAN for safety's sake (bmy, 10/16/06)
      NCS = NCSURBAN

      ! Return to calling program      
      END SUBROUTINE JV_INDEX

        !COMPILER-GENERATED INTERFACE MODULE: Mon Mar 26 14:19:48 2012
        MODULE FINDPAR_ERROR__genmod
          INTERFACE 
            SUBROUTINE FINDPAR_ERROR(ERROR,PAR_STR,STATUS,NM,MESSAGES,  &
     &ACTIONS)
              LOGICAL(KIND=4), INTENT(IN) :: ERROR
              CHARACTER(*), INTENT(IN) :: PAR_STR
              INTEGER(KIND=4), INTENT(INOUT) :: STATUS
              INTEGER(KIND=4), INTENT(INOUT) :: NM
              CHARACTER(*), INTENT(INOUT) :: MESSAGES(0:25)
              CHARACTER(*), INTENT(INOUT) :: ACTIONS(0:25)
            END SUBROUTINE FINDPAR_ERROR
          END INTERFACE 
        END MODULE FINDPAR_ERROR__genmod

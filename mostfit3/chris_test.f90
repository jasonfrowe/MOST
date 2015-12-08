PROGRAM RMORB
!ATTEMPTS TO PHASE AND CLIP THE ORBITAL CONTRIBUTION FROM MOST DATA.
!TAKES TIME COUNTS ERRORS AS INPUT.

USE CONSTANTS, ONLY: R8, I4, LGT
USE RINTER1
USE UTIL, ONLY: GETU, MEAN, SDEV,REALLOCAT
IMPLICIT NONE


REAL(R8), ALLOCATABLE, DIMENSION(:)  :: TIME,COUNTS,ERR,BX,BY,AVX,AVY,BPHASE
REAL(R8), POINTER, DIMENSION(:)            :: TT, TC, TBIN, CBIN
INTEGER(I4)                                                   :: IONML, IOOUT, NDAT,BI
REAL(R8)                :: SCUT, TM, TE, SUMM, MC, SC, ORBP,PBIN,PSHIFT,Y,DBIN,PHASE
INTEGER(I4)           :: HEADER, I, NFOLD,MXPHAS, J, K,NB,THRES,BINSZ
CHARACTER (500) :: FNAME, OUTNAME
CHARACTER (500) :: INFILE,OFILE
LOGICAL (LGT) :: GOOD

NAMELIST/ORBCUT/  HEADER,FNAME,OUTNAME,SCUT,ORBP,NFOLD,THRES,DBIN

!OPEN AND READ NAMELIST VARIABLES
IONML = 2004
OPEN  (IONML, FILE = 'rmorbs.nml', STATUS = 'OLD')
READ  (IONML, NML=ORBCUT)
CLOSE(IONML)

!GET THE DATA
INFILE = TRIM(FNAME)
CALL READDAT (INFILE, TIME, COUNTS, HEADER, ERR)

NDAT  = SIZE(TIME)
TM      = MINVAL(TIME)


!!****** BIN DATA HERE. DBIN IS THE BIN SIZE IN DAYS*********************
!BIN THE DATA FIRST IN DEFINED TIME BINS.
GOOD=.TRUE.
SUMM = 0
K = 1
ALLOCATE(TBIN(NDAT),CBIN(NDAT))
TBIN = -99.99D60
CBIN = -99.99D60
IF (DBIN > 0) THEN
DO WHILE(GOOD)
TE = TM + DBIN ! BIN WIDTH
!COUNT
BI = COUNT(((TM <= TIME).AND.(TIME < TE)))
IF ( (BI < THRES).AND.(TE < TIME(NDAT)) ) THEN
100  TE = TE + DBIN
   BI = COUNT(((TM <= TIME).AND.(TIME < TE)))
    IF ( ((BI < THRES).AND.(TE < TIME(NDAT)))  ) GOTO 100
ENDIF !PACK BINS -- brute force
ALLOCATE(TT(BI),TC(BI))
SUMM = 1
DO I = 1, NDAT
IF ((TM <= TIME(I)).AND.(TIME(I) < TE)) THEN
 TT(SUMM) = TIME(I)
  TC(SUMM) = COUNTS(I)
  SUMM = SUMM + 1
ENDIF
ENDDO

 TBIN(K) =   MINVAL(TT) + (MAXVAL(TT) - MINVAL(TT))/2.0 !TAKE CENTER TIME AS THE BIN TIME
 CBIN(K) =   MEAN(TC)
 K = K + 1
DEALLOCATE(TT,TC)
TM = TE
IF (TM > TIME(NDAT)) GOOD = .FALSE.
ENDDO

BINSZ = K - 1
K = 0
TBIN => REALLOCAT(TBIN,BINSZ)
CBIN => REALLOCAT(CBIN,BINSZ)
ELSE
 BINSZ = NDAT
 TBIN   = TIME
 CBIN   = COUNTS
ENDIF
!!******************************************************************************



!STEP 2 PHASE THE BINNED DATA TO THE ORBITAL PERIOD OF MOST = 0.0705044
!DAYS
TM      = MINVAL(TIME)
ALLOCATE(BX(BINSZ),BY(BINSZ),BPHASE(BINSZ))
BPHASE =(TBIN - MINVAL(TBIN))/ORBP

MXPHAS = MAXVAL(INT(BPHASE))
WRITE(*,'(1X,A9,1X,I6,1X,A26)') "THERE ARE", MXPHAS, "ORBITAL PHASES IN THE DATA"

IOOUT = GETU()
OPEN (IOOUT, FILE = TRIM(OUTNAME),STATUS="UNKNOWN")
!CHECK THAT THERE AREN'T MORE PHASE FOLDS THAN THERE ARE PHASES
IF (NFOLD > MXPHAS) GOTO 110
PBIN = 0.005 ! KEEP THIS FIXED -- SCAN THROUGH PHASES AT THIS STEP

!ADAPTED FROM JASON ROWE
DO I =1, NDAT
!     pivot point is currect phase
!     we want all data ± nfold/2
       PHASE = (TIME(I) - TM)/ORBP
        NB=0
        DO J=1,BINSZ
           IF((BPHASE(J) > PHASE-REAL(NFOLD)/2.).AND.&
               (BPHASE(J) < PHASE+REAL(NFOLD)/2.)) THEN
              NB=NB+1
              BX(NB)=BPHASE(J)-INT(BPHASE(J))
              BY(NB)=CBIN(J)
           ENDIF
        ENDDO
      !       now shift data so important phase is in the middle
         PSHIFT=PHASE - INT(PHASE) - 0.5
         DO J=1,NB
           BX(J)=BX(J)-PSHIFT
           IF(BX(J) < 0.0) BX(J)=BX(J)+1.0
           IF(BX(J)  >1.0) BX(J)=BX(J)-1.0
         ENDDO

         ALLOCATE(AVX(NB), AVY(NB))
         AVX = 0.0
         AVY = 0.0

!         GATHER DATA INTO PBINS

               K = 0
          DO J =1,NB
           IF(ABS(0.5-BX(J)) < PBIN) THEN
               K = K +1
               AVX(K)=BX(J)
               AVY(K)=BY(J)
            ENDIF
          ENDDO
! K IS THE CURRENT SIZE OF THE PBIN BINS

!!!!!!!!!!!! WHEN BINNING SOMETIMES K = 0
  ! SIGMA CLIP THE DATA.
!sigma clip from ROWE.
        call sigclip(K,avx,avy,SCUT)
                Y=0
        DO J=1,K
           Y=Y+AVY(J)
         ENDDO
!USE THIS AVERAGE CORRECTION TO "SMOOTH" THE ORBITAL/NEAR ORBITAL
!COMPONENTS            Y=Y/REAL(K)
      WRITE(IOOUT,*) TIME(I), COUNTS (I) - Y, ERR(I)
     ! write(*,*) TIME(I), COUNTS(I), Y
     !  pause
      DEALLOCATE(AVX,AVY)
ENDDO

CLOSE(IOOUT)
DEALLOCATE(TIME, COUNTS, ERR,BX,BY,TBIN,CBIN,BPHASE)
STOP
110 WRITE(*,*) "THERE ARE MORE PHASE FOLDS THAN THERE ARE PHASES"
      WRITE(*,*) "LOWER NFOLD"
DEALLOCATE(TIME, COUNTS, ERR,BX,BY,TBIN,CBIN,BPHASE)
STOP
END PROGRAM RMORB 
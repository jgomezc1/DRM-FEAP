CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   ----------------P R O G R A M   P I R O--------------------------  C
C                    Finite Element Method                             C
C                                                                      C
C  GENERATES EFFECTIVE LOADS TO BE APPLIED IN THIRD PARTY SOFTWARE     C
C  BASED UPON A DOMAIN REDUCTION METHOD APPROACH                       C
C                                                                      C
C  DYNAMIC MEMORY ALLOCATION                                           C
C  EAFIT UNIVERSITY                                                    C
C  APPLIED MECHANICS LAB                                               C
C  VIVIANA DIAZ-VELEZ; JUAN GOMEZ                                      C
C  MARCH 02/2016                                                       C
C                                                                      C
C  UNIX VERSION                                                        C
C                                                                      C
C            P R O B L E M   P A R A M E T E R S                       C
C                                                                      C
C  NUMNP     :NUMBER OF NODAL POINTS                                   C
C  NUMEL     :NUMBER OF ELEMENTS                                       C
C  NUMAT     :NUMBER OF MATERIAL PROFILES                              C
C  TMAX      :SIZE OF THE TIME WINDOW.                                 C
C  NINCR     :NUMBER OF INCREMENTS                                     C
C  NDOFDIM   :PROBLEM DIMENSIONALITY                                   C
C  NMNE      :MAXIMUM NUMBER OF NODES PER ELEMENT                      C
C  NMDOFEL   :MAXIMUM NUMBER OF DOF PER ELEMENT                        C
C  NMDOFNP   :MAXIMUM NUMBER OF DOF PER NODE                           C
C  NMPR      :MAXIMUM NUMBER OF MATERIAL PROPERTIES IN A PROFILE       C
C  NIPR      :MAXIMUM NUMBER OF INTEGER MAT PROPERTIES IN A PROFILE    C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C                                                                      C
C                      MAIN PROGRAM STARTS                             C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      PROGRAM PIRO
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION AWAVE(5)
C
      ALLOCATABLE ID(:,:),NDOFN(:),COORD(:,:),MATP(:),NMATP(:),NIMTP(:),
     1            AMATE(:,:),IMTEI(:,:),NNE(:),NDOFEL(:),
     2            INDE(:),IELCON(:,:),LM(:,:),RHSG(:,:)
C
      CHARACTER*80 HED
      CHARACTER*10 FILENAME
C
C     *****************************************************************C
C     ***          P R O B L E M  F I L E S                          **C
C     *****************************************************************C
C
      IIN=5
      IOUT=4
      IDAT=3
      WRITE(*,*) 'INPUT THE JOB NAME(max 10 characters):'
      READ(*,*) FILENAME
c      FILENAME='IPIRO'
      LST=LEN_TRIM(FILENAME)
      OPEN(UNIT=IIN,FILE =FILENAME(1:LST)//".inp",FORM='FORMATTED')
      OPEN(UNIT=IOUT,FILE=FILENAME(1:LST)//".dat",FORM='FORMATTED')
C
C     *****************************************************************C
C     ***                 I N P U T   P H A S E                      **C
C     *****************************************************************C
C
C     READS PROBLEM DEFINITION PARAMETERS.
C
      READ(IIN,*,IOSTAT=INOUTSTATUS) HED
      IF(INOUTSTATUS.LT.0) STOP "***COULD NOT OPEN FILE"
C
      READ(IIN,     *) NUMNP,NUMEL,NUMAT,TMAX,DT,NINCR,NDOFDIM,NMNE,
     1                 NMDOFEL,NMDOFNP,NMPR,NIPR,ISOW
C
      WRITE(IOUT,1900) HED
      WRITE(IOUT,2000) NUMNP,NUMEL,NUMAT,TMAX,DT,NINCR,NDOFDIM,NMNE,
     1                 NMDOFEL,NMDOFNP,NMPR,NIPR
C
      ALLOCATE(ID(NMDOFNP,NUMNP),NDOFN(NUMNP),COORD(NDOFDIM,NUMNP),
     1         MATP(NUMEL),NMATP(NUMAT),NIMTP(NUMAT),AMATE(NMPR,NUMAT),
     2         IMTEI(NIPR,NUMAT),NNE(NUMEL),NDOFEL(NUMEL),
     3         INDE(NUMNP),IELCON(NMNE,NUMEL),LM(NMDOFEL,NUMEL))
C
      CALL CLEARIM(ID,NMDOFNP,NUMNP)
      CALL CLEARIV(NDOFN,NUMNP)
      CALL CLEAR(COORD,NDOFDIM,NUMNP)
C
      CALL NODINP(NUMNP,ID,NMDOFNP,NDOFDIM,COORD,NDOFN,NEQ,INDE,IIN,
     &            IOUT)
      ALLOCATE(RHSG(NEQ,NINCR))
      CALL CLEAR(RHSG,NEQ,NINCR)
C
C     CLEARS STORAGE.
C
      CALL CLEARIV(MATP,NUMEL)
      CALL CLEARIV(NMATP,NUMAT)
      CALL CLEARIV(NIMTP,NUMAT)
      CALL CLEAR(AMATE,NMPR,NUMAT)
      CALL CLEARIM(IMTEI,NIPR,NUMAT)
      CALL CLEARIV(NNE,NUMEL)
      CALL CLEARIV(NDOFEL,NUMEL)
      CALL CLEARIM(IELCON,NMNE,NUMEL)
      CALL CLEARIM(LM,NMDOFEL,NUMEL)
C
C     READS MODEL.
C
      CALL MATINP(NUMAT,NMPR,NMATP,AMATE,NIPR,NIMTP,IMTEI,AWAVE,IWAVE,
     1            IIN,IOUT)
      CALL ELEINP(NUMEL,NNE,NDOFEL,NMNE,MATP,IELCON,IIN,IOUT)
C
C     CREATES ASSEMBLY LIST (DME OPERATOR).
C
      CALL ASSEMLIS(NUMNP,NUMEL,NMNE,NMDOFEL,NMDOFNP,NNE,NDOFN,
     1              IELCON,LM,ID,INDE,IOUT)
C
      CALL GSTFASSEM(NUMNP,NUMEL,NUMAT,NNE,NMNE,NMDOFEL,NDOFDIM,
     1               IELCON,NDOFN,NDOFEL,MATP,NMATP,NMPR,NIMTP,
     2               NIPR,INDE,AMATE,AWAVE,IWAVE,ISOW,IMTEI,COORD,
     3               LM,NEQ,NINCR,RHSG,DT,IOUT)
C
C     PRINT THE OUTPUT FILE
C
      SELECT CASE (IWAVE)
        CASE(1:2)
            CALL GENLOADPS(NUMNP,NINCR,NMDOFNP,INDE,RHSG)
        CASE(3)
            CALL GENLOADSH(NUMNP,NINCR,NMDOFNP,INDE,RHSG)
      END SELECT
!
      CALL cpu_time(T1)
!
      WRITE(*,*) 'TIME RUN PIRO',T1
C
      WRITE(*,3050), FILENAME
C
      STOP
C
C     *****************************************************************C
C     ***         S L N  O U T P U T  P H A S E                      **C
C     *****************************************************************C
C
 1900 FORMAT(//,10X,'P R O B L E M   N A M E',10X,A80,//)
 2000 FORMAT (///,
     1    ' C O N T R O L  I N F O R M A T I O N',//,
     2    '   NUMBER OF NODAL POINTS',10(' .'),' (NUMNP)=',I5,//,
     3    '   NUMBER OF ELEMENTS',12(' .'),'(NUMEL)=',I5,//,
     4    '   NUMBER OF MATERIALS       ',9(' .'),'(NUMAT)=',I5,//,
     5    '   SIZE OF THE TIME WINDOW',9(' .'),'(TMAX)=',F10.5,//,
     6    '   TIME STEP',9(' .'),'(DT)=',F10.5,//,
     7    '   NUMBER OF LOAD INCREMENTS',9(' .'),'(NINCR)=',I5,//,
     8    '   DEGREE OF FREEDOM DIMENSION',6(' .'),'(NDOFDIM)=',I1,//,
     9    '   MAX.NUMBER OF NODES IN AN ELE.',6(' .'),'(NMNE)=',I2,//,
     1    '   MAX.NUMBER OF DOF IN AN ELE.',6(' .'),'(NMDOFEL)=',I2,//,
     2    '   MAX.NUMBER OF DOF IN AN NODE',6(' .'),'(NMDOFNP)=',I2,//,
     3    '   MAX.NUMBER OF MAT PROPERTIES.',6(' .'),'(NMPR)=',I2,//,
     4    '   MAX.NUMBER OF INT MATP ROPERTIES.',6(' .'),'(NIPR)=',I2)

C
 3050 FORMAT(/,'JOB=',A8,1X,'COMPLETED')
C
      END

C      MODULE PREPROC
C      CONTAINS
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C SUBROUTINE NODINP                                                    C
C                                                                      C
C   NUMNP      :Numberof nodal points                                  C 
C   N          :Nodal identifier                                       C
C   NDOFN(N)   :Number of degrees of freedom at the current node       C
C   ID(II,I)   :Boundary condition identifiers at the current node     C
C   COORD(JJ,I):Coordinates of the current node                        C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE NODINP(NUMNP,ID,MXDOFNP,MXDOFDIM,COORD,NDOFN,NEQ,INDE,
     &                  IIN,IOUT)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION ID(MXDOFNP,NUMNP),COORD(MXDOFDIM,NUMNP),NDOFN(NUMNP),
     1          INDE(NUMNP)
C
      CALL CLEAR(COORD,MXDOFDIM,NUMNP)
      CALL CLEARIM(ID,MXDOFNP,NUMNP)
C
C     READ AND GENERATE NODAL POINT DATA
C
      WRITE(IOUT,1000)
      WRITE(IOUT,1005)
C
      DO I=1,NUMNP
        READ(IIN,     *) INDE(I),NDOFN(I),(COORD(JJ,I),JJ=1,MXDOFDIM)

        WRITE(IOUT,   *) INDE(I),NDOFN(I),(COORD(JJ,I),JJ=1,MXDOFDIM)

      END DO
C
C     ASSIGN EQUATION NUMBERS TO ACTIVE DOF
C
      ICOUNT=1
      DO K1=1,NUMNP
        DO K2=1,NDOFN(K1)
          IF(ID(K2,K1).EQ.0) THEN
            ID(K2,K1)=ICOUNT
            ICOUNT=ICOUNT+1
          ELSE
            ID(K2,K1)=0
          END IF
        END DO
      END DO
      NEQ=ICOUNT-1
C
C     WRITE EQUATION NUMBERS
C
      WRITE(IOUT,1020)
      WRITE(IOUT,1025)
      WRITE(IOUT,1030) (INDE(N),(ID(I,N),I=1,MXDOFNP),N=1,NUMNP)

C
 1000 FORMAT(//,6X,' N O D A L  D A T A',//)
 1005 FORMAT('ID',2X,'NDOF',2X,'BC-X',2X,'BC-Y',2X,6X,'COORD-X',6X,
     1       'COORD-Y',/)
 1010 FORMAT(I5,4X,I2,4X,I2,4X,I2,6X,F5.2,6X,F5.2)
 1020 FORMAT(/,'E Q U A T I O N  N U M B E R S',/)
 1025 FORMAT(/,'NODE',6X,'DEGREES OF FREEDOM',/,
     1      'NUMBER',8X,'U',9X,'V'/)
 1030 FORMAT(I5,5X,I5,5X,I5)
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C SUBROUTINE MATINP                                                    C
C                                                                      C
C   N       :ID of material profile                                    C
C   NMATP() :Number of real material parameters at the current profile C
C   NIMTP() :Number of integer material parmt´s at the current profile C
C   AMATE(,):Array with real material parameters for each profile      C
C   IMTEI(,):Array with integer material parameters for each profile   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE MATINP(NUMAT,NMPR,NMATP,AMATE,NIPR,NIMTP,IMTEI,AWAVE,
     1                  IWAVE,IIN,IOUT)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION AMATE(NMPR,NUMAT),NMATP(NUMAT),IMTEI(NIPR,NUMAT),
     1          NIMTP(NUMAT),AWAVE(5)
C
      CALL CLEAR( AMATE,NMPR,NUMAT)
      CALL CLEARIM(IMTEI,NIPR,NUMAT)
C
C     READS AND GENERATES MATERIAL DATA
C
      WRITE (IOUT,1000)
      WRITE (IOUT,3000)
      DO I=1,NUMAT
        READ  (IIN, *) N,NMATP(N),NIMTP(N),(AMATE(JJ,N),JJ=1,NMATP(N)),
     1                                     (IMTEI(JJ,N),JJ=1,NIMTP(N))
        WRITE (IOUT,3020) N,NMATP(N),NIMTP(N),(AMATE(JJ,N),JJ=1,
     1                    NMATP(N))
      END DO
C
      WRITE(IOUT,3025)
      WRITE(IOUT,3030)
      READ(IIN,     *) IWAVE, AWAVE(1),AWAVE(2),AWAVE(3),AWAVE(4),
     1                 AWAVE(5)
      WRITE (IOUT,3035) IWAVE,AWAVE(1),AWAVE(2),AWAVE(3),AWAVE(4),
     1                  AWAVE(5)
C
 1000 FORMAT(//,4X,'M A T E R I A L  D A T A',/)
 3000 FORMAT('Mat-Id',2X,'R-Prop',5X,'I-Prop',20X,
     1'Vp',20X,'Vs',12X,'Density',12X,'Alfa',16X,'Beta'//)
 3020 FORMAT(I5,3X,I5,3X,I5,12X,5(5X,F13.3))
 3025 FORMAT(//,4X,'W A V E  D A T A',/)
 3030 FORMAT('Wav-Id',2X,'W-P1',5X,'W-P2',5X,'W-P3',5X,'W-P4',5X,
     1       'W-P5'//)
 3035 FORMAT(I5,3X,5(5X,F10.5))
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE ELEINP                                                C
C                                                                      C
C     NUMNP     :Number of nodal points                                C
C     NUMEL     :Number of elements                                    C
C     NNE()     :Number of nodes at the current element                C
C     IELT()    :Array of element type identifiers                     C
C     NDOFEL()  :Number of degrees of freedom at the current element   C
C     NMNE      :Maximum number of nodes per element                   C
C     MATP()    :Material pofile at the current element                C
C     IELCON(,) :Element connectivity array                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE ELEINP(NUMEL,NNE,NDOFEL,NMNE,MATP,IELCON,
     1                  IIN,IOUT)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION NNE(NUMEL),MATP(NUMEL),IELCON(NMNE,NUMEL),
     1          NDOFEL(NUMEL)
C
C     Creates IELCON()
C
      WRITE(IOUT,1000)
      WRITE(IOUT,1001)
      WRITE(IOUT,1002)
      DO I=1,NUMEL
        READ(IIN,*) M,NDOFEL(I),MATP(I),NNE(I),
     1  (IELCON(J,I),J=1,NNE(I))
        WRITE(IOUT,1010) M,NDOFEL(I),MATP(I),NNE(I),
     1  (IELCON(J,I),J=1,NNE(I))
      END DO
C
 1000 FORMAT(///,8X,'E L E M E N T  I N F O R M A T I O N',//)
 1001 FORMAT('ID  TYPE NDOF MAT NNEL NODE NODE NODE NODE NODE NODE
     1 NODE NODE')
 1002 FORMAT(25X'1    2    3    4    5    6   7    8',/)
 1010 FORMAT(5(2X,I2),X,8(2X,I3))
C
      RETURN
C
      END
C

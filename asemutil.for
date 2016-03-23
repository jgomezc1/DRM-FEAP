cc      MODULE ASEMUTIL
cc      CONTAINS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C         --G E N E R A L  A S S E M B L Y  S U B R O U T I N E S--    C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE ASSEMLIS                                              C
C                                                                      C
C     Creates the DME() operator LM()                                  C
C                                                                      C
C     NUMNP    :Number of nodal points                                 C
C     NUMEL    :Number of elements                                     C
C     MXNE     :Maximum number of nodes per element                    C
C     MXDOFDIM :Maximum degree of freedom dimension                    C
C     NNE()    :Number of nodes at the current element                 C
C     NDOFN()  :Number of degrees of freedom at the current node       C
C     IELCON(,):Nodal connectivity array                               C
C     LM(,)    :DME opeartor                                           C               
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE ASSEMLIS(NUMNP,NUMEL,MXNE,MXDOFEL,MXDOFNP,NNE,
     1                    NDOFN,IELCON,LM,ID,INDE,IOUT)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION NNE(NUMEL),IELCON(MXNE,NUMEL),LM(MXDOFEL,NUMEL),
     1          ID(MXDOFNP,NUMNP),NDOFN(NUMNP),INDE(NUMNP)
C
      IFLAG=0
      DO I=1,NUMEL
        K3=0
        DO J=1,NNE(I)
          IVAL=IELCON(J,I)
          CALL SEARCHPOS(INDE,NUMNP,IVAL,IFLAG,IPOS)
          K2I=NDOFN(IPOS)
          DO K2=1,K2I
            LM(K3+K2,I)=ID(K2,IPOS)
          END DO
          K3=K3+K2I
        END DO
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE GSTFASEM                                              C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE GSTFASSEM(NUMNP,NUMEL,NUMAT,NNE,NMNE,MXDOFEL,MXDOFDIM,
     1                     IELCON,NDOFN,NDOFEL,MATP,NMATP,NMPR,NIMTP,
     2                     NIPR,INDE,AMATE,AWAVE,IWAVE,ISOW,IMTEI,
     3                     COORD,LM,NEQ,NINCR,RHSG,DT,IOUT)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ZERO=0.D0)
C
C     INTEGER ARRAYS
C
      DIMENSION NNE(NUMEL),IELCON(NMNE,NUMEL),MATP(NUMEL),
     1          LM(MXDOFEL,NUMEL),LML(MXDOFEL),NDOFN(NUMNP),
     2          NDOFEL(NUMEL), NMATP(NUMAT), INDE(NUMNP),NIMTP(NUMAT),
     3          IMTEI(NIPR,NUMAT)
C
C     REAL ARRAYS
C
      DIMENSION AMATE(NMPR,NUMAT),COORD(MXDOFDIM,NUMNP),
     1          ELCOOR(MXDOFDIM,NMNE),RHSG(NEQ,NINCR),PROPS(NMPR),
     2          IPROPS(NIPR),RHS(MXDOFEL,NINCR),AWAVE(5)
C
C     LOCAL ARRAYS
C
      CALL CLEAR(RHSG,NEQ,NINCR)
      CALL CLEAR(RHS,MXDOFEL,NINCR)
      CALL CLEARIV(LML,MXDOFEL)
C
C     RETRIEVES ELEMENT MATERIAL PROPERTIES
C
      DO I=1,NUMEL
C
        NN=NNE(I)
        NPROPS=NMATP(MATP(I))
        NIPROPS=NIMTP(MATP(I))
        DO K1=1,NPROPS
          PROPS(K1)=AMATE(K1,MATP(I))
        END DO
C
        DO K1=1,NIPROPS
          IPROPS(K1)=IMTEI(K1,MATP(I))
        END DO
C
        ND=NDOFEL(I)
C
C       RETRIEVES ELEMENT NODAL COORDINATES
C       AND LOCAL ASSEMBLY LIST
C
        K3=0
        DO J=1,NNE(I)
          IVAL=IELCON(J,I)
          IFLAG=0
          CALL SEARCHPOS(INDE,NUMNP,IVAL,IFLAG,IPOS)
          ELCOOR(1,J)=COORD(1,IPOS)
          ELCOOR(2,J)=COORD(2,IPOS)
C
          K1=NDOFN(IPOS)
          DO K2=1,K1
C
            LML(K3+K2)=LM(K3+K2,I)
          END DO
          K3=K3+K1
        END DO
C
C       CALLS EMBEDED USER ELEMENT SUBROUTINE UEL.
C
        SELECT CASE (IWAVE)
            CASE(1:2)
                CALL UELINCOTPS(RHS,ND,PROPS,NPROPS,IPROPS,NIPROPS,
     1                          AWAVE,IWAVE,ISOW,ELCOOR,MXDOFDIM,NN,
     2                          NINCR,DT,IOUT)
            CASE (3)
                CALL UELINCOTSH(RHS,ND,PROPS,NPROPS,IPROPS,NIPROPS,
     1                          AWAVE,IWAVE,ISOW,ELCOOR,MXDOFDIM,NN,
     2                          NINCR,DT,IOUT)
        END SELECT
C
C       ASSEMBLE GLOBAL FORCE VECTOR
C
        DO II=1,ND
          KK=LML(II)
          DO JJ=1,NINCR
            RHSG(KK,JJ)=RHSG(KK,JJ)+RHS(II,JJ)
          END DO
        END DO
C
        IF(I.EQ.NUMEL) WRITE(*,*)'END ELEMENT',I
C
      END DO
C
      RETURN
C
      END
C
cc      END MODULE ASEMUTIL

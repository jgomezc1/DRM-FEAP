cc      MODULE UELUTIL
cc      CONTAINS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C       --U S E R   E L E M E N T    S U B R O U T I N E S---          C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE UELINCOTPS(RHS,NDOFEL,PROPS,NPROPS,IPROP,IPROPS,AWAVE,C
C    1                      IWAVE,COORDS,MCRD,NNODE,NINCR,DT,IOUT)     C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE UELINCOTPS(RHS,NDOFEL,PROPS,NPROPS,IPROP,IPROPS,AWAVE,
     1                      IWAVE,ISOW,COORDS,MCRD,NNODE,NINCR,DT,IOUT)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,HALF=0.5D0,ONE=1.0D0,TWO=2.D0,THREE=3.D0,
     1           NGPTS=9,NDI=3,NTENS=4)
C
      DIMENSION RHS(NDOFEL,NINCR),PROPS(NPROPS),AWAVE(5),
     1          COORDS(MCRD,NNODE),DDSDDE(NTENS,NTENS),B(NTENS,NDOFEL),
     2          BT(NDOFEL,NTENS),FRST1(NDOFEL,NDOFEL),XP(2,NGPTS),
     3          XW(NGPTS),AUX1(NTENS,NDOFEL),R00(NDOFEL,NINCR),
     4          AMATRX(NDOFEL,NDOFEL)
C
      DIMENSION IPROP(IPROPS)
C
      IDF=IPROP(1)
C
C     Clears arrays.
C
      CALL CLEAR(RHS,NDOFEL,NINCR)
      CALL CLEAR(DDSDDE,NTENS,NTENS)
      CALL CLEAR(B,NTENS,NDOFEL)
      CALL CLEAR(BT,NDOFEL,NTENS)
      CALL CLEAR(FRST1,NDOFEL,NDOFEL)
      CALL CLEAR(XP,2,NGPTS)
      CALL CLEARV(XW,NGPTS)
      CALL CLEAR(AUX1,NTENS,NDOFEL)
      CALL CLEAR(R00,NDOFEL,NINCR)
      CALL CLEAR(AMATRX,NDOFEL,NDOFEL)
C
C     Generates Gauss points and weights
C     and loops around all Gauss points.
C
      CALL GPOINTS3X3(XP,XW)
      NGPT=9
C
      DO NN=1,NGPT
C
        RII=XP(1,NN)
        SII=XP(2,NN)
        ALF=XW(NN  )
C
C       Computes B matrix.
C
        SELECT CASE (NNODE)
            CASE(4)
                CALL STDM4(NNODE,NDOFEL,NTENS,COORDS,MCRD,B,DDET,RII,
     1                     SII,XBAR)
            CASE(8)
                CALL STDM8(NNODE,NDOFEL,NTENS,COORDS,MCRD,B,DDET,RII,
     1                     SII,XBAR)
            CASE(9)
                CALL STDM9(NNODE,NDOFEL,NTENS,COORDS,MCRD,B,DDET,RII,
     1                     SII,XBAR)
        END SELECT
C
C       Assembles stiffness matrix.
C
        CALL UMAT(DDSDDE,NDI,NTENS,PROPS,NPROPS)
        CALL MMULT(DDSDDE,NTENS,NTENS,B,NTENS,NDOFEL,AUX1)
        CALL MTRAN(B,NTENS,NDOFEL,BT)
        CALL MMULT(BT,NDOFEL,NTENS,AUX1,NTENS,NDOFEL,FRST1)
C
C       Considers Gauss weight and Jacobian determinant
C       into partial stiffness matrix.
C
        CALL SMULT(FRST1,NDOFEL,NDOFEL,ALF*DDET*XBAR)
C
C       Updates stiffness matrix.
C
        CALL UPDMAT(AMATRX,NDOFEL,NDOFEL,FRST1)
C
C       Clears arrays for new Gauss point
C
        CALL CLEAR(DDSDDE,NTENS,NTENS)
        CALL CLEAR(FRST1,NDOFEL,NDOFEL)
C
      END DO
C
C     Computes the incoming wave field and makes it into an effective
C     forces vector R00(n).
C
      CALL INCOMINGDISPS(IDF,AMATRX,COORDS,PROPS,NPROPS,MCRD,AWAVE,
     1                   IWAVE,ISOW,NNODE,NDOFEL,R00,NINCR,DT,IOUT)
C
C     Updates the effective loads compatible with
C     the incoming wave field.
C
      DO I=1,NDOFEL
        DO J=1,NINCR
          RHS(I,J)=R00(J,I)
        END DO
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE UELINCOTSH(RHS,NDOFEL,PROPS,NPROPS,IPROP,IPROPS,AWAVE,  C
C    1                      IWAVE,COORDS,MCRD,NNODE,NINCR,DT,IOUT)     C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE UELINCOTSH(RHS,NDOFEL,PROPS,NPROPS,IPROP,IPROPS,AWAVE,
     1                      IWAVE,ISOW,COORDS,MCRD,NNODE,NINCR,DT,IOUT)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ZERO=0.D0,HALF=0.5D0,ONE=1.0D0,TWO=2.D0,THREE=3.D0,
     1           NGPTS=9,NDI=2,NTENS=2)
C
      DIMENSION RHS(NDOFEL,NINCR),PROPS(NPROPS),AWAVE(5),
     1          COORDS(MCRD,NNODE),DDSDDE(NTENS,NTENS),B(NTENS,NDOFEL),
     2          BT(NDOFEL,NTENS),FRST1(NDOFEL,NDOFEL),XP(2,NGPTS),
     3          XW(NGPTS),AUX1(NTENS,NDOFEL),R00(NDOFEL,NINCR),
     4          AMATRX(NDOFEL,NDOFEL)
C
      DIMENSION IPROP(IPROPS)
C
      IDF=IPROP(1)
C
C     Clears arrays.
C
      CALL CLEAR(RHS,NDOFEL,NINCR)
      CALL CLEAR(DDSDDE,NTENS,NTENS)
      CALL CLEAR(B,NTENS,NDOFEL)
      CALL CLEAR(BT,NDOFEL,NTENS)
      CALL CLEAR(FRST1,NDOFEL,NDOFEL)
      CALL CLEAR(XP,2,NGPTS)
      CALL CLEARV(XW,NGPTS)
      CALL CLEAR(AUX1,NTENS,NDOFEL)
      CALL CLEAR(R00,NDOFEL,NINCR)
      CALL CLEAR(AMATRX,NDOFEL,NDOFEL)
C
C     Generates Gauss points and weights
C     and loops around all Gauss points.
C
      CALL GPOINTS3X3(XP,XW)
      NGPT=9
C
      DO NN=1,NGPT
C
        RII=XP(1,NN)
        SII=XP(2,NN)
        ALF=XW(NN  )
C
C       Computes B matrix.
        SELECT CASE (NNODE)
            CASE(4)
                CALL STDM4SH(NNODE,NDOFEL,NTENS,COORDS,MCRD,B,DDET,RII,
     1                       SII,XBAR)
            CASE(9)
                CALL STDM9SH(NNODE,NDOFEL,NTENS,COORDS,MCRD,B,DDET,RII,
     1                       SII,XBAR)
        END SELECT
C
C       Assembles stiffness matrix.
C
        CALL UMATSH(DDSDDE,NDI,NTENS,PROPS,NPROPS)
        CALL MMULT(DDSDDE,NTENS,NTENS,B,NTENS,NDOFEL,AUX1)
        CALL MTRAN(B,NTENS,NDOFEL,BT)
        CALL MMULT(BT,NDOFEL,NTENS,AUX1,NTENS,NDOFEL,FRST1)
C
C       Considers Gauss weight and Jacobian determinant
C       into partial stiffness matrix.
C
        CALL SMULT(FRST1,NDOFEL,NDOFEL,ALF*DDET*XBAR)
C
C       Updates stiffness matrix.
C
        CALL UPDMAT(AMATRX,NDOFEL,NDOFEL,FRST1)
C
C       Clears arrays for new Gauss point
C
        CALL CLEAR(DDSDDE,NTENS,NTENS)
        CALL CLEAR(FRST1,NDOFEL,NDOFEL)
C
      END DO
C
C     Computes the incoming wave field and makes it into an effective
C     forces vector R00(n).
C
      CALL INCOMINGDISSH(IDF,AMATRX,COORDS,PROPS,NPROPS,MCRD,AWAVE,
     1                   IWAVE,ISOW,NNODE,NDOFEL,R00,NINCR,DT,IOUT)
C
C     Updates the effective loads compatible with
C     the incoming wave field.
C
      DO I=1,NDOFEL
        DO J=1,NINCR
          RHS(I,J)=R00(J,I)
        END DO
      END DO
C
      RETURN
C
      END
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE INCOMINGDISPS(IDF,AMATRX,COORDS,PROPS,NPROPS,MCRD,    C
C    1                         AWAVE,IWAVE,NNODE,NDOFEL,R00,NINCR,     C
C    2                         DT,IOUT)                                C
C                                                                      C
C     Computes effective forces consistent with the incoming wave      C
C     field at time t.                                                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE INCOMINGDISPS(IDF,AMATRX,COORDS,PROPS,NPROPS,MCRD,
     1                         AWAVE,IWAVE,ISOW,NNODE,NDOFEL,R00,NINCR,
     2                         DT,IOUT)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(NGPTS=6,ONE=1.0D0,ONENEG=-1.0D0)
C
      DIMENSION AMATRX(NDOFEL,NDOFEL),COORDS(MCRD,NNODE),PROPS(NPROPS),
     1          R00(NDOFEL,NINCR),U0(NDOFEL,NINCR),
     2          AWAVE(5),IELCON(4,3),LEE(6),
     3          U00(NDOFEL,NINCR),RP(NDOFEL,NINCR),UT(NDOFEL),
     4          RT(NDOFEL)
C
      ALLOCATABLE LB(:),LE(:)
C
      CALL CLEAR(U0,NDOFEL,NINCR)
C
C     Evaluate incoming displacements at time t.
C
      SELECT CASE (ISOW)
        CASE(0)
             CALL DISPINCTW(U0,COORDS,NNODE,PROPS,NPROPS,MCRD,NDOFEL,
     1                      AWAVE,IWAVE,DT,NINCR,IOUT)
        CASE(1)
             CALL DISPINCTS(U0,COORDS,NNODE,PROPS,NPROPS,MCRD,NDOFEL,
     1                      AWAVE,IWAVE,DT,NINCR,IOUT)
      END SELECT
C
C     Convert displacements and velocities into
C     effective forces at time t.
C
      DO I=1,NDOFEL
        AMATRX(I,I)=0.0D0
        DO JJ=1,NINCR
          U00(I,JJ)=U0(I,JJ)
        END DO
      END DO
C
      CALL LOCIEL(IELCON)
      CALL SELLEE(LEE,IDF)
      CALL SELNNF(NNODE,NNB,NNE)
      ALLOCATE(LB(2*NNB),LE(2*NNE))
C
C     Creates list of exterior DOFs
C
      DO I=1,NNE
        J=LEE(I)
        LE(2*I-1)=2*J-1
        LE(2*I)=2*J
      END DO
C
C     Creates list of boundary DOFs
C
      DO I=1,NNB
        J=IELCON(IDF,I)
        LB(2*I-1)=2*J-1
        LB(2*I)=2*J
      END DO
C
C     Identifies DOF in LB.
C
      DO J=1,NDOFEL
        IFLAG=0
        CALL BELONGS(LB,2*NNB,J,IFLAG)
        IF (IFLAG.EQ.1) THEN
           DO I=1,NINCR
             U00(J,I)=0.D0
           END DO
        END IF
      END DO
C
      CALL CLEAR(RP,NDOFEL,NINCR)
C
      DO II=1,NINCR
        DO JJ=1,NDOFEL
          UT(JJ)=U00(JJ,II)
        END DO
        CALL MAVEC(AMATRX,NDOFEL,NDOFEL,UT,RT)
        DO JJ=1,NDOFEL
          RP(JJ,II)=RT(JJ)
          RT(JJ)=0.D0
        END DO
      END DO
C
C     Assembles effective forces for the boundary nodes.
C
      DO J=1,NDOFEL
        IFLAG=0
        CALL BELONGS(LB,2*NNB,J,IFLAG)
        IF (IFLAG.EQ.1) THEN
           DO II=1,NINCR
             R00(J,II)=-RP(J,II)
           END DO
        END IF
      END DO

      CALL CLEAR(U00,NDOFEL,NINCR)
C
      DO I=1,NDOFEL
        DO JJ=1,NINCR
          U00(I,JJ)=U0(I,JJ)
        END DO
      END DO
C
C     Identifies DOF in LE.
C
      DO J=1,NDOFEL
        IFLAG=0
        CALL BELONGS(LE,2*NNE,J,IFLAG)
        IF (IFLAG.EQ.1) THEN
           DO I=1,NINCR
             U00(J,I)=0.D0
           END DO
        END IF
      END DO
C
      CALL CLEAR(RP,NDOFEL,NINCR)
C
      DO II=1,NINCR
        DO JJ=1,NDOFEL
          UT(JJ)=U00(JJ,II)
        END DO
        CALL MAVEC(AMATRX,NDOFEL,NDOFEL,UT,RT)
        DO JJ=1,NDOFEL
          RP(JJ,II)=RT(JJ)
          RT(JJ)=0.D0
        END DO
      END DO
C
C     Assembles effective forces for the exterior nodes.
C
      DO J=1,NDOFEL
        IFLAG=0
        CALL BELONGS(LE,2*NNE,J,IFLAG)
        IF (IFLAG.EQ.1) THEN
           DO II=1,NINCR
             R00(J,II)=RP(II,J)
           END DO
        END IF
      END DO
C
      RETURN
C
      END
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C      SUBROUTINE INCOMINGDISSH(IDF,AMATRX,COORDS,PROPS,NPROPS,MCRD,   C
C     1                         AWAVE,IWAVE,NNODE,NDOFEL,R00,NINCR,DT, C
C     2                         IOUT)                                  C
C                                                                      C
C   Computes effective forces consistent with the incoming wave field  C
C   at time t.                                                         C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE INCOMINGDISSH(IDF,AMATRX,COORDS,PROPS,NPROPS,MCRD,
     1                         AWAVE,IWAVE,ISOW,NNODE,NDOFEL,R00,NINCR,
     2                         DT,IOUT)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(NGPTS=6,ONE=1.0D0,ONENEG=-1.0D0)
C
      DIMENSION AMATRX(NDOFEL,NDOFEL),COORDS(MCRD,NNODE),PROPS(NPROPS),
     1          R00(NDOFEL,NINCR),U0(NDOFEL,NINCR),
     2          AWAVE(5),IELCON(4,3),LEE(6),
     3          U00(NDOFEL,NINCR),RP(NDOFEL,NINCR),UT(NDOFEL),
     4          RT(NDOFEL)
C
      ALLOCATABLE LB(:),LE(:)
C
      CALL CLEAR(U0,NDOFEL,NINCR)
C
C     Evaluate incoming displacements at time t.
      SELECT CASE (ISOW)
        CASE(0)
             CALL DISPINCTW(U0,COORDS,NNODE,PROPS,NPROPS,MCRD,NDOFEL,
     1                      AWAVE,IWAVE,DT,NINCR,IOUT)
        CASE(1)
             CALL DISPINCTS(U0,COORDS,NNODE,PROPS,NPROPS,MCRD,NDOFEL,
     1                      AWAVE,IWAVE,DT,NINCR,IOUT)
      END SELECT
C
C     Convert displacements and velocities into
C     effective forces at time t.
C
      DO I=1,NDOFEL
        AMATRX(I,I)=0.0D0
        DO JJ=1,NINCR
          U00(I,JJ)=U0(I,JJ)
        END DO
      END DO
C
      CALL LOCIEL(IELCON)
      CALL SELLEE(LEE,IDF)
      CALL SELNNF(NNODE,NNB,NNE)
      ALLOCATE(LB(NNB),LE(NNE))
C
C     Creates list of exterior DOFs
C
      DO I=1,NNE
        J=LEE(I)
        LE(I)=J
      END DO
C
C     Creates list of boundary DOFs
C
      DO I=1,NNB
        J=IELCON(IDF,I)
        LB(I)=J
      END DO
C
C     Identifies DOF in LB.
C
      DO J=1,NDOFEL
        IFLAG=0
        CALL BELONGS(LB,NNB,J,IFLAG)
        IF (IFLAG.EQ.1) THEN
           DO I=1,NINCR
             U00(J,I)=0.D0
           END DO
        END IF
      END DO
C
      CALL CLEAR(RP,NDOFEL,NINCR)
C
      DO II=1,NINCR
        DO JJ=1,NDOFEL
          UT(JJ)=U00(JJ,II)
        END DO
        CALL MAVEC(AMATRX,NDOFEL,NDOFEL,UT,RT)
        DO JJ=1,NDOFEL
          RP(JJ,II)=RT(JJ)
          RT(JJ)=0.D0
        END DO
      END DO
C
C     Assembles effective forces for the boundary nodes.
C
      DO J=1,NDOFEL
        IFLAG=0
        CALL BELONGS(LB,NNB,J,IFLAG)
        IF (IFLAG.EQ.1) THEN
           DO II=1,NINCR
             R00(J,II)=-RP(J,II)
           END DO
        END IF
      END DO

      CALL CLEAR(U00,NDOFEL,NINCR)
C
      DO I=1,NDOFEL
        DO JJ=1,NINCR
          U00(I,JJ)=U0(I,JJ)
        END DO
      END DO
C
C     Identifies DOF in LE.
C
      DO J=1,NDOFEL
        IFLAG=0
        CALL BELONGS(LE,NNE,J,IFLAG)
        IF (IFLAG.EQ.1) THEN
           DO I=1,NINCR
             U00(J,I)=0.D0
           END DO
        END IF
      END DO
C
      CALL CLEAR(RP,NDOFEL,NINCR)
C
      DO II=1,NINCR
        DO JJ=1,NDOFEL
          UT(JJ)=U00(JJ,II)
        END DO
        CALL MAVEC(AMATRX,NDOFEL,NDOFEL,UT,RT)
        DO JJ=1,NDOFEL
          RP(JJ,II)=RT(JJ)
          RT(JJ)=0.D0
        END DO
      END DO
C
C     Assembles effective forces for the exterior nodes.
C
      DO J=1,NDOFEL
        IFLAG=0
        CALL BELONGS(LE,NNE,J,IFLAG)
        IF (IFLAG.EQ.1) THEN
           DO II=1,NINCR
             R00(J,II)=RP(J,II)
           END DO
        END IF
      END DO
C
      RETURN
C
      END
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE DISPINCTW(U0,COORDS,NNODE,PROPS,NPROPS,MCRD,NDOFEL,   C
c                          AWAVE,IWAVE,DT,NINCR,IOUT)                  C
C                                                                      C
C     Evaluates and assembles the displacement field.                  C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE DISPINCTW(U0,COORDS,NNODE,PROPS,NPROPS,MCRD,NDOFEL,
     1                     AWAVE,IWAVE,DT,NINCR,IOUT)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION U0(NDOFEL,NINCR),US1(NINCR),US2(NINCR),AWAVE(5),
     1          COORDS(MCRD,NNODE),PROPS(NPROPS)
C
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0,TWO=2.0D0,FOUR=4.0D0)
C
      CALL CLEARV(US1,NINCR)
      CALL CLEARV(US2,NINCR)
C
C     Incident wave properties.
C
      FC=AWAVE(1)
      TMAX=AWAVE(2)
      TINI=AWAVE(3)
      AMPR=AWAVE(4)
      PHI=AWAVE(5)
C
C     Elastic properties
C
      AALFA=PROPS(1)
      ABETA=PROPS(2)
      RO=PROPS(3)
C
C     Compute incoming dispalcements for all the nodes in
C     an element.
C
      DO I=1,NNODE
        X1=COORDS(1,I)
        X3=COORDS(2,I)
        II=I
C
        SELECT CASE (IWAVE)
          CASE (1)
            CALL PWDISPT(FC,TINI,AMPR,PHI,AALFA,ABETA,X1,X3,US1,US2,
     1                   NINCR,DT)
C
              DO JJ=1,NINCR
                U0(2*II-1,JJ)=US1(JJ)
                U0(2*II  ,JJ)=US2(JJ)
              END DO
C
          CASE (2)
            CALL SVWDISPT(FC,TINI,AMPR,PHI,AALFA,ABETA,X1,X3,US1,US2,
     1                    NINCR,DT)
C
              DO JJ=1,NINCR
                U0(2*II-1,JJ)=US1(JJ)
                U0(2*II  ,JJ)=US2(JJ)
              END DO
C
          CASE (3)
            CALL SHWDISPT(FC,TINI,AMPR,PHI,ABETA,X1,X3,US1,NINCR,DT)
C
              DO JJ=1,NINCR
                U0(II,JJ)=US1(JJ)
              END DO
        END SELECT
C
      END DO
C
      RETURN
C
      END
C
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE DISPINCTS(U0,COORDS,NNODE,PROPS,NPROPS,MCRD,            C
C                       AWAVE,IWAVE,DT,NINCR,IOUT)                     C
C                                                                      C
C   Evaluates and assembles the displacement field.                    C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE DISPINCTS(U0,COORDS,NNODE,PROPS,NPROPS,MCRD,NDOFEL,
     1                     AWAVE,IWAVE,DT,NINCR,IOUT)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION U0(NDOFEL,NINCR),US1(NINCR),US2(NINCR),AWAVE(5),
     1          COORDS(MCRD,NNODE),PROPS(NPROPS)
C
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0,TWO=2.0D0,FOUR=4.0D0)
C
      CALL CLEARV(US1,NINCR)
      CALL CLEARV(US2,NINCR)
C
C     Compute incoming dispalcements for all the nodes in
C     an element.
C
      DO I=1,NNODE
        X1=COORDS(1,I)
        X3=COORDS(2,I)
        CALL SEARCHBEM(X1,X3,NINCR,DT,US1,US2)
C
        SELECT CASE (IWAVE)
            CASE(1:2)
              DO JJ=1,NINCR
                U0(2*I-1,JJ)=US1(JJ)
                U0(2*I  ,JJ)=US2(JJ)
              END DO
            CASE (3)
              DO JJ=1,NINCR
                U0(I,JJ)=US1(JJ)
              END DO
        END SELECT
C
      END DO

      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE SEARCHBEM(X1,X3,NINCR,DT,US1,US2)                     C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SEARCHBEM(X1,X3,NTFEM,DTFEM,US1,US2)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION US1(NTFEM),US2(NTFEM)
C
      ALLOCATABLE UBEM1(:),UBEM2(:)
C
      CHARACTER*10 FILENAME
C
      IIN=1
      FILENAME='DATABEM'
      LST=LEN_TRIM(FILENAME)
      OPEN(UNIT=IIN,FILE =FILENAME(1:LST)//".inp",FORM='FORMATTED')
C
      READ(IIN,     *) NTBEM,DTBEM,NNOBS,TE
C
      ALLOCATE(UBEM1(NTBEM),UBEM2(NTBEM))
      
      CALL READBEM(X1,X3,TE,NTBEM,UBEM1,UBEM2,IIN)
      CALL INTERPO(NTBEM,NTFEM,DTFEM,DTBEM,UBEM1,UBEM2,US1,US2)
C
      CLOSE(IIN)
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE READBEM(X1,X3,TE,NNOBS,NTBEM,UBEM1,UBEM2,IIN)         C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE READBEM(X1,X3,TE,NTBEM,UBEM1,UBEM2,IIN)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION UBEM1(NTBEM),UBEM2(NTBEM)
      CHARACTER*2 TITLE

      TOL=TE/10.0D0
C
      C=0
      DO WHILE (C .EQ. 0)
        READ(IIN,*) TITLE
        IF (TITLE .EQ. 'PU') THEN
            READ(IIN,*)X,Y
            IF (X1.LE.X+TOL.AND.X1.GE.X-TOL) THEN
                IF (X3.LE.Y+TOL.AND.X3.GE.Y-TOL) THEN
                    C=1
                END IF
            END IF
        END IF
      END DO

      DO II=1,NTBEM
        READ(IIN,*)T,UX,UY
        UBEM1(II)=UX
        UBEM2(II)=UY
      END DO

C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE INTERPO(NTBEM,NTFEM,DTFEM,DTBEM,UBEM1,UBEM2,US1,US2)  C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE INTERPO(NTBEM,NTFEM,DTFEM,DTBEM,UBEM1,UBEM2,US1,US2)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION US1(NTFEM),US2(NTFEM),UBEM1(NTBEM),UBEM2(NTBEM)
C

      TOL=DTFEM/10.0D0 ! ESTA TOLERANCIA DEPENDE DEL DT DE BEM

      DO I=1,NTFEM
        T=DTFEM*(I-1)
        DO II=1,(NTBEM-1)
            T1=DTBEM*(II-1)
            T2=DTBEM*(II)
            A11=UBEM1(II)
            A12=UBEM1(II+1)
            A21=UBEM2(II)
            A22=UBEM2(II+1)
C
            IF (T.LE.T2.AND.T.GE.T1) THEN
                P1=(A12-A11)/(T2-T1)
                B1=A12-P1*T2
                A1=P1*T+B1
                US1(I)=A1
C
                P2=(A22-A21)/(T2-T1)
                B2=A22-P2*T2
                A2=P2*T+B2
                US2(I)=A2
                EXIT
            END IF
        END DO
C
        DO II=1,NTBEM
            T1=DTBEM*(II-1)
            T2=DTBEM*(II)
            IF (T.LE.T1+TOL.AND.T.GE.T1-TOL) THEN
                US1(I)=UBEM1(II)
                US2(I)=UBEM2(II)
                EXIT
            END IF
        END DO
      END DO

      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE PWDISPT(FC,TMAX,TINI,AMP,PHI,AALFA,ABETA,RO,X1,X3,US1,C
C    1                   US2,VS1,VS2,AS1,AS2,NF)                       C
C                                                                      C
C     Computes the displacements vector for a P-Wave incident at a     C
C     point X1,X3.                                                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE PWDISPT(FC,TINI,AMP,PHI,AALFA,ABETA,X1,X3,US1,US2,
     1                   NINCR,dtao)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION UPPI(NINCR),WPPI(NINCR),UPPR(NINCR),WPPR(NINCR),
     1          UPSR(NINCR),WPSR(NINCR),US1(NINCR),US2(NINCR)
C
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0,TWO=2.0D0,FOUR=4.0D0)
C
      PI=4.0D0*DATAN(ONE)
C
      PHIR=PHI*PI/180.D0
      PL=DSIN(PHIR)/AALFA
      ZHIR=DASIN(ABETA*PL)
C
      SNX=DSIN(PHIR)/AALFA*X1
      CNZ=DCOS(PHIR)/AALFA*X3
      ENX=DSIN(ZHIR)/ABETA*X1
      ENZ=DCOS(ZHIR)/ABETA*X3
C
      F1=ONE/ABETA/ABETA-TWO*PL*PL
      F1S=F1**2
      F2=FOUR*PL*PL*(DCOS(PHIR)/AALFA)*(DCOS(ZHIR)/ABETA)
C
C     Compute reflection coefficients (PP,PS) and related constants.
C
      PP=(F2-F1S)/(F2+F1S)
      PS=(FOUR*AALFA/ABETA*PL*(DCOS(PHIR)/AALFA)*F1)/(F1S+F2)
C
C     Incident P-wave field.
C
      do ii=1,NINCR
        tao=(SNX-CNZ)-((dtao*(ii-1)-TINI))
        upft=(2.*(PI*tao*fc)**2-1.)*DEXP(-(PI*tao*fc)**2)
        uppi(ii)= AMP*DSIN(PHIR)*upft
        wppi(ii)=-AMP*DCOS(PHIR)*upft
      end do
C
C     Reflected P-wave field.
C
      do ii=1,NINCR
        tao=(SNX+CNZ)-((dtao*(ii-1)-TINI))
        upft=(2.*(PI*tao*fc)**2-1.)*DEXP(-(PI*tao*fc)**2)
        uppr(ii)=AMP*DSIN(PHIR)*PP*upft
        wppr(ii)=AMP*DCOS(PHIR)*PP*upft
      end do
C
C     Reflected S-wave field.
C
      do ii=1,NINCR
        tao=(ENX+ENZ)-((dtao*(ii-1)-TINI))
        upft=(2.*(PI*tao*fc)**2-1.)*DEXP(-(PI*tao*fc)**2)
        upsr(ii)= AMP*DCOS(ZHIR)*PS*upft
        wpsr(ii)=-AMP*DSIN(ZHIR)*PS*upft
      end do
C
C     Total wave field.
C
      do ii=1,NINCR
        US1(ii)=UPPI(ii)+UPPR(ii)+UPSR(ii)
        US2(ii)=WPPI(ii)+WPPR(ii)+WPSR(ii)
      end do
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C    SUBROUTINE SVWDISPT(FC,TMAX,TINI,AMP,ZHI,AALFA,ABETA,RO,X1,X3,US1,C
C    1                   US2,VS1,VS2,AS1,AS2,NF)                       C
C                                                                      C
C     Computes the displacements vector for a P-Wave incident at a     C
C     point X1,X3.                                                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SVWDISPT(FC,TINI,AMP,ZHI,AALFA,ABETA,X1,X3,US1,US2,
     1                    NINCR,dtao)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION USVI(NINCR),WSVI(NINCR),USPR(NINCR),WSPR(NINCR),
     1          USVR(NINCR),WSVR(NINCR),US1(NINCR),US2(NINCR)
C
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0,TWO=2.0D0,FOUR=4.0D0)
C
      PI=4.0D0*DATAN(ONE)
C
      ZHIR=ZHI*PI/180.D0
      PL=DSIN(ZHIR)/ABETA
      PHIR=DASIN(AALFA*PL)
C
      SNX=DSIN(PHIR)/AALFA*X1
      CNZ=DCOS(PHIR)/AALFA*X3
      ENX=DSIN(ZHIR)/ABETA*X1
      ENZ=DCOS(ZHIR)/ABETA*X3
C
      F1=ONE/ABETA/ABETA-TWO*PL*PL
      F1S=F1**2
      F2=FOUR*PL*PL*(DCOS(PHIR)/AALFA)*(DCOS(ZHIR)/ABETA)
C
C     Compute reflection coefficients PP,PS and related constants.
C
      SS=(F1S-F2)/(F2+F1S)
      SP=(FOUR*ABETA/AALFA*PL*(DCOS(ZHIR)/ABETA)*F1)/(F1S+F2)
C
C     Incident SV-wave field.
C
      do ii=1,NINCR
        tao=(ENX-ENZ)-((dtao*(ii-1)-TINI))
        usft=(2.*(PI*tao*fc)**2-1.)*DEXP(-(PI*tao*fc)**2)
        usvi(ii)= AMP*DCOS(ZHIR)*usft
        wsvi(ii)= AMP*DSIN(ZHIR)*usft
      end do
C
C     Reflected P-wave field.
C
      do ii=1,NINCR
        tao=(SNX+CNZ)-((dtao*(ii-1)-TINI))
        usft=(2.*(PI*tao*fc)**2-1.)*DEXP(-(PI*tao*fc)**2)
        uspr(ii)= AMP*DSIN(PHIR)*SP*usft
        wspr(ii)= AMP*DCOS(PHIR)*SP*usft
      end do
C
C     Reflected SV-wave field.
C
      do ii=1,NINCR
        tao=(ENX+ENZ)-((dtao*(ii-1)-TINI))
        usft=(2.*(PI*tao*fc)**2-1.)*DEXP(-(PI*tao*fc)**2)
        usvr(ii)= AMP*DCOS(ZHIR)*SS*usft
        wsvr(ii)=-AMP*DSIN(ZHIR)*SS*usft
      end do
C
C     Total wave field.
C
      do ii=1,NINCR
        US1(ii)=usvi(ii)+uspr(ii)+usvr(ii)
        US2(ii)=wsvi(ii)+wspr(ii)+wsvr(ii)
      end do
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE SHWDISPT(FC,TINI,AMP,ZHI,AALFA,ABETA,X1,X3,US1,       C
C    1                    NINCR,dtao)                                  C
C                                                                      C
C     Computes the displacements vector for a P-Wave incident at a     C
C     point X1,X3.                                                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SHWDISPT(FC,TINI,AMP,ZHI,ABETA,X1,X3,US1,NINCR,dtao)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION USHI(NINCR),USHR(NINCR),US1(NINCR)
C
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0,TWO=2.0D0,FOUR=4.0D0)
C
      PI=4.0D0*DATAN(ONE)
C
      ZHIR=ZHI*PI/180.D0
C
      ENX=DSIN(ZHIR)/ABETA*X1
      ENZ=DCOS(ZHIR)/ABETA*X3

C     Incident SH-wave field.
C
      do ii=1,NINCR
        tao=(ENX-ENZ)-((dtao*(ii-1)-TINI))
        usft=(2.*(PI*tao*fc)**2-1.)*DEXP(-(PI*tao*fc)**2)
        ushi(ii)= AMP*usft
      end do
C
C     Reflected SH-wave field.
C
      do ii=1,NINCR
        tao=(ENX+ENZ)-((dtao*(ii-1)-TINI))
        usft=(2.*(PI*tao*fc)**2-1.)*DEXP(-(PI*tao*fc)**2)
        ushr(ii)= AMP*usft
      end do
C
C     Total wave field.
C
      do ii=1,NINCR
        US1(ii)=ushi(ii)+ushr(ii)
      end do
C
      RETURN
C
      END
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE SELNNF(NNODE,NNI,NNE)                                 C
C                                                                      C
C     Creates the number of nodes at internal and external face        C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SELNNF(NNODE,NNI,NNE)
C
      IMPLICIT REAL*8(A-H,O-Z)
C

      SELECT CASE (NNODE)
C
        CASE(4)
            NNI=2
            NNE=2
C
        CASE(8)
            NNI=3
            NNE=5
C
        CASE(9)
            NNI=3
            NNE=6
C
      END SELECT
C
      RETURN
C
      END
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE SELLEE(LEE,IDF)                                       C
C                                                                      C
C     Creates the list of external nodes                               C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SELLEE(LEE,IDF)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION LEE(6)

      SELECT CASE (IDF)
C
        CASE (1)
         LEE(1)=3
         LEE(2)=4
         LEE(3)=6
         LEE(4)=7
         LEE(5)=8
         LEE(6)=9
C
        CASE (2)
         LEE(1)=1
         LEE(2)=4
         LEE(3)=5
         LEE(4)=7
         LEE(5)=8
         LEE(6)=9
C
        CASE (3)
         LEE(1)=1
         LEE(2)=2
         LEE(3)=5
         LEE(4)=6
         LEE(5)=8
         LEE(6)=9
C
        CASE (4)
         LEE(1)=2
         LEE(2)=3
         LEE(3)=5
         LEE(4)=6
         LEE(5)=7
         LEE(6)=9
C
      END SELECT
C
      RETURN
C
      END
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE LOCIEL(IELCON)                                        C
C                                                                      C
C     Creates the connectivity array for an 8-noded 2D element which   C
C     is considered like an assembly of 4 3-noded 1D elements          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE LOCIEL(IELCON)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION IELCON(4,3)
C
      IELCON(1,1)=1
      IELCON(1,2)=2
      IELCON(1,3)=5
C
      IELCON(2,1)=2
      IELCON(2,2)=3
      IELCON(2,3)=6
C
      IELCON(3,1)=3
      IELCON(3,2)=4
      IELCON(3,3)=7
C
      IELCON(4,1)=4
      IELCON(4,2)=1
      IELCON(4,3)=8
C
      RETURN
C
      END
C
cc      END MODULE UELUTIL

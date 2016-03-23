cc      MODULE FEMUTIL
cc      CONTAINS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C11111111122222222223333333333444444444455555555556666666666777777777777
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C         --G E N E R A L  F E M  S U B R O U T I N E S--              C
C                                                                      C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C STRAIN DISPLACEMENT INTERPOLATORS                                    C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE STDM4:GENERATES THE STRAIN-DISPLACEMENT MATRIX B      C
C     AND JACOBIAN DETERMINANT DDET AT THE POINT r ,s                  C
C                                                 i  j                 C
C     FOR AN 4-NODED 2D ELEMENT-PLANE STRAIN                           C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE STDM4(NNE,NDOFEL,NTENS,XX,MCRD,B,DDET,R,S,XBAR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=0.5D0)
C
      DIMENSION XX(MCRD,NNE),B(NTENS,NDOFEL),P(2,NNE),AUX1(2,NNE),
     1          XJ(2,2),XJI(2,2)
C
C     Initialize arrays
C
      XBAR=ONE
      CALL CLEAR(B,NTENS,NDOFEL)
      CALL CLEAR(XJ,2,2)
      CALL CLEAR(XJI,2,2)
      CALL CLEAR(P,2,NNE)
C
C     Shape functions derivatives w.r.t natural coordinates
C
      CALL SFDER4(NNE,R,S,P)
C
C     Computes the Jacobian operator and its determinant at point (r,s)
C
      CALL JACOPER(NNE,XJ,XX,P)
      DDET=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1)
C
C     Computes the inverse of the Jacobiam operator
C
      CALL JACINVE(XJ,DDET,XJI)
C
C     Jacobian Inverse times Shape Functions derivatives w.r.t natural coordinates
C     to produce shape function derivatives w.r.t x,y coordinates.
C
      CALL MMULT(XJI,2,2,P,2,NNE,AUX1)
C
C     Assembles B matrix for a
C     Cosserat element.
C
      DO I=1,NNE
        II=2*(I-1)+1
        B(1,II)=AUX1(1,I)
        B(2,II+1)=AUX1(2,I)
        B(4,II)=AUX1(2,I)
        B(4,II+1)=AUX1(1,I)
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE STDM4SH:GENERATES THE STRAIN-DISPLACEMENT MATRIX B    C
C     AND JACOBIAN DETERMINANT DDET AT THE POINT r ,s                  C
C                                                 i  j                 C
C     FOR AN 4-NODED 2D ELEMENT- SH WAVE                               C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE STDM4SH(NNE,NDOFEL,NTENS,XX,MCRD,B,DDET,R,S,XBAR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=0.5D0)
C
      DIMENSION XX(MCRD,NNE),B(NTENS,NDOFEL),P(2,NNE),AUX1(2,NNE),
     1          XJ(2,2),XJI(2,2)
C
C     Initialize arrays
C
      XBAR=ONE
      CALL CLEAR(B,NTENS,NDOFEL)
      CALL CLEAR(XJ,2,2)
      CALL CLEAR(XJI,2,2)
      CALL CLEAR(P,2,NNE)
C
C     Shape functions derivatives w.r.t natural coordinates
C
      CALL SFDER4(NNE,R,S,P)
C
C     Computes the Jacobian operator and its determinant at point (r,s)
C
      CALL JACOPER(NNE,XJ,XX,P)
      DDET=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1)
C
C     Computes the inverse of the Jacobiam operator
C
      CALL JACINVE(XJ,DDET,XJI)
C
C     Jacobian Inverse times Shape Functions derivatives w.r.t natural coordinates
C     to produce shape function derivatives w.r.t x,y coordinates.
C
      CALL MMULT(XJI,2,2,P,2,NNE,AUX1)
C
C     Assembles B matrix for a
C     Cosserat element.
C
      DO I=1,NNE
        B(1,I)=AUX1(1,I)
        B(2,I)=AUX1(2,I)
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE STDM6:GENERATES THE STRAIN-DISPLACEMENT MATRIX B      C
C     AND JACOBIAN DETERMINANT DDET AT THE POINT r ,s                  C
C                                                 i  j                 C
C     FOR AN 6-NODED 2D ELEMENT-PLANE STRAIN                           C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE STDM6(NNE,NDOFEL,NTENS,XX,MCRD,B,DDET,R,S,XBAR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=0.5D0)
C
      DIMENSION XX(MCRD,NNE),B(NTENS,NDOFEL),P(2,NNE),AUX1(2,NNE),
     1          XJ(2,2),XJI(2,2)
C
C     Initialize arrays
C
      XBAR=ONE
      CALL CLEAR(B,NTENS,NDOFEL)
      CALL CLEAR(XJ,2,2)
      CALL CLEAR(XJI,2,2)
      CALL CLEAR(P,2,NNE)
C
C     Shape functions derivatives w.r.t natural coordinates
C
      CALL SFDER6(NNE,R,S,P)
C
C     Computes the Jacobian operator and its determinant at point (r,s)
C
      CALL JACOPER(NNE,XJ,XX,P)
      DDET=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1)
C
C     Computes the inverse of the Jacobiam operator
C
      CALL JACINVE(XJ,DDET,XJI)
C
C     Jacobian Inverse times Shape Functions derivatives w.r.t natural coordinates
C     to produce shape function derivatives w.r.t x,y coordinates.
C          
      CALL MMULT(XJI,2,2,P,2,NNE,AUX1)
C
C     Assembles B matrix for a
C     Cosserat element.
C
      DO I=1,NNE
        II=2*(I-1)+1
        B(1,II)=AUX1(1,I)
        B(2,II+1)=AUX1(2,I)
        B(4,II)=AUX1(2,I)
        B(4,II+1)=AUX1(1,I)
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE STDM8:GENERATES THE STRAIN-DISPLACEMENT MATRIX B      C
C     AND JACOBIAN DETERMINANT DDET AT THE POINT r ,s                  C
C                                                 i  j                 C
C     FOR AN 8-NODED 2D ELEMENT-PLANE STRAIN                           C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE STDM8(NNE,NDOFEL,NTENS,XX,MCRD,B,DDET,R,S,XBAR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=0.5D0)
C
      DIMENSION XX(MCRD,NNE),B(NTENS,NDOFEL),P(2,NNE),AUX1(2,NNE),
     1          XJ(2,2),XJI(2,2)
C
C     Initialize arrays
C
      XBAR=ONE
      CALL CLEAR(B,NTENS,NDOFEL)
      CALL CLEAR(XJ,2,2)
      CALL CLEAR(XJI,2,2)
      CALL CLEAR(P,2,NNE)
C
C     Shape functions derivatives w.r.t natural coordinates
C
      CALL SFDER8(NNE,R,S,P)
C
C     Computes the Jacobian operator and its determinant at point (r,s)
C
      CALL JACOPER(NNE,XJ,XX,P)
      DDET=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1)
C
C     Computes the inverse of the Jacobiam operator
C
      CALL JACINVE(XJ,DDET,XJI)
C
C     Jacobian Inverse times Shape Functions derivatives w.r.t natural coordinates
C     to produce shape function derivatives w.r.t x,y coordinates.
C          
      CALL MMULT(XJI,2,2,P,2,NNE,AUX1)
C
C     Assembles B matrix for a
C     Cosserat element.
C
      DO I=1,NNE
        II=2*(I-1)+1
        B(1,II)=AUX1(1,I)
        B(2,II+1)=AUX1(2,I)
        B(4,II)=AUX1(2,I)
        B(4,II+1)=AUX1(1,I)
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE STDM9:GENERATES THE STRAIN-DISPLACEMENT MATRIX B      C
C     AND JACOBIAN DETERMINANT DDET AT THE POINT r ,s                  C
C                                                 i  j                 C
C     FOR AN 9-NODED 2D ELEMENT-PLANE STRAIN                           C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE STDM9(NNE,NDOFEL,NTENS,XX,MCRD,B,DDET,R,S,XBAR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=0.5D0)
C
      DIMENSION XX(MCRD,NNE),B(NTENS,NDOFEL),P(2,NNE),AUX1(2,NNE),
     1          XJ(2,2),XJI(2,2)
C
C     Initialize arrays
C
      XBAR=ONE
      CALL CLEAR(B,NTENS,NDOFEL)
      CALL CLEAR(XJ,2,2)
      CALL CLEAR(XJI,2,2)
      CALL CLEAR(P,2,NNE)
C
C     Shape functions derivatives w.r.t natural coordinates
C
      CALL SFDER9(NNE,R,S,P)
C
C     Computes the Jacobian operator and its determinant at point (r,s)
C
      CALL JACOPER(NNE,XJ,XX,P)
      DDET=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1)
C
C     Computes the inverse of the Jacobiam operator
C
      CALL JACINVE(XJ,DDET,XJI)
C
C     Jacobian Inverse times Shape Functions derivatives w.r.t natural coordinates
C     to produce shape function derivatives w.r.t x,y coordinates.
C          
      CALL MMULT(XJI,2,2,P,2,NNE,AUX1)
C
C     Assembles B matrix for a
C     Cosserat element.
C
      DO I=1,NNE
        II=2*(I-1)+1
        B(1,II)=AUX1(1,I)
        B(2,II+1)=AUX1(2,I)
        B(4,II)=AUX1(2,I)
        B(4,II+1)=AUX1(1,I)
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE STDM9SH:GENERATES THE STRAIN-DISPLACEMENT MATRIX B    C
C     AND JACOBIAN DETERMINANT DDET AT THE POINT r ,s                  C
C                                                 i  j                 C
C     FOR AN 9-NODED 2D ELEMENT- SH WAVE                               C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE STDM9SH(NNE,NDOFEL,NTENS,XX,MCRD,B,DDET,R,S,XBAR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ONE=1.D0,TWO=2.D0,THREE=3.D0,HALF=0.5D0)
C
      DIMENSION XX(MCRD,NNE),B(NTENS,NDOFEL),P(2,NNE),AUX1(2,NNE),
     1          XJ(2,2),XJI(2,2)
C
C     Initialize arrays
C
      XBAR=ONE
      CALL CLEAR(B,NTENS,NDOFEL)
      CALL CLEAR(XJ,2,2)
      CALL CLEAR(XJI,2,2)
      CALL CLEAR(P,2,NNE)
C
C     Shape functions derivatives w.r.t natural coordinates
C
      CALL SFDER9(NNE,R,S,P)
C
C     Computes the Jacobian operator and its determinant at point (r,s)
C
      CALL JACOPER(NNE,XJ,XX,P)
      DDET=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1)
C
C     Computes the inverse of the Jacobiam operator
C
      CALL JACINVE(XJ,DDET,XJI)
C
C     Jacobian Inverse times Shape Functions derivatives w.r.t natural coordinates
C     to produce shape function derivatives w.r.t x,y coordinates.
C
      CALL MMULT(XJI,2,2,P,2,NNE,AUX1)
C
C     Assembles B matrix for a
C     Cosserat element.
C
      DO I=1,NNE
        B(1,I)=AUX1(1,I)
        B(2,I)=AUX1(2,I)
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SHAPE FUNCTIONS AND THEIR DERIVATIVES                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE SFDER4:GENERATES THE SHAPE FUNCTION DERIVATIVES       C
C     ACCORDING TO ELEMENT TYPE AT THE POINT r ,s                      C
C                                             i  j                     C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SFDER4(NNE,R,S,P)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ONE=1.D0,TWO=2.D0, HALF=0.5D0,QUARTER=0.25D0,FOUR=4.D0)
C
      DIMENSION P(2,NNE)
C
      RP= ONE+R
      SP= ONE+S
      RM= ONE-R
      SM= ONE-S
C
C     4-NODED ELEMENT
C
C     Derivatives w.r.t the natural coordinates
C     w.r.t.r
C
      P(1,2)= QUARTER*SM
      P(1,1)=-QUARTER*SM
      P(1,4)=-QUARTER*SP
      P(1,3)= QUARTER*SP
C
C     w.r.t.s
C
      P(2,2)=-QUARTER*RP
      P(2,1)=-QUARTER*RM
      P(2,4)= QUARTER*RM
      P(2,3)= QUARTER*RP
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE SFDER6:GENERATES THE SHAPE FUNCTION DERIVATIVES       C
C     ACCORDING TO ELEMENT TYPE AT THE POINT r ,s                      C
C                                             i  j                     C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SFDER6(NNE,R,S,P)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ONE=1.D0,TWO=2.D0,THREE=3.0D0, HALF=0.5D0,
     1          QUARTER=0.25D0,FOUR=4.D0,EIGHT=8.0D0, ZERO=0.0D0)
C
      DIMENSION P(2,NNE)
C
C     6-NODED ELEMENT
C
C     Derivatives w.r.t the natural coordinates
C     w.r.t.r
C
      P(1,6)=-FOUR*S
      P(1,5)= FOUR*S
      P(1,4)= FOUR-EIGHT*R-FOUR*S
      P(1,3)= ZERO
      P(1,2)= FOUR*R-ONE
      P(1,1)=-THREE+FOUR*R+FOUR*S
C
C     w.r.t.s
C
      P(2,6)= FOUR-FOUR*R-EIGHT*S
      P(2,5)= FOUR*R
      P(2,4)=-FOUR*R
      P(2,3)= FOUR*S-ONE
      P(2,2)= ZERO
      P(2,1)=-THREE+FOUR*R+FOUR*S
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE SFDER8:GENERATES THE SHAPE FUNCTION DERIVATIVES       C
C     ACCORDING TO ELEMENT TYPE AT THE POINT r ,s                      C
C                                             i  j                     C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SFDER8(NNE,R,S,P)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ONE=1.D0,TWO=2.D0, HALF=0.5D0,QUARTER=0.25D0,FOUR=4.D0)
C
      DIMENSION P(2,NNE)
C
      RP= ONE+R
      SP= ONE+S
      RM= ONE-R
      SM= ONE-S
      RMS=ONE-R**TWO
      SMS=ONE-S**TWO
C
C     8-NODED ELEMENT
C
C     Derivatives w.r.t the natural coordinates
C     w.r.t.r
C      
      P(1,6)= HALF*SMS
      P(1,5)=-R*SM
      P(1,8)=-HALF*SMS
      P(1,7)=-R*SP
      P(1,2)= QUARTER*SM-HALF*P(1,5)-HALF*P(1,6)
      P(1,1)=-QUARTER*SM-HALF*P(1,8)-HALF*P(1,5)
      P(1,4)=-QUARTER*SP-HALF*P(1,7)-HALF*P(1,8)
      P(1,3)= QUARTER*SP-HALF*P(1,7)-HALF*P(1,6)
C
C     w.r.t.s
C
      P(2,6)=-S*RP
      P(2,5)=-HALF*RMS
      P(2,8)=-S*RM
      P(2,7)= HALF*RMS
      P(2,2)=-QUARTER*RP-HALF*P(2,5)-HALF*P(2,6)
      P(2,1)=-QUARTER*RM-HALF*P(2,8)-HALF*P(2,5)
      P(2,4)= QUARTER*RM-HALF*P(2,7)-HALF*P(2,8)
      P(2,3)= QUARTER*RP-HALF*P(2,7)-HALF*P(2,6)
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     SUBROUTINE SFDER9:GENERATES THE SHAPE FUNCTION DERIVATIVES       C
C     ACCORDING TO ELEMENT TYPE AT THE POINT r ,s                      C
C                                             i  j                     C
C     B    =STRAIN-DISPLACEMENT MATRIX                                 C
C     DDET =JACOBIAN DETERMINANT                                       C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SFDER9(NNE,R,S,P)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER(ONE=1.D0,TWO=2.D0, HALF=0.5D0,QUARTER=0.25D0,FOUR=4.D0)
C
      DIMENSION P(2,NNE)
C
      RP= ONE+R
      SP= ONE+S
      RM= ONE-R
      SM= ONE-S
      RMS=ONE-R**TWO
      SMS=ONE-S**TWO
C
C     9-NODED ELEMENT
C     Derivatives w.r.t the natural coordinates
C     w.r.t.r
C
      P(1,9)=-TWO*R*SMS
      P(1,8)=-HALF*SMS-HALF*P(1,9)
      P(1,7)=-R*SP-HALF*P(1,9)
      P(1,6)= HALF*SMS-HALF*P(1,9)
      P(1,5)=-R*SM-HALF*P(1,9)
      P(1,4)=-QUARTER*SP-HALF*P(1,7)-HALF*P(1,8)-QUARTER*P(1,9)
      P(1,3)= QUARTER*SP-HALF*P(1,7)-HALF*P(1,6)-QUARTER*P(1,9)
      P(1,2)= QUARTER*SM-HALF*P(1,5)-HALF*P(1,6)-QUARTER*P(1,9)
      P(1,1)=-QUARTER*SM-HALF*P(1,8)-HALF*P(1,5)-QUARTER*P(1,9)
C
C        w.r.t.s
C
      P(2,9)=-TWO*S*RMS
      P(2,8)=-S*RM-HALF*P(2,9)
      P(2,7)= HALF*RMS-HALF*P(2,9)
      P(2,6)=-S*RP-HALF*P(2,9)
      P(2,5)=-HALF*RMS-HALF*P(2,9)
      P(2,4)= QUARTER*RM-HALF*P(2,7)-HALF*P(2,8)-QUARTER*P(2,9)
      P(2,3)= QUARTER*RP-HALF*P(2,7)-HALF*P(2,6)-QUARTER*P(2,9)
      P(2,2)=-QUARTER*RP-HALF*P(2,5)-HALF*P(2,6)-QUARTER*P(2,9)
      P(2,1)=-QUARTER*RM-HALF*P(2,8)-HALF*P(2,5)-QUARTER*P(2,9)
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE SHAPEF6(SN,R,S)                                      C
C      6-noded triangular element shape function at point R,S          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SHAPEF6(SN,R,S)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ONE=1.0D0,HALF=0.5D0,FOUR=4.0D0, TWO=2.0D0)
C
      DIMENSION SN(6)
C
      CALL CLEARV(SN,6)
C
      OMRMS=ONE-R-S
C
      SN(6)=FOUR*S*OMRMS
      SN(5)=FOUR*R*S
      SN(4)=FOUR*R*OMRMS
      SN(3)=TWO*S*(S-HALF)
      SN(2)=TWO*R*(R-HALF)
      SN(1)=TWO*(HALF-R-S)*OMRMS
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE SHAPEF8(SN,R,S)                                      C
C      8-noded element shape function at point R,S                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SHAPEF8(SN,R,S)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ONE=1.0D0,HALF=0.5D0,QUART=0.25D0)
C
      DIMENSION SN(8)
C
      CALL CLEARV(SN,8)
C
      RP =ONE+R
      RM =ONE-R
      RMS=ONE-R*R
      SP =ONE+S
      SM =ONE-S
      SMS=ONE-S*S
C
      SN(8)=HALF*SMS*RM
      SN(7)=HALF*RMS*SP
      SN(6)=HALF*SMS*RP
      SN(5)=HALF*RMS*SM
      SN(1)=QUART*RM*SM-HALF*SN(8)-HALF*SN(5)
      SN(2)=QUART*RP*SM-HALF*SN(6)-HALF*SN(5)
      SN(3)=QUART*RP*SP-HALF*SN(6)-HALF*SN(7)
      SN(4)=QUART*RM*SP-HALF*SN(8)-HALF*SN(7)
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE SHAPEF9(SN,R,S)                                      C
C      9-noded element shape function at point R,S                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SHAPEF9(SN,R,S)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ONE=1.0D0,HALF=0.5D0,QUART=0.25D0)
C
      DIMENSION SN(9)
C
      CALL CLEARV(SN,9)
C
      RP =ONE+R
      RM =ONE-R
      RMS=ONE-R*R
      SP =ONE+S
      SM =ONE-S
      SMS=ONE-S*S
C
      SN(9)=RMS*SMS
      SN(8)=HALF*SMS*RM-HALF*SN(9)
      SN(7)=HALF*RMS*SP-HALF*SN(9)
      SN(6)=HALF*SMS*RP-HALF*SN(9)
      SN(5)=HALF*RMS*SM-HALF*SN(9)
      SN(1)=QUART*RM*SM-HALF*SN(8)-HALF*SN(5)-QUART*SN(9)
      SN(2)=QUART*RP*SM-HALF*SN(6)-HALF*SN(5)-QUART*SN(9)
      SN(3)=QUART*RP*SP-HALF*SN(6)-HALF*SN(7)-QUART*SN(9)
      SN(4)=QUART*RM*SP-HALF*SN(8)-HALF*SN(7)-QUART*SN(9)
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C                     MASS- MATRICES                                   C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE CMASS6(MASS,NDOFEL,RII,SII)                          C
C      Computes the mass matrix for an 6-noded triangular element      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE AMASS6(AMASS,NDOFEL,RII,SII)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (NTENS=2)
C
      DIMENSION SF(NTENS,NDOFEL),SFT(NDOFEL,NTENS),SN(6),
     1          AUX1(NDOFEL,NDOFEL),AMASS(NDOFEL,NDOFEL)
C
      CALL CLEAR(SF,NTENS,NDOFEL)
      CALL CLEAR(SFT,NDOFEL,NTENS)
      CALL CLEARV(SN,6)
      CALL CLEAR(AUX1,NDOFEL,NDOFEL)
      CALL CLEAR(AMASS,NDOFEL,NDOFEL)
C
      CALL SHAPEF6(SN,RII,SII)
C
      DO I=1,6
        KK=2*I-1
        SF(1,  KK)=SN(I)
        SF(2,KK+1)=SN(I)
      END DO
C
      CALL MTRAN(SF,NTENS,NDOFEL,SFT)
      CALL MMULT(SFT,NDOFEL,NTENS,SF,NTENS,NDOFEL,AUX1)
C
      DO I=1,NDOFEL
        DO J=1,NDOFEL
          AMASS(I,J)=DCMPLX(AUX1(I,J))
        END DO
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE CMASS8(MASS,NDOFEL,RII,SII)                          C
C      Computes the mass matrix for an 8-noded element                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE AMASS8(AMASS,NDOFEL,RII,SII)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (NTENS=2)
C
      DIMENSION SF(NTENS,NDOFEL),SFT(NDOFEL,NTENS),SN(8),
     1          AUX1(NDOFEL,NDOFEL),AMASS(NDOFEL,NDOFEL)
C
      CALL CLEAR(SF,NTENS,NDOFEL)
      CALL CLEAR(SFT,NDOFEL,NTENS)
      CALL CLEARV(SN,8)
      CALL CLEAR(AUX1,NDOFEL,NDOFEL)
      CALL CLEAR(AMASS,NDOFEL,NDOFEL)
C
      CALL SHAPEF8(SN,RII,SII)
C
      DO I=1,8
        KK=2*I-1
        SF(1,  KK)=SN(I)
        SF(2,KK+1)=SN(I)
      END DO
C
      CALL MTRAN(SF,NTENS,NDOFEL,SFT)
      CALL MMULT(SFT,NDOFEL,NTENS,SF,NTENS,NDOFEL,AUX1)
C
      DO I=1,NDOFEL
        DO J=1,NDOFEL
          AMASS(I,J)=AUX1(I,J)
        END DO
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      SUBROUTINE AMASS9(MASS,NDOFEL,RII,SII)                          C
C      Computes the mass matrix for an 8-noded element                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE AMASS9(AMASS,NDOFEL,RII,SII)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (NTENS=2)
C
      DIMENSION SF(NTENS,NDOFEL),SFT(NDOFEL,NTENS),SN(9),
     1          AUX1(NDOFEL,NDOFEL),AMASS(NDOFEL,NDOFEL)
C
      CALL CLEAR(SF,NTENS,NDOFEL)
      CALL CLEAR(SFT,NDOFEL,NTENS)
      CALL CLEARV(SN,9)
      CALL CLEAR(AUX1,NDOFEL,NDOFEL)
      CALL CLEAR(AMASS,NDOFEL,NDOFEL)
C
      CALL SHAPEF9(SN,RII,SII)
C
      DO I=1,9
        KK=2*I-1
        SF(1,  KK)=SN(I)
        SF(2,KK+1)=SN(I)
      END DO
C
      CALL MTRAN(SF,NTENS,NDOFEL,SFT)
      CALL MMULT(SFT,NDOFEL,NTENS,SF,NTENS,NDOFEL,AUX1)
C
      DO I=1,NDOFEL
        DO J=1,NDOFEL
          AMASS(I,J)=AUX1(I,J)
        END DO
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE KJACOPER                                                C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE JACOPER(NNE,XJA,XCORD,RDER)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(ZERO=0.0D0)
C
      DIMENSION XJA(2,2),XCORD(2,NNE),RDER(2,NNE)
C
      CALL CLEAR(XJA,2,2)
C
      DUM=ZERO
      DO K=1,2
        DO J=1,2
          DO I=1,NNE
            DUM=DUM+RDER(K,I)*XCORD(J,I)
          END DO
          XJA(K,J)=DUM
          DUM=ZERO
        END DO
      END DO
C
      RETURN
C
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C                                                                      C
C      2X2 Jacobian inverse                                            C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE JACINVE(XJA,DD,XJAI)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION XJA(2,2),XJAI(2,2),XADJ(2,2)
C
      XADJ(1,1)=XJA(2,2)
      XADJ(1,2)=XJA(1,2)
      XADJ(2,1)=XJA(2,1)
      XADJ(2,2)=XJA(1,1)
      DO J=1,2
        DO K=1,2
          COFA=((-1)**(J+K))*XADJ(J,K)
          XJAI(J,K)=COFA/DD
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
C   SUBROUTINE SHAPE1D                                                 C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SHAPE1D(SM,ETA)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ONE=1.0D0,HALF=0.5D0)
C
      DIMENSION SM(6,2)
C
      CALL CLEAR(SM,6,2)
C
      EP=ONE+ETA
      EM=ONE-ETA
      EMS=ONE-ETA**2
C
      S11=HALF*EM-HALF*EMS
      S51=HALF*EP-HALF*EMS
      S31=EMS
C
      SM(1,1)=S11
      SM(2,2)=SM(1,1)
      SM(3,1)=S31
      SM(4,2)=SM(3,1)
      SM(5,1)=S51
      SM(6,2)=SM(5,1)
C
      RETURN
C
      END
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE SHAPE1D                                                 C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SHAPE1DM(SM,ETA)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ONE=1.0D0,HALF=0.5D0)
C
      DIMENSION SM(6,2)
C
      CALL CLEAR(SM,6,2)
C
      EP=ONE+ETA
      EM=ONE-ETA
      EMS=ONE-ETA**2
C
      S1=HALF*EM-HALF*EMS
      S2=HALF*EP-HALF*EMS
      S3=EMS
C
      SM(1,1)=S1
      SM(2,2)=SM(1,1)
      SM(3,1)=S2
      SM(4,2)=SM(3,1)
      SM(5,1)=S3
      SM(6,2)=SM(5,1)
C
      RETURN
C
      END
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE SHAPE1D                                                 C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SHAPEF31D(SM,ETA)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (ONE=1.0D0,HALF=0.5D0)
C
      DIMENSION SM(3)
C
      CALL CLEARV(SM,3)
C
      EP=ONE+ETA
      EM=ONE-ETA
      EMS=ONE-ETA**2
C
      SM(1)=HALF*EM-HALF*EMS
      SM(2)=HALF*EP-HALF*EMS
      SM(3)=EMS
C
      RETURN
C
      END
C
C23456789012345678901234567890123456789012345678901234567890123456789012
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C   SUBROUTINE SHAPEMAT8                                               C
C                                                                      C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      SUBROUTINE SHAPEMAT8(SNN,SN)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION SNN(2,16),SN(8)
C
      DO II=1,8
        SNN(1,2*II-1)=SN(II)
        SNN(1,2*II  )=0.0
        SNN(2,2*II-1)=0.0
        SNN(2,2*II  )=SN(II)
      END DO
C
      RETURN
C
      END
C
cc      END MODULE FEMUTIL

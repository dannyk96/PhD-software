C-----------------------------------------------------------------------
      SUBROUTINE VMPL (E,V,STRESS,PL)
C
C      This forms the plastic matrix for a Von-Mises material in 2D
C
      REAL STRESS(*),TERM(4),PL(4,4)
      SX = STRESS(1)
      SY = STRESS(2)
      TXY= STRESS(3)
      SZ = STRESS(4)
      DSBAR = SQRT((SX-SY)**2+(SY-SZ)**2+(SZ-SX)**2+6.*TXY**2)/SQRT(2.)
      EE = 1.5*E/((1.+V)*DSBAR*DSBAR)
      TERM(1) = (2.*SX-SY-SZ)/3.
      TERM(2) = (2.*SY-SZ-SX)/3.
      TERM(3) = TXY
      TERM(4) = (2.*SZ-SX-SY)/3.
      DO I = 1,4
        DO J = I,4
          PL(I,J) = TERM(I)*TERM(J)*EE
          PL(J,I) = PL(I,J)   !-- symmetry
        ENDDO
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE MOCOPL (PHI,PSI,E,V,STRESS,PL,IPL)
C
C      This forms the plastic stress/strain matrix        -- 2D
C      for a Mohr-Coulomb material  (phi,psi in degrees)
C
C     hmm 3D is fairly easy as 'uncoupled' :
C     just abstract more of the stresses and loop more of ROW and COL
C
      REAL STRESS(4), ROW(6), COL(6), PL(IPL.IPL)

      SX = STRESS(1)
      SY = STRESS(2)
      TXY = STRESS(3)
      SZ = STRESS(4)

      PI = 4.*ATAN(1.)
      PHIR = PHI*PI/180.
      PSIR = PSI*PI/180.
      SNPH = SIN(PHIR)
      SNPS = SIN(PSIR)
      SQ3 = SQRT(3.)
      CC = 1.-2.*V
      DX = (2.*SX-SY-SZ)/3.
      DY = (2.*SY-SZ-SX)/3.
      DZ = (2.*SZ-SX-SY)/3.
      D2 = SQRT(-DX*DY-DY*DZ-DZ*DX+TXY*TXY)
      D3 = DX*DY*DZ-DZ*TXY*TXY
      TH = MIN(MAX(-1.,-3.*SQ3*D3/(2.*D2**3),1.)

C      IF (TH.GT. 1.) TH = 1.
C      IF (TH.LT.-1.) TH = -1.
C      TH = ASIN(TH)/3. 

      SNTH = SIN(ASIN(TH)/3.)
      IF (ABS(SNTH).GT..49) THEN     !-- close to the 'Kink'
        SIG = -1.
        IF (SNTH.LT.0.) SIG = 1.
        RPH = SNPH*(1.+V)/3.
        RPS = SNPS*(1.+V)/3.
        CPS = .25*SQ3/D2*(1.+SIG*SNPS/3.)
        CPH = .25*SQ3/D2*(1.+SIG*SNPH/3.)
        COL(1) = RPH+CPH*((1.-V)*DX+V*(DY+DZ))
        COL(2) = RPH+CPH*((1.-V)*DY+V*(DZ+DX))
        COL(3) = CPH*CC*TXY
        COL(4) = RPH+CPH*((1.-V)*DZ+V*(DX+DY))
        ROW(1) = RPS+CPS*((1.-V)*DX+V*(DY+DZ))
        ROW(2) = RPS+CPS*((1.-V)*DY+V*(DZ+DX))
        ROW(3) = CPS*CC*TXY
        ROW(4) = RPS+CPS*((1.-V)*DZ+V*(DX+DY))
        EE = E/((1.+V)*CC*(RPH*SNPS+2.*CPH*CPS*D2*D2*CC))
      ELSE                       !------ 'normal' state of affairs
        ALP = ATAN(ABS((SX-SY)/(2.*TXY)))
        CA = COS(ALP)
        SA = SIN(ALP)
        DD = CC*SA
        S1 = 1.
        S2 = 1.
        IF ((SX-SY).LT.0.) S1 = -1.
        IF ( TXY.   LT.0.) S2 = -1.
        COL(1) = SNPH+S1*DD
        COL(2) = SNPH-S1*DD
        COL(3) = S2*CC*CA
        COL(4) = 2.*V*SNPH
        ROW(1) = SNPS+S1*DD
        ROW(2) = SNPS-S1*DD
        ROW(3) = S2*CC*CA
        ROW(4) = 2.*V*SNPS
        EE = E/(2.*(1.+V)*CC*(SNPH*SNPS+CC))
      END IF
      DO I = 1,4
        DO J = 1,4
          PL(I,J) = EE*ROW(I)*COL(J)
        ENDDO
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE MOCOUF (PHI,C,SIGM,DSBAR,THETA,F)
C
C     This calculates the value of the yield function
C     for a mohr-coulomb material (phi in degrees)
C
      PHIR = PHI*4.*ATAN(1.)/180.
      SNPH = SIN(PHIR)
      CSPH = COS(PHIR)
      CSTH = COS(THETA)
      SNTH = SIN(THETA)
      F = SNPH*SIGM+DSBAR* (CSTH/SQRT(3.)-SNTH*SNPH/3.) -C*CSPH
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE MOCOUQ (PSI,DSBAR,THETA,DQ1,DQ2,DQ3)
C
C      THIS SUBROUTINE FORMS THE DERIVATIVES OF A MOHR-COULOMB
C      POTENTIAL FUNCTION WITH RESPECT TO THE THREE INVARIANTS
C      PSI IN DEGREES
C
      PSIR = PSI*4.*ATAN(1.)/180.
      SNTH = SIN(THETA)
      SNPS = SIN(PSIR)
      SQ3 = SQRT(3.)
      DQ1 = SNPS
      IF(ABS(SNTH).GT..49)THEN
      C1 = 1.
      IF(SNTH.LT.0.)C1 = -1.
      DQ2 = (SQ3*.5-C1*SNPS*.5/SQ3)*SQ3*.5/DSBAR
      DQ3 = 0.
      ELSE
      CSTH = COS(THETA)
      CS3TH = COS(3.*THETA)
      TN3TH = TAN(3.*THETA)
      TNTH = SNTH/CSTH
      DQ2 = SQ3*CSTH/DSBAR*((1.+TNTH*TN3TH)+SNPS*(TN3TH-TNTH)/SQ3)*.5
      DQ3 = 1.5*(SQ3*SNTH+SNPS*CSTH)/(CS3TH*DSBAR*DSBAR)
      END IF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FORMM (STRESS,M1,M2,M3)
C
C     This forms the derivatives of the invariants   - 2D
C     with respect to the stresses
C
      REAL STRESS(*),M1(4,4),M2(4,4),M3(4,4)

      CALL NULL (M1,4,4,4)
      CALL NULL (M2,4,4,4)
      CALL NULL (M3,4,4,4)

      SX = STRESS(1)
      SY = STRESS(2)
      TXY= STRESS(3)
      SZ = STRESS(4)
      DX = (2.*SX-SY-SZ)/3.
      DY = (2.*SY-SZ-SX)/3.
      DZ = (2.*SZ-SX-SY)/3.
      SIGM = ( SX+SY+SZ)/3.

      M1(1,1) = 1.       !-- hmm could really set all to 1.
      M1(1,2) = 1.       !-- then null the 3rd row & col
      M1(2,1) = 1.
      M1(1,4) = 1.
      M1(4,1) = 1.
      M1(2,2) = 1.
      M1(2,4) = 1.
      M1(4,2) = 1.       !- note: only 7/16 zeros 
      M1(4,4) = 1.
      DO I = 1,4
        DO J = 1,4
          M1(I,J) = M1(I,J)/(9.*SIGM)
        ENDDO
      ENDDO
      M2(1,1) = .6666666666666666     !- diagonal    (FORMM3 does /3.
      M2(2,2) = .6666666666666666                 !   at the end :-)
      M2(4,4) = .6666666666666666
      M2(3,3) = 2.                    !- shear term
      M2(2,4) = -.3333333333333333
      M2(4,2) = -.3333333333333333    !- off-diagonals
      M2(1,2) = -.3333333333333333
      M2(2,1) = -.3333333333333333
      M2(1,4) = -.3333333333333333    !- note: only 6/16 zeros !
      M2(4,1) = -.3333333333333333

      M3(1,1) = DX/3.
      M3(2,4) = DX/3.
      M3(4,2) = DX/3.
      M3(2,2) = DY/3.
      M3(1,4) = DY/3.
      M3(4,1) = DY/3.
      M3(4,4) = DZ/3.
      M3(1,2) = DZ/3.
      M3(2,1) = DZ/3.
      M3(3,3) =-DZ
      M3(3,4) = TXY*-2./3.
      M3(4,3) = TXY*-2./3.
      M3(1,3) = TXY/3.
      M3(3,1) = TXY/3.
      M3(2,3) = TXY/3.                 !- note: only 6/16 zeros
      M3(3,2) = TXY/3.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FORMM3 (STRESS,M1,M2,M3)
C
C      This forms the derivatives of the invariants
C      with respect to the stresses (3-d)
C
C  ---> OK for 2d/3d... replace '3' with NODOF (careful with axisym)
C
C
      REAL STRESS(*), M1(6,6), M2(6,6), M3(6,6)

      CALL NULL (M1,6,6,6)
      CALL NULL (M2,6,6,6)
C     CALL NULL (M3,6,6,6)

      SX  = STRESS(1)
      SY  = STRESS(2)
      SZ  = STRESS(3)
      TXY = STRESS(4)
      TYZ = STRESS(5)
      TZX = STRESS(6)
      SIGM= (SX+SY+SZ)/3.
      DX = SX-SIGM         !--- hmm only the deviatorics are relevant
      DY = SY-SIGM
      DZ = SZ-SIGM

      DO I = 1,3
        DO J = 1,3         !-- only J=I:3 is needed (symmetry)
          M1(I,J) = 1./(3.*SIGM)   !- direct-stres terms
        ENDDO
      ENDDO
      DO I = 1,3
        M2(I  ,I  ) = 2.           !- direct-stress terms
        M2(I+3,I+3) = 6.           !- shear -stress terms
      ENDDO
      M2(1,2) = -1.                !- off-diagonal terms
      M2(1,3) = -1.
      M2(2,3) = -1.

      M3(1,1) = DX                 !- loop ?
      M3(1,2) = DZ
      M3(1,3) = DY
      M3(1,4) = TXY
      M3(1,5) = -2.*TYZ
      M3(1,6) = TZX                !- or use IF round the 3d bits
      M3(2,2) = DY
      M3(2,3) = DX
      M3(2,4) = TXY
      M3(2,5) = TYZ
      M3(2,6) = -2.*TZX
      M3(3,3) = DZ
      M3(3,4) = -2.*TXY
      M3(3,5) = TYZ
      M3(3,6) = TZX
      M3(4,4) = -3.*DZ
      M3(4,5) =  3.*TZX
      M3(4,6) =  3.*TYZ
      M3(5,5) = -3.*DX
      M3(5,6) =  3.*TXY
      M3(6,6) = -3.*DY
      DO I = 1,6
        DO J = I,6
          M1(I,J) = M1(I,J)/3.
          M1(J,I) = M1(I,J)       !-- symmetry
          M2(I,J) = M2(I,J)/3.    !-- and scaling 
          M2(J,I) = M2(I,J)       
          M3(I,J) = M3(I,J)/3.
          M3(J,I) = M3(I,J)       
        ENDDO
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------


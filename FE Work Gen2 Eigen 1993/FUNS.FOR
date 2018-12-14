C-----------------------------------------------------------------------
C        - These are ALL the Shape Function subroutines
C
C   - The names have all been made consistent! eg. FM8NQ   :-)
C   - Also the Sampling-points have been made consistent (+ no weights)
C   - About half of the subroutines were also in FE5LIB.. the others
C     (eg. 1D elements, 12N-Quad,10N-triangle and the 14 N-brick are 
C      all by Dan Kidger )
C
C   - Also included are 2 general-purpose routines: 
C      SHAPE  which calls up the appropriate element subroutine and
C      WTHATN which calls SHAPE to get the local-coords of an element
C
C-----------------------------------------------------------------------
      SUBROUTINE SHAPE (NODOF,NOD,ITYPE,DER,IDER,FUN,SMP,ISMP,I)
C
C     A macro header routine to call up the shape functions
C     and derivatives of ANY known element.
C     ITYPE = optional element type
C                                       Dan Kidger   April '91
C
C   (.. samp order reversed for 4nq,8nq,8nb,11nb etc. )
C    so that the 'first' coord is ALWAYS local-X !         8-4-92
C
      REAL SMP(ISMP,*), DER(IDER,*), FUN(*), SAMP

      FUN(1) = 230964.
      IF (NODOF.EQ.1)THEN
        IF(NOD.EQ. 2) CALL FM2NL  (DER,IDER,FUN,SMP,ISMP,I)
        IF(NOD.EQ. 3) CALL FM3NL  (DER,IDER,FUN,SMP,ISMP,I)
      ELSEIF (NODOF.EQ.2)THEN
        IF(NOD.EQ. 3) CALL FM3NT  (DER,IDER,FUN,SMP,ISMP,I)
        IF(NOD.EQ. 6) CALL FM6NT  (DER,IDER,FUN,SMP,ISMP,I)
        IF(NOD.EQ.10) CALL FM10NT (DER,IDER,FUN,SMP,ISMP,I)
        IF(NOD.EQ.15) CALL FM15NT (DER,IDER,FUN,SMP,ISMP,I)

        IF(NOD.EQ. 4) CALL FM4NQ  (DER,IDER,FUN,SMP,ISMP,I)
        IF(NOD.EQ. 5) CALL FM5NQ  (DER,IDER,FUN,SMP,ISMP,I)
        IF(NOD.EQ. 8) CALL FM8NQ  (DER,IDER,FUN,SMP,ISMP,I)
        IF(NOD.EQ. 9) CALL FM9NQ  (DER,IDER,FUN,SMP,ISMP,I)
        IF(NOD.EQ.12) CALL FM12NQ (DER,IDER,FUN,SMP,ISMP,I) 
        IF(NOD.EQ.17) CALL FM17NQ (DER,IDER,FUN,SMP,ISMP,I) 
 
      ELSEIF (NODOF.EQ.3)THEN
        IF(NOD.EQ. 4) CALL FM4NP  (DER,IDER,FUN,SMP, ISMP,I)

        IF(NOD.EQ. 8) CALL FM8NB  (DER,IDER,FUN,SMP, ISMP,I)
        IF(NOD.EQ.11) CALL FM11NB (DER,IDER,FUN,SMP, ISMP,I)
        IF(NOD.EQ.14) CALL FM14NB (DER,IDER,FUN,SMP, ISMP,I)
        IF(NOD.EQ.20) CALL FM20NB (DER,IDER,FUN,SMP, ISMP,I)
        IF(NOD.EQ.26) CALL FM26NB (DER,IDER,FUN,SMP, ISMP,I)
        IF(NOD.EQ.27) CALL FM27NB (DER,IDER,FUN,SMP, ISMP,I)
      ENDIF
      IF (ABS(FUN(1)-230964.).LT..001)
     +   WRITE(*,'(A,I2,A,I3,A,I3,a,i2)')
     + '** WARNING: element not found: nodof=',nodof,' nod=',nod
     +,' type=',itype
      RETURN
      END
C***********************************************************************
c-----------------------------------------------------------------------
      SUBROUTINE WTHATN (NEN,NDIME,ITYPE,LCOORD,ILCOORD)
C
C     This recursively searches coords in the range -1 to +1 searching
C     for the local coords of an element usng SHAPE  - Dan Kidger 1992
C
      PARAMETER (MAX_LEVEL=12,IDER=5,ISMP=1, TOL=1.E-5 )
      REAL LCOORD(ILCOORD,*)
      REAL FUN(30), DER(IDER,30), SMP(ISMP,IDER)
      INTEGER POINT(5)
      DO I=1,NEN             !- null FUN incase SHAPE can't find it :-)
        FUN(I) = 0.
      ENDDO
      INODE = 0
      DO N=1,MAX_LEVEL       !- loop the 'recursion' level 
        DO J=1,NDIME   
          POINT(J) = 0       ! set-up the 'initial values'
        ENDDO
c-------------- loop through the local co-ords -------------------------
      DO LEVEL=1,32760
        IF (LEVEL.GT.1) THEN
c------------------ 'toggle' up the POINT pointers ---------------------
    2     IPT = NDIME
    1     POINT(IPT) = POINT(IPT) + 1
          IF (POINT(IPT).GT.N)THEN
            POINT(IPT) = 0
            IPT = IPT - 1              
            IF (IPT.GT.0) GOTO 1        !- push left
            GOTO 99                     !- finished all possiblilities
          ENDIF
        ENDIF
c--------- Sieve of Erasmthus (sp?) search of the prime factors --------
      IF (N.GE.2) THEN
        DO J=2,MAX_LEVEL/2
          IF (MOD(N,J).EQ.0)THEN
            IFLAG=0
            DO IPT=1,NDIME
              IF (MOD(POINT(IPT),J).NE.0) IFLAG=1
            ENDDO  
            IF (IFLAG.EQ.0) GOTO 2      !- try another point
          ENDIF
        ENDDO
      ENDIF
c--------------- form the sampling point into SMP ----------------------
        DO J=1,NDIME
          SMP(1,J) = -1. + 2. * REAL(POINT(J)) / REAL(N) 
        ENDDO
        CALL SHAPE (NDIME,NEN,ITYPE,DER,IDER,FUN,SMP,ISMP,1)

C------------ loop the nodes to find which has FUN = 1.0 ---------------
        MAXPOS = 0
        NZERO = 0
        DO I=1,NEN
          IF (ABS(FUN(I)-1.).LT.TOL)THEN           !- 'own node'
            MAXPOS = I                        
          ELSEIF (ABS(FUN(I)).LT.TOL)THEN          !- 'other node'
            NZERO = NZERO +1
          ENDIF
        ENDDO

        IF (MAXPOS.NE.0.AND.NZERO+1.EQ.NEN) THEN    !- a 'hit'
          INODE = INODE + 1
          DO J=1,NDIME
            LCOORD(MAXPOS,J) = SMP(1,J)
          ENDDO
        ENDIF
        IF (INODE.EQ.NEN)  RETURN                   !- all done

      ENDDO        !- end of the 'pseudo' local co-ord permutation loop
   99 CONTINUE
      ENDDO        !- end of the recusion 'level'
      PRINT*,'*** ERROR: coords not found after',max_level,' iterations'
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FM2NL (DER,IDER,FUN,SAMP,ISAMP,I)
C
C     This forms the shape functions and
C     their derivatives for a 2-noded line element
C             by  Dan Kidger, Jan 1993
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      X=SAMP(I,1)
      FUN(1) = (1.-X)/2.
      FUN(2) = (1.+X)/2.
      DER(1,1) =  -1./2.
      DER(1,2) =   1./2.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FM3NL (DER,IDER,FUN,SAMP,ISAMP,I)
C
C     This forms the shape functions and
C     their derivatives for a 3-noded line element
C             by  Dan Kidger, Jan 1993
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      X=SAMP(I,1)
      FUN(1) = (1.-X)*X/2.
      FUN(2) = (1.+X)*(1.+X)/4.
      FUN(3) = (1.+X)*X/2.
      DER(1,1) =  -1./2.    !-
      DER(1,2) =   1./2.    !-  fix these !
      DER(1,3) =   1./2.    !-
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FM3NT (DER,IDER,FUN,SAMP,ISAMP,I)
C
C     This forms the shape functions and
C     their derivatives for a 3-noded triangular element
C
      REAL DER(IDER,*),SAMP(ISAMP,*),FUN(*)
      FUN(1)=SAMP(I,1)
      FUN(2)=SAMP(I,2)
      FUN(3)=1.-FUN(1)-FUN(2)
      DER(1,1)=1.
      DER(1,2)=0.
      DER(1,3)=-1.
      DER(2,1)=0.
      DER(2,2)=1.
      DER(2,3)=-1.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FM6NT(DER,IDER,FUN,SAMP,ISAMP,I)
C
C     This forms the shape functions and
C     their derivatives for 6-noded triangular elements
C
      REAL DER(IDER,*),SAMP(ISAMP,*),FUN(*)
      C1=SAMP(I,1)
      C2=SAMP(I,2)
      C3=1.-C1-C2
      FUN(1)=(2.*C1-1.)*C1
      FUN(2)=4.*C1*C2
      FUN(3)=(2.*C2-1.)*C2
      FUN(4)=4.*C2*C3
      FUN(5)=(2.*C3-1.)*C3
      FUN(6)=4.*C3*C1
      DER(1,1)=4.*C1-1.
      DER(1,2)=4.*C2
      DER(1,3)=0.
      DER(1,4)=-4.*C2
      DER(1,5)=-(4.*C3-1.)
      DER(1,6)=4.*(C3-C1)
      DER(2,1)=0.
      DER(2,2)=4.*C1
      DER(2,3)=4.*C2-1.
      DER(2,4)=4.*(C3-C2)
      DER(2,5)=-(4.*C3-1.)
      DER(2,6)=-4.*C1
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FM10NT (DER,IDER,FUN,SAMP,ISAMP,I)
C
C     This forms the shape functions and
C     their derivatives for 10-noded triangular elements
C          Created using REDUCE by  Dan Kidger  MAY 1990
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      X=SAMP(I,1)
      Y=SAMP(I,2)
      FUN(1)= ((3.*X-1.)*(3.*X-2.)*X)/2.
      FUN(2)=  (9.*(3.*X-1.)*X*Y)/2.
      FUN(3)=  (9.*(3.*Y-1.)*X*Y)/2.
      FUN(4)= ((3.*Y-1.)*(3.*Y-2.)*Y)/2.
      FUN(5)= -(9.*(X+Y-1.)*(3.*Y-1.)*Y)/2.
      FUN(6)=  (9.*(3.*X+3.*Y-2.)*(X+Y-1.)*Y)/2.
      FUN(7)=-((3.*X+3.*Y-1.)*(3.*X+3.*Y-2.)*(X+Y-1.))/2.
      FUN(8)=  (9.*(3.*X+3.*Y-2.)*(X+Y-1.)*X)/2.
      FUN(9)= -(9.*(3.*X-1.)*(X+Y-1.)*X)/2.
      FUN(10)=-27.*((Y-1.)+X)*X*Y
      DER(1,1)=(27.*X**2-18.*X+2.)/2.
      DER(1,2)=(9.*(6.*X-1.)*Y)/2.
      DER(1,3)=(9.*(3.*Y-1.)*Y)/2.
      DER(1,4)=0.
      DER(1,5)= -(9.*(3.*Y-1.)*Y)/2.
      DER(1,6)=  (9.*(6.*X+6.*Y-5.)*Y)/2.
      DER(1,7)=-(27.*X**2+54.*X*Y-36.*X+27.*Y**2-36.*Y+11.)/2.
      DER(1,8)=  (9.*(9.*X**2+12.*X*Y-10.*X+3.*Y**2-5.*Y+2.))/2.
      DER(1,9)= -(9.*(9.*X**2+6.*X*Y-8.*X-Y+1.))/2.
      DER(1,10)=-27.*(((Y-1.)+X)+X)*Y
      DER(2,1)=0.
      DER(2,2)=  (9.*(3.*X-1.)*X)/2.
      DER(2,3)=  (9.*(6.*Y-1.)*X)/2.
      DER(2,4)= (27.*Y**2-18.*Y+2.)/2.
      DER(2,5)= -(9.*((X+Y-1.)*(6.*Y-1.)+(3.*Y-1.)*Y))/2.
      DER(2,6)=  (9.*(3.*X**2+12.*X*Y-5.*X+9.*Y**2-10.*Y+2.))/2.
      DER(2,7)=-(27.*X**2+54.*X*Y-36.*X+27.*Y**2-36.*Y+11.)/2.
      DER(2,8)=  (9.*(6.*X+6.*Y-5.)*X)/2.
      DER(2,9)= -(9.*(3.*X-1.)*X)/2.
      DER(2,10)=-27.*(((Y-1.)+X)+Y)*X
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FM15NT (DER,IDER,FUN,SAMP,ISAMP,I)
C
C     This forms the shape functions and
C     their derivatives for a 15-noded triangular element
C
      REAL DER(IDER,*),SAMP(ISAMP,*),FUN(*)
      C1=SAMP(I,1)
      C2=SAMP(I,2)
      C3=1.-C1-C2
      T1=C1-.25
      T2=C1-.5
      T3=C1-.75
      T4=C2-.25
      T5=C2-.5
      T6=C2-.75
      T7=C3-.25
      T8=C3-.5
      T9=C3-.75
      FUN(1)=32./3.*C1*T1*T2*T3
      FUN(2)=128./3.*C1*C2*T1*T2
      FUN(3)=64.*C1*C2*T1*T4
      FUN(4)=128./3.*C1*C2*T4*T5
      FUN(5)=32./3.*C2*T4*T5*T6
      FUN(6)=128./3.*C2*C3*T4*T5
      FUN(7)=64.*C2*C3*T4*T7
      FUN(8)=128./3.*C2*C3*T7*T8
      FUN(9)=32./3.*C3*T7*T8*T9
      FUN(10)=128./3.*C3*C1*T7*T8
      FUN(11)=64.*C3*C1*T1*T7
      FUN(12)=128./3.*C3*C1*T1*T2
      FUN(13)=128.*C1*C2*T1*C3
      FUN(14)=128.*C1*C2*C3*T4
      FUN(15)=128.*C1*C2*C3*T7
      DER(1,1)=32./3.*(T2*T3*(T1+C1)+C1*T1*(T3+T2))
      DER(1,2)=128./3.*C2*(T2*(T1+C1)+C1*T1)
      DER(1,3)=64.*C2*T4*(T1+C1)
      DER(1,4)=128./3.*C2*T4*T5
      DER(1,5)=0.
      DER(1,6)=-128./3.*C2*T4*T5
      DER(1,7)=-64.*C2*T4*(T7+C3)
      DER(1,8)=-128./3.*C2*(T8*(T7+C3)+C3*T7)
      DER(1,9)=-32./3.*(T8*T9*(T7+C3)+C3*T7*(T8+T9))
      DER(1,10)=128./3.*(C3*T7*T8-C1*(T8*(T7+C3)+C3*T7))
      DER(1,11)=64.*(C3*T7*(T1+C1)-C1*T1*(T7+C3))
      DER(1,12)=128./3.*(C3*(T2*(T1+C1)+C1*T1)-C1*T1*T2)
      DER(1,13)=128.*C2*(C3*(T1+C1)-C1*T1)
      DER(1,14)=128.*C2*T4*(C3-C1)
      DER(1,15)=128.*C2*(C3*T7-C1*(T7+C3))
      DER(2,1)=0.0
      DER(2,2)=128./3.*C1*T1*T2
      DER(2,3)=64.*C1*T1*(T4+C2)
      DER(2,4)=128./3.*C1*(T5*(T4+C2)+C2*T4)
      DER(2,5)=32./3.*(T5*T6*(T4+C2)+C2*T4*(T6+T5))
      DER(2,6)=128./3.*((C3*(T5*(T4+C2)+C2*T4))-C2*T4*T5)
      DER(2,7)=64.*(C3*T7*(T4+C2)-C2*T4*(T7+C3))
      DER(2,8)=128./3.*(C3*T7*T8-C2*(T8*(T7+C3)+C3*T7))
      DER(2,9)=-32./3.*(T8*T9*(T7+C3)+C3*T7*(T8+T9))
      DER(2,10)=-128./3.*C1*(T8*(T7+C3)+C3*T7)
      DER(2,11)=-64.*C1*T1*(T7+C3)
      DER(2,12)=-128./3.*C1*T1*T2
      DER(2,13)=128.*C1*T1*(C3-C2)
      DER(2,14)=128.*C1*(C3*(T4+C2)-C2*T4)
      DER(2,15)=128.*C1*(C3*T7-C2*(C3+T7))
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FM4NQ (DER,IDER,FUN,SAMP,ISAMP,I)
C
C     This forms the shape functions and
C     their derivatives for the 4-noded quadrilateral element
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      ETA=SAMP(I,1)
      XI =SAMP(I,2)
      ETAM=.25*(1.-ETA)
      ETAP=.25*(1.+ETA)
      XIM=.25*(1.-XI)
      XIP=.25*(1.+XI)
      FUN(1)=4.*XIM*ETAM
      FUN(2)=4.*XIM*ETAP
      FUN(3)=4.*XIP*ETAP
      FUN(4)=4.*XIP*ETAM
      DER(1,1)=-ETAM
      DER(1,2)=-ETAP
      DER(1,3)=ETAP
      DER(1,4)=ETAM
      DER(2,1)=-XIM
      DER(2,2)=XIM
      DER(2,3)=XIP
      DER(2,4)=-XIP
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FM5NQ (DER,IDER,FUN,SAMP,ISAMP,I)
C
C     This forms the shape functions and
C     their derivatives for a 5-noded quadrilateral element
C         by  Dan Kidger         13th July 1992
C        .. the centre is a 'bubble function'
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      X = SAMP(I,1)
      Y = SAMP(I,2)
      FUN(1)=-((4.*X*Y+4.*X+4.*Y+3.)*(X-1.)*(Y-1.))/4.
      FUN(2)=-((4.*X*Y-4.*X+4.*Y-3.)*(X-1.)*(Y+1.))/4.
      FUN(3)=-((4.*X*Y-4.*X-4.*Y+3.)*(X+1.)*(Y+1.))/4.
      FUN(4)=-((4.*X*Y+4.*X-4.*Y-3.)*(X+1.)*(Y-1.))/4.
      FUN(5)=(X+1.)*(X-1.)*(Y+1.)*(Y-1.)
      DER(1,1)=-((8.*X*Y+8.*X-1.)*(Y-1.))/4.
      DER(1,2)=-((8.*X*Y-8.*X+1.)*(Y+1.))/4.
      DER(1,3)=-((8.*X*Y-8.*X-1.)*(Y+1.))/4.
      DER(1,4)=-((8.*X*Y+8.*X+1.)*(Y-1.))/4.
      DER(1,5)=2.*(Y+1.)*(Y-1.)*X
      DER(2,1)=-((8.*X*Y+8.*Y-1.)*(X-1.))/4.
      DER(2,2)=-((8.*X*Y+8.*Y+1.)*(X-1.))/4.
      DER(2,3)=-((8.*X*Y-8.*Y-1.)*(X+1.))/4.
      DER(2,4)=-((8.*X*Y-8.*Y+1.)*(X+1.))/4.
      DER(2,5)=2.*(X+1.)*(X-1.)*Y
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FM8NQ (DER,IDER,FUN,SAMP,ISAMP,I)
C
C     This forms the shape functions and
C     their derivatives for 8-noded quadrilateral elements
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      ETA=SAMP(I,1)
      XI =SAMP(I,2)
      ETAM=.25*(1.-ETA)
      ETAP=.25*(1.+ETA)
      XIM=.25*(1.-XI)
      XIP=.25*(1.+XI)
      FUN(1)=4.*ETAM*XIM*(-XI-ETA-1.)
      FUN(2)=32.*ETAM*XIM*ETAP
      FUN(3)=4.*ETAP*XIM*(-XI+ETA-1.)
      FUN(4)=32.*XIM*XIP*ETAP
      FUN(5)=4.*ETAP*XIP*(XI+ETA-1.)
      FUN(6)=32.*ETAP*XIP*ETAM
      FUN(7)=4.*XIP*ETAM*(XI-ETA-1.)
      FUN(8)=32.*XIM*XIP*ETAM
      DER(1,1)=ETAM*(2.*XI+ETA)
      DER(1,2)=-8.*ETAM*ETAP
      DER(1,3)=ETAP*(2.*XI-ETA)
      DER(1,4)=-4.*ETAP*XI
      DER(1,5)=ETAP*(2.*XI+ETA)
      DER(1,6)=8.*ETAP*ETAM
      DER(1,7)=ETAM*(2.*XI-ETA)
      DER(1,8)=-4.*ETAM*XI
      DER(2,1)=XIM*(XI+2.*ETA)
      DER(2,2)=-4.*XIM*ETA
      DER(2,3)=XIM*(2.*ETA-XI)
      DER(2,4)=8.*XIM*XIP
      DER(2,5)=XIP*(XI+2.*ETA)
      DER(2,6)=-4.*XIP*ETA
      DER(2,7)=XIP*(2.*ETA-XI)
      DER(2,8)=-8.*XIM*XIP
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FM9NQ (DER,IDER,FUN,SAMP,ISAMP,I)
C
C     This forms the shape functions and
C     their derivatives for 9-noded quadrilateral elements
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      ETA=SAMP(I,1)
      XI =SAMP(I,2)
      ETAM=ETA-1.
      XIM=XI-1.
      ETAP=ETA+1.
      XIP=XI+1.
      X2P1=2.*XI+1.
      X2M1=2.*XI-1.
      E2P1=2.*ETA+1.
      E2M1=2.*ETA-1.
      FUN(1)=.25*XI*XIM*ETA*ETAM
      FUN(2)=-.5*XI*XIM*ETAP*ETAM
      FUN(3)=.25*XI*XIM*ETA*ETAP
      FUN(4)=-.5*XIP*XIM*ETA*ETAP
      FUN(5)=.25*XI*XIP*ETA*ETAP
      FUN(6)=-.5*XI*XIP*ETAP*ETAM
      FUN(7)=.25*XI*XIP*ETA*ETAM
      FUN(8)=-.5*XIP*XIM*ETA*ETAM
      FUN(9)=XIP*XIM*ETAP*ETAM
      DER(1,1)=.25*X2M1*ETA*ETAM
      DER(1,2)=-.5*X2M1*ETAP*ETAM
      DER(1,3)=.25*X2M1*ETA*ETAP
      DER(1,4)=-XI*ETA*ETAP
      DER(1,5)=.25*X2P1*ETA*ETAP
      DER(1,6)=-.5*X2P1*ETAP*ETAM
      DER(1,7)=.25*X2P1*ETA*ETAM
      DER(1,8)=-XI*ETA*ETAM
      DER(1,9)=2.*XI*ETAP*ETAM
      DER(2,1)=.25*XI*XIM*E2M1
      DER(2,2)=-XI*XIM*ETA
      DER(2,3)=.25*XI*XIM*E2P1
      DER(2,4)=-.5*XIP*XIM*E2P1
      DER(2,5)=.25*XI*XIP*E2P1
      DER(2,6)=-XI*XIP*ETA
      DER(2,7)=.25*XI*XIP*E2M1
      DER(2,8)=-.5*XIP*XIM*E2M1
      DER(2,9)=2.*XIP*XIM*ETA
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FM12NQ (DER,IDER,FUN,SAMP,ISAMP,I)
C
C     This forms the shape functions and
C     their derivatives for 12-noded quadrilateral elements
C         by  Dan Kidger         8th April 1992
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      X = SAMP(I,1)
      Y = SAMP(I,2)
      FUN(1)= ((9.*X**2+9.*Y**2-10.)*(X-1.)*(Y-1.))/32.
      FUN(2)=-(9.*(X-1.)*(3.*Y-1.)*(Y+1.)*(Y-1.))/32.
      FUN(3)= (9.*(X-1.)*(3.*Y+1.)*(Y+1.)*(Y-1.))/32.
      FUN(4)=-((9.*X**2+9.*Y**2-10.)*(X-1.)*(Y+1.))/32.
      FUN(5)= (9.*(3.*X-1.)*(X+1.)*(X-1.)*(Y+1.))/32.
      FUN(6)=-(9.*(3.*X+1.)*(X+1.)*(X-1.)*(Y+1.))/32.
      FUN(7)= ((9.*X**2+9.*Y**2-10.)*(X+1.)*(Y+1.))/32.
      FUN(8)=-(9.*(X+1.)*(3.*Y+1.)*(Y+1.)*(Y-1.))/32.
      FUN(9)= (9.*(X+1.)*(3.*Y-1.)*(Y+1.)*(Y-1.))/32.
      FUN(10)=-((9.*X**2+9.*Y**2-10.)*(X+1.)*(Y-1.))/32.
      FUN(11)= (9.*(3.*X+1.)*(X+1.)*(X-1.)*(Y-1.))/32.
      FUN(12)=-(9.*(3.*X-1.)*(X+1.)*(X-1.)*(Y-1.))/32.
      DER(1,1)= ((27.*X**2-18.*X+9.*Y**2-10.)*(Y-1.))/32.
      DER(1,2)=-(9.*(3.*Y-1.)*(Y+1.)*(Y-1.))/32.
      DER(1,3)= (9.*(3.*Y+1.)*(Y+1.)*(Y-1.))/32.
      DER(1,4)=-((27.*X**2-18.*X+9.*Y**2-10.)*(Y+1.))/32.
      DER(1,5)= (9.*(9.*X**2-2.*X-3.)*(Y+1.))/32.
      DER(1,6)=-(9.*(9.*X**2+2.*X-3.)*(Y+1.))/32.
      DER(1,7)= ((27.*X**2+18.*X+9.*Y**2-10.)*(Y+1.))/32.
      DER(1,8)=-(9.*(3.*Y+1.)*(Y+1.)*(Y-1.))/32.
      DER(1,9)= (9.*(3.*Y-1.)*(Y+1.)*(Y-1.))/32.
      DER(1,10)=-((27.*X**2+18.*X+9.*Y**2-10.)*(Y-1.))/32.
      DER(1,11)= (9.*(9.*X**2+2.*X-3.)*(Y-1.))/32.
      DER(1,12)=-(9.*(9.*X**2-2.*X-3.)*(Y-1.))/32.
      DER(2,1)= ((9.*X**2+27.*Y**2-18.*Y-10.)*(X-1.))/32.
      DER(2,2)=-(9.*(X-1.)*(9.*Y**2-2.*Y-3.))/32.
      DER(2,3)= (9.*(X-1.)*(9.*Y**2+2.*Y-3.))/32.
      DER(2,4)=-((9.*X**2+27.*Y**2+18.*Y-10.)*(X-1.))/32.
      DER(2,5)= (9.*(3.*X-1.)*(X+1.)*(X-1.))/32.
      DER(2,6)=-(9.*(3.*X+1.)*(X+1.)*(X-1.))/32.
      DER(2,7)= ((9.*X**2+27.*Y**2+18.*Y-10.)*(X+1.))/32.
      DER(2,8)=-(9.*(X+1.)*(9.*Y**2+2.*Y-3.))/32.
      DER(2,9)= (9.*(X+1.)*(9.*Y**2-2.*Y-3.))/32.
      DER(2,10)=-((9.*X**2+27.*Y**2-18.*Y-10.)*(X+1.))/32.
      DER(2,11)= (9.*(3.*X+1.)*(X+1.)*(X-1.))/32.
      DER(2,12)=-(9.*(3.*X-1.)*(X+1.)*(X-1.))/32.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FM17NQ (DER,IDER,FUN,SAMP,ISAMP,I,J)
C
C     This forms the shape functions and
C     their derivatives for 17-noded quadrilateral elements
C         by  Dan Kidger      April 1992
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      X = SAMP(I,1)
      Y = SAMP(I,2)
      FUN(1)=-((4.*X**3-3.*X*Y-4.*X+4.*Y**3-4.*Y)*(X-1.)*(Y-1.))/12.
      FUN(2)= (2.*(X-1.)*(2.*Y-1.)*(Y+1.)*(Y-1.)*Y)/3.
      FUN(3)=-((X+4.*Y**2)*(X-1.)*(Y+1.)*(Y-1.))/2.
      FUN(4)= (2.*(X-1.)*(2.*Y+1.)*(Y+1.)*(Y-1.)*Y)/3.
      FUN(5)= ((4.*X**3+3.*X*Y-4.*X-4.*Y**3+4.*Y)*(X-1.)*(Y+1.))/12.
      FUN(6)=-(2.*(2.*X-1.)*(X+1.)*(X-1.)*(Y+1.)*X)/3.
      FUN(7)= ((4.*X**2-Y)*(X+1.)*(X-1.)*(Y+1.))/2.
      FUN(8)=-(2.*(2.*X+1.)*(X+1.)*(X-1.)*(Y+1.)*X)/3.
      FUN(9)= ((4.*X**3+3.*X*Y-4.*X+4.*Y**3-4.*Y)*(X+1.)*(Y+1.))/12.
      FUN(10)=-(2.*(X+1.)*(2.*Y+1.)*(Y+1.)*(Y-1.)*Y)/3.
      FUN(11)=-((X-4.*Y**2)*(X+1.)*(Y+1.)*(Y-1.))/2.
      FUN(12)=-(2.*(X+1.)*(2.*Y-1.)*(Y+1.)*(Y-1.)*Y)/3.
      FUN(13)=-((4.*X**3-3.*X*Y-4.*X-4.*Y**3+4.*Y)*(X+1.)*(Y-1.))/
     . 12.
      FUN(14)= (2.*(2.*X+1.)*(X+1.)*(X-1.)*(Y-1.)*X)/3.
      FUN(15)=-((4.*X**2+Y)*(X+1.)*(X-1.)*(Y-1.))/2.
      FUN(16)= (2.*(2.*X-1.)*(X+1.)*(X-1.)*(Y-1.)*X)/3.
      FUN(17)= (X+1.)*(X-1.)*(Y+1.)*(Y-1.)
      DER(1,1)=-((16.*X**3-12.*X**2-6.*X*Y-8.*X+4.*Y**3-Y+4.)*(Y-1.)
     . )/12.
      DER(1,2)= (2.*(2.*Y-1.)*(Y+1.)*(Y-1.)*Y)/3.
      DER(1,3)=-(((X+4.*Y**2)+(X-1.))*(Y+1.)*(Y-1.))/2.
      DER(1,4)= (2.*(2.*Y+1.)*(Y+1.)*(Y-1.)*Y)/3.
      DER(1,5)= ((16.*X**3-12.*X**2+6.*X*Y-8.*X-4.*Y**3+Y+4.)*(Y+1.))
     . /12.
      DER(1,6)=-(2.*(8.*X**3-3.*X**2-4.*X+1.)*(Y+1.))/3.
      DER(1,7)= (8.*X**2-Y-4.)*(Y+1.)*X
      DER(1,8)=-(2.*(8.*X**3+3.*X**2-4.*X-1.)*(Y+1.))/3.
      DER(1,9)= ((16.*X**3+12.*X**2+6.*X*Y-8.*X+4.*Y**3-Y-4.)*(Y+1.))
     . /12.
      DER(1,10)=-(2.*(2.*Y+1.)*(Y+1.)*(Y-1.)*Y)/3.
      DER(1,11)=-(((X-4.*Y**2)+(X+1.))*(Y+1.)*(Y-1.))/2.
      DER(1,12)=-(2.*(2.*Y-1.)*(Y+1.)*(Y-1.)*Y)/3.
      DER(1,13)=-((16.*X**3+12.*X**2-6.*X*Y-8.*X-4.*Y**3+Y-4.)*(Y-1.
     . ))/12.
      DER(1,14)= (2.*(8.*X**3+3.*X**2-4.*X-1.)*(Y-1.))/3.
      DER(1,15)=-(8.*X**2+Y-4.)*(Y-1.)*X
      DER(1,16)= (2.*(8.*X**3-3.*X**2-4.*X+1.)*(Y-1.))/3.
      DER(1,17)=2.*(Y+1.)*(Y-1.)*X
      DER(2,1)=-((4.*X**3-6.*X*Y-X+16.*Y**3-12.*Y**2-8.*Y+4.)*(X-1.)
     . )/12.
      DER(2,2)= (2.*(X-1.)*(8.*Y**3-3.*Y**2-4.*Y+1.))/3.
      DER(2,3)=-(X+8.*Y**2-4.)*(X-1.)*Y
      DER(2,4)= (2.*(X-1.)*(8.*Y**3+3.*Y**2-4.*Y-1.))/3.
      DER(2,5)= ((4.*X**3+6.*X*Y-X-16.*Y**3-12.*Y**2+8.*Y+4.)*(X-1.))
     . /12.
      DER(2,6)=-(2.*(2.*X-1.)*(X+1.)*(X-1.)*X)/3.
      DER(2,7)= (((4.*X**2-Y)-(Y+1.))*(X+1.)*(X-1.))/2.
      DER(2,8)=-(2.*(2.*X+1.)*(X+1.)*(X-1.)*X)/3.
      DER(2,9)= ((4.*X**3+6.*X*Y-X+16.*Y**3+12.*Y**2-8.*Y-4.)*(X+1.))
     . /12.
      DER(2,10)=-(2.*(X+1.)*(8.*Y**3+3.*Y**2-4.*Y-1.))/3.
      DER(2,11)=-(X-8.*Y**2+4.)*(X+1.)*Y
      DER(2,12)=-(2.*(X+1.)*(8.*Y**3-3.*Y**2-4.*Y+1.))/3.
      DER(2,13)=-((4.*X**3-6.*X*Y-X-16.*Y**3+12.*Y**2+8.*Y-4.)*(X+1.
     . ))/12.
      DER(2,14)= (2.*(2.*X+1.)*(X+1.)*(X-1.)*X)/3.
      DER(2,15)=-(((4.*X**2+Y)+(Y-1.))*(X+1.)*(X-1.))/2.
      DER(2,16)= (2.*(2.*X-1.)*(X+1.)*(X-1.)*X)/3.
      DER(2,17)=2.*(X+1.)*(X-1.)*Y
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE FM4NP (DER,IDER,FUN,SAMP,ISAMP,I)
C
C     This forms the shape functions and
C     their derivatives for 4-node tetrahedron elements
C
      REAL DER(IDER,*),SAMP(ISAMP,*),FUN(*)
      FUN(1)=SAMP(I,1)
      FUN(2)=SAMP(I,2)
      FUN(3)=SAMP(I,3)
      FUN(4)=1.-FUN(1)-FUN(2)-FUN(3)
      DO 1 M=1,3
      DO 1 N=1,4
    1 DER(M,N)=0.
      DER(1,1)=1.
      DER(2,2)=1.
      DER(3,3)=1.
      DER(1,4)=-1.
      DER(2,4)=-1.
      DER(3,4)=-1.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FM8NB (DER,IDER,FUN,SAMP,ISAMP,I)
C
C     This forms the shape functions and their
C     derivatives for 8-noded brick elements
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      ETA = SAMP(I,1)
      XI  = SAMP(I,2)
      ZETA= SAMP(I,3)
      ETAM=1.-ETA
      XIM=1.-XI
      ZETAM=1.-ZETA
      ETAP=ETA+1.
      XIP=XI+1.
      ZETAP=ZETA+1.
      FUN(1)=.125*XIM*ETAM*ZETAM
      FUN(2)=.125*XIM*ETAM*ZETAP
      FUN(3)=.125*XIP*ETAM*ZETAP
      FUN(4)=.125*XIP*ETAM*ZETAM
      FUN(5)=.125*XIM*ETAP*ZETAM
      FUN(6)=.125*XIM*ETAP*ZETAP
      FUN(7)=.125*XIP*ETAP*ZETAP
      FUN(8)=.125*XIP*ETAP*ZETAM
      DER(1,1)=-.125*ETAM*ZETAM
      DER(1,2)=-.125*ETAM*ZETAP
      DER(1,3)=.125*ETAM*ZETAP
      DER(1,4)=.125*ETAM*ZETAM
      DER(1,5)=-.125*ETAP*ZETAM
      DER(1,6)=-.125*ETAP*ZETAP
      DER(1,7)=.125*ETAP*ZETAP
      DER(1,8)=.125*ETAP*ZETAM
      DER(2,1)=-.125*XIM*ZETAM
      DER(2,2)=-.125*XIM*ZETAP
      DER(2,3)=-.125*XIP*ZETAP
      DER(2,4)=-.125*XIP*ZETAM
      DER(2,5)=.125*XIM*ZETAM
      DER(2,6)=.125*XIM*ZETAP
      DER(2,7)=.125*XIP*ZETAP
      DER(2,8)=.125*XIP*ZETAM
      DER(3,1)=-.125*XIM*ETAM
      DER(3,2)=.125*XIM*ETAM
      DER(3,3)=.125*XIP*ETAM
      DER(3,4)=-.125*XIP*ETAM
      DER(3,5)=-.125*XIM*ETAP
      DER(3,6)=.125*XIM*ETAP
      DER(3,7)=.125*XIP*ETAP
      DER(3,8)=-.125*XIP*ETAP
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FM11NB (DER,IDER,FUN,SAMP,ISAMP,I)
C
C     This forms the shape functions and their
C     derivatives for 11-noded brick elements
C        Created by Dan Kidger 1990
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)

      CALL FM8NB (DER,IDER,FUN,SAMP, ISAMP,I,J,K)
      ETA = SAMP(I,1)
      XI  = SAMP(I,2)
      ZETA= SAMP(I,3)
      FUN( 9) = 1.0-  XI*XI
      FUN(10) = 1.0- ETA*ETA
      FUN(11) = 1.0-ZETA*ZETA
      DO 1,I=1,3
      DO 1,J=9,11
    1 DER(I, J) =  0.0
      DER(1, 9) = -2.0*XI
      DER(2,10) = -2.0*ETA
      DER(3,11) = -2.0*ZETA
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FM14NB (DER,IDER,FUN,SAMP,ISAMP,I,J,K)
C
C     This forms the shape functions and their
c     derivatives for 14-noded brick elements
C        by Dan Kidger MAY 1989          Type A:   X*Y*Y etc terms
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      X=1+SAMP(I,1)
      Y=1+SAMP(J,1)
      Z=1+SAMP(K,1)
      FUN(1)=-(2.*((Z-4.)*Y+Z**2-4.*Z+7.)*X+(2.*Y-3.)*X**2
     . -2.*(4.*Z-7.)*Y+(2.*Z-3.)*Y**2-3.*Z**2+14.*Z-16.)/16.
      FUN(2)=(2.*((Z-4.)*Y-Z**2+2.*Z+1.)*X+(2.*Y-1.)*X**2-
     . (2.*Z-3.)*Y**2+2.*Y+3.*Z**2-6.*Z)/16.
      FUN(3)=(2.*(Y*Z+Z**2-4.*Z+1.)*X-(2.*Y-3.)*X**2-(Y-2.
     . )*(2.*Z-3.)*Y-Z**2+2.*Z)/16.
      FUN(4)=-(2.*(Y*Z-Z**2+2.*Z-1.)*X-(2.*Y-1.)*X**2+(Y-
     . 2.)*(2.*Z-3.)*Y+Z**2-2.*Z)/16.
      FUN(5)=(2.*((Z+2.)*Y-Z**2-3.)*X-(2.*Y-3.)*X**2-2.*(
     . 4.*Z-1.)*Y+(2.*Z-1.)*Y**2+3.*Z**2+2.*Z)/16.
      FUN(6)=-(2.*((Z+2.)*Y+Z**2-2.*Z-1.)*X-(2.*Y-1.)*X**2
     . -(2.*Z-1.)*Y**2-2.*Y-3.*Z**2+6.*Z)/16.
      FUN(7)=-(2.*((Z-2.)*Y-Z**2+3.)*X+(2.*Y-3.)*X**2-(Y-
     . 2.)*(2.*Z-1.)*Y+Z**2-2.*Z)/16.
      FUN(8)=(2.*((Z-2.)*Y+Z**2-2.*Z+1.)*X+(2.*Y-1.)*X**2+
     . (Y-2.)*(2.*Z-1.)*Y-Z**2+2.*Z)/16.
      FUN(9)=((Y-2.)*(2.*Z-3.)*Y-X**2+2.*X+Z**2-2.*Z)/4.
      FUN(10)=-((Y-2.)*(2.*Z-1.)*Y+X**2-2.*X-Z**2+2.*Z)/4.
      FUN(11)=((2.*Y-3.)*X**2-2.*(2.*Y-3.)*X+Y**2-2.*Y-Z**2+2.*Z)/4.
      FUN(12)=-((2.*Y-1.)*X**2-2.*(2.*Y-1.)*X-Y**2+2.*Y+Z**2-2.*Z)/4.
      FUN(13)=(2.*(Z**2-2.*Z-1.)*X+X**2-Y**2+2.*Y-3.*Z**2+6.*Z)/4.
      FUN(14)=-(2.*(Z-1.)**2*X-X**2+Y**2-2.*Y-Z**2+2.*Z)/4.
      DER(1,1)=-(((Z-4.)*Y+Z**2-4.*Z+7.)+(2.*Y-3.)*X)/8.
      DER(1,2)=(((Z-4.)*Y-Z**2+2.*Z+1.)+(2.*Y-1.)*X)/8.
      DER(1,3)=((Y*Z+Z**2-4.*Z+1.)-(2.*Y-3.)*X)/8.
      DER(1,4)=-((Y*Z-Z**2+2.*Z-1.)-(2.*Y-1.)*X)/8.
      DER(1,5)=(((Z+2.)*Y-Z**2-3.)-(2.*Y-3.)*X)/8.
      DER(1,6)=-(((Z+2.)*Y+Z**2-2.*Z-1.)-(2.*Y-1.)*X)/8.
      DER(1,7)=-(((Z-2.)*Y-Z**2+3.)+(2.*Y-3.)*X)/8.
      DER(1,8)=(((Z-2.)*Y+Z**2-2.*Z+1.)+(2.*Y-1.)*X)/8.
      DER(1,9) =-(X-1.)/2.
      DER(1,10)=-(X-1.)/2.
      DER(1,11)=((X-1.)*(2.*Y-3.))/2.
      DER(1,12)=-((X-1.)*(2.*Y-1.))/2.
      DER(1,13)=((Z**2-2.*Z-1.)+X)/2.
      DER(1,14)=-((Z-1.)**2-X)/2.
      DER(2,1)=((4.*Z-7.)-(2.*Z-3.)*Y-(Z-4.)*X-X**2)/8.
      DER(2,2)=-((2.*Z-3.)*Y-(Z-4.)*X-X**2-1.)/8.
      DER(2,3)=-(X**2-X*Z+2.*Y*Z-3.*Y-2.*Z+3.)/8.
      DER(2,4)=(X**2-X*Z-2.*Y*Z+3.*Y+2.*Z-3.)/8.
      DER(2,5)=-((4.*Z-1.)-(2.*Z-1.)*Y-(Z+2.)*X+X**2)/8.
      DER(2,6)=((2.*Z-1.)*Y-(Z+2.)*X+X**2+1.)/8.
      DER(2,7)=-(X**2+X*Z-2.*X-2.*Y*Z+Y+2.*Z-1.)/8.
      DER(2,8)=(X**2+X*Z-2.*X+2.*Y*Z-Y-2.*Z+1.)/8.
      DER(2,9)=((Y-1.)*(2.*Z-3.))/2.
      DER(2,10)=-((Y-1.)*(2.*Z-1.))/2.
      DER(2,11)=(X**2-2.*X+Y-1.)/2.
      DER(2,12)=-(X**2-2.*X-Y+1.)/2.
      DER(2,13)=-(Y-1.)/2.
      DER(2,14)=-(Y-1.)/2.
      DER(3,1)=-((Y+2.*Z-4.)*X+Y**2-4.*Y-3.*Z+7.)/8.
      DER(3,2)=((Y-2.*Z+2.)*X-Y**2+3.*Z-3.)/8.
      DER(3,3)=((Y+2.*Z-4.)*X-(Y-2.)*Y-Z+1.)/8.
      DER(3,4)=-((Y-2.*Z+2.)*X+(Y-2.)*Y+Z-1.)/8.
      DER(3,5)=((Y-2.*Z)*X+Y**2-4.*Y+3.*Z+1.)/8.
      DER(3,6)=-((Y+2.*Z-2.)*X-Y**2-3.*Z+3.)/8.
      DER(3,7)=-((Y-2.*Z)*X-(Y-2.)*Y+Z-1.)/8.
      DER(3,8)=((Y+2.*Z-2.)*X+(Y-2.)*Y-Z+1.)/8.
      DER(3,9)=((Y-2.)*Y+Z-1.)/2.
      DER(3,10)=-((Y-2.)*Y-Z+1.)/2.
      DER(3,11)=-(Z-1.)/2.
      DER(3,12)=-(Z-1.)/2.
      DER(3,13)=((2.*X-3.)*(Z-1.))/2.
      DER(3,14)=-((2.*X-1.)*(Z-1.))/2.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FM20NB (DER,IDER,FUN,SAMP,ISAMP,I)
C
C     This forms the shape functions and their
c     derivatives for 20-noded brick elements
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      INTEGER XII(20),ETAI(20),ZETAI(20)
      ETA = SAMP(I,1)
      XI  = SAMP(I,2)
      ZETA= SAMP(I,3)
      XII(1)=-1
      XII(2)=-1
      XII(3)=-1
      XII(9)=-1
      XII(10)=-1
      XII(13)=-1
      XII(14)=-1
      XII(15)=-1
      XII(4)=0
      XII(8)=0
      XII(16)=0
      XII(20)=0
      XII(5)=1
      XII(6)=1
      XII(7)=1
      XII(11)=1
      XII(12)=1
      XII(17)=1
      XII(18)=1
      XII(19)=1
      DO 1 L=1,8
    1 ETAI(L)=-1
      DO 2 L=9,12
    2 ETAI(L)=0
      DO 3 L=13,20
    3 ETAI(L)=1
      ZETAI(1)=-1
      ZETAI(7)=-1
      ZETAI(8)=-1
      ZETAI(9)=-1
      ZETAI(12)=-1
      ZETAI(13)=-1
      ZETAI(19)=-1
      ZETAI(20)=-1
      ZETAI(2)=0
      ZETAI(6)=0
      ZETAI(14)=0
      ZETAI(18)=0
      ZETAI(3)=1
      ZETAI(4)=1
      ZETAI(5)=1
      ZETAI(10)=1
      ZETAI(11)=1
      ZETAI(15)=1
      ZETAI(16)=1
      ZETAI(17)=1
      DO 4 L=1,20
      XIO=XI*XII(L)
      ETAO=ETA*ETAI(L)
      ZETAO=ZETA*ZETAI(L)
      IF(L.EQ.4.OR.L.EQ.8.OR.L.EQ.16.OR.L.EQ.20)THEN
      FUN(L)=.25*(1.-XI*XI)*(1.+ETAO)*(1.+ZETAO)
      DER(1,L)=-.5*XI*(1.+ETAO)*(1.+ZETAO)
      DER(2,L)=.25*ETAI(L)*(1.-XI*XI)*(1.+ZETAO)
      DER(3,L)=.25*ZETAI(L)*(1.-XI*XI)*(1.+ETAO)
      ELSE IF(L.GE.9.AND.L.LE.12)THEN
      FUN(L)=.25*(1.+XIO)*(1.-ETA*ETA)*(1.+ZETAO)
      DER(1,L)=.25*XII(L)*(1.-ETA*ETA)*(1.+ZETAO)
      DER(2,L)=-.5*ETA*(1.+XIO)*(1.+ZETAO)
      DER(3,L)=.25*ZETAI(L)*(1.+XIO)*(1.-ETA*ETA)
      ELSE IF(L.EQ.2.OR.L.EQ.6.OR.L.EQ.14.OR.L.EQ.18)THEN
      FUN(L)=.25*(1.+XIO)*(1.+ETAO)*(1.-ZETA*ZETA)
      DER(1,L)=.25*XII(L)*(1.+ETAO)*(1.-ZETA*ZETA)
      DER(2,L)=.25*ETAI(L)*(1.+XIO)*(1.-ZETA*ZETA)
      DER(3,L)=-.5*ZETA*(1.+XIO)*(1.+ETAO)
      ELSE
      FUN(L)=.125*(1.+XIO)*(1.+ETAO)*(1.+ZETAO)*(XIO+ETAO+ZETAO-2.)
      DER(1,L)=.125*XII(L)*(1.+ETAO)*(1.+ZETAO)*(2.*XIO+ETAO+ZETAO-1.)
      DER(2,L)=.125*ETAI(L)*(1.+XIO)*(1.+ZETAO)*(XIO+2.*ETAO+ZETAO-1.)
      DER(3,L)=.125*ZETAI(L)*(1.+XIO)*(1.+ETAO)*(XIO+ETAO+2.*ZETAO-1.)
      END IF
    4 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FM20NB_A (DER,IDER,FUN,SAMP,ISAMP,I)
C
C     This forms the shape functions and their
C     derivatives for 20-noded brick elements
C         Created using REDUCE by Dan Kidger,  March  1992
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      X = SAMP(I,1)
      Y = SAMP(I,2)
      Z = SAMP(I,3)
      FUN(1)= ((X+Y+Z+2.)*(X-1.)*(Y-1.)*(Z-1.))/8.
      FUN(2)=-((X-1.)*(Y+1.)*(Y-1.)*(Z-1.))/4.
      FUN(3)=-((X-Y+Z+2.)*(X-1.)*(Y+1.)*(Z-1.))/8.
      FUN(4)= ((X+1.)*(X-1.)*(Y+1.)*(Z-1.))/4.
      FUN(5)=-((X+Y-Z-2.)*(X+1.)*(Y+1.)*(Z-1.))/8.
      FUN(6)= ((X+1.)*(Y+1.)*(Y-1.)*(Z-1.))/4.
      FUN(7)= ((X-Y-Z-2.)*(X+1.)*(Y-1.)*(Z-1.))/8.
      FUN(8)=-((X+1.)*(X-1.)*(Y-1.)*(Z-1.))/4.
      FUN(9)=-((X-1.)*(Y-1.)*(Z+1.)*(Z-1.))/4.
      FUN(10)= ((X-1.)*(Y+1.)*(Z+1.)*(Z-1.))/4.
      FUN(11)=-((X+1.)*(Y+1.)*(Z+1.)*(Z-1.))/4.
      FUN(12)= ((X+1.)*(Y-1.)*(Z+1.)*(Z-1.))/4.
      FUN(13)=-((X+Y-Z+2.)*(X-1.)*(Y-1.)*(Z+1.))/8.
      FUN(14)= ((X-1.)*(Y+1.)*(Y-1.)*(Z+1.))/4.
      FUN(15)= ((X-Y-Z+2.)*(X-1.)*(Y+1.)*(Z+1.))/8.
      FUN(16)=-((X+1.)*(X-1.)*(Y+1.)*(Z+1.))/4.
      FUN(17)= ((X+Y+Z-2.)*(X+1.)*(Y+1.)*(Z+1.))/8.
      FUN(18)=-((X+1.)*(Y+1.)*(Y-1.)*(Z+1.))/4.
      FUN(19)=-((X-Y+Z-2.)*(X+1.)*(Y-1.)*(Z+1.))/8.
      FUN(20)= ((X+1.)*(X-1.)*(Y-1.)*(Z+1.))/4.
      DER(1,1)= ((2.*X+Y+Z+1.)*(Y-1.)*(Z-1.))/8.
      DER(1,2)=-((Y+1.)*(Y-1.)*(Z-1.))/4.
      DER(1,3)=-((2.*X-Y+Z+1.)*(Y+1.)*(Z-1.))/8.
      DER(1,4)= ((Y+1.)*(Z-1.)*X)/2.
      DER(1,5)=-((2.*X+Y-Z-1.)*(Y+1.)*(Z-1.))/8.
      DER(1,6)= ((Y+1.)*(Y-1.)*(Z-1.))/4.
      DER(1,7)= ((2.*X-Y-Z-1.)*(Y-1.)*(Z-1.))/8.
      DER(1,8)=-((Y-1.)*(Z-1.)*X)/2.
      DER(1,9)=-((Y-1.)*(Z+1.)*(Z-1.))/4.
      DER(1,10)= ((Y+1.)*(Z+1.)*(Z-1.))/4.
      DER(1,11)=-((Y+1.)*(Z+1.)*(Z-1.))/4.
      DER(1,12)= ((Y-1.)*(Z+1.)*(Z-1.))/4.
      DER(1,13)=-((2.*X+Y-Z+1.)*(Y-1.)*(Z+1.))/8.
      DER(1,14)= ((Y+1.)*(Y-1.)*(Z+1.))/4.
      DER(1,15)= ((2.*X-Y-Z+1.)*(Y+1.)*(Z+1.))/8.
      DER(1,16)=-((Y+1.)*(Z+1.)*X)/2.
      DER(1,17)= ((2.*X+Y+Z-1.)*(Y+1.)*(Z+1.))/8.
      DER(1,18)=-((Y+1.)*(Y-1.)*(Z+1.))/4.
      DER(1,19)=-((2.*X-Y+Z-1.)*(Y-1.)*(Z+1.))/8.
      DER(1,20)= ((Y-1.)*(Z+1.)*X)/2.
      DER(2,1)= ((X+2.*Y+Z+1.)*(X-1.)*(Z-1.))/8.
      DER(2,2)=-((X-1.)*(Z-1.)*Y)/2.
      DER(2,3)=-((X-2.*Y+Z+1.)*(X-1.)*(Z-1.))/8.
      DER(2,4)= ((X+1.)*(X-1.)*(Z-1.))/4.
      DER(2,5)=-((X+2.*Y-Z-1.)*(X+1.)*(Z-1.))/8.
      DER(2,6)= ((X+1.)*(Z-1.)*Y)/2.
      DER(2,7)= ((X-2.*Y-Z-1.)*(X+1.)*(Z-1.))/8.
      DER(2,8)=-((X+1.)*(X-1.)*(Z-1.))/4.
      DER(2,9)=-((X-1.)*(Z+1.)*(Z-1.))/4.
      DER(2,10)= ((X-1.)*(Z+1.)*(Z-1.))/4.
      DER(2,11)=-((X+1.)*(Z+1.)*(Z-1.))/4.
      DER(2,12)= ((X+1.)*(Z+1.)*(Z-1.))/4.
      DER(2,13)=-((X+2.*Y-Z+1.)*(X-1.)*(Z+1.))/8.
      DER(2,14)= ((X-1.)*(Z+1.)*Y)/2.
      DER(2,15)= ((X-2.*Y-Z+1.)*(X-1.)*(Z+1.))/8.
      DER(2,16)=-((X+1.)*(X-1.)*(Z+1.))/4.
      DER(2,17)= ((X+2.*Y+Z-1.)*(X+1.)*(Z+1.))/8.
      DER(2,18)=-((X+1.)*(Z+1.)*Y)/2.
      DER(2,19)=-((X-2.*Y+Z-1.)*(X+1.)*(Z+1.))/8.
      DER(2,20)= ((X+1.)*(X-1.)*(Z+1.))/4.
      DER(3,1)= ((X+Y+2.*Z+1.)*(X-1.)*(Y-1.))/8.
      DER(3,2)=-((X-1.)*(Y+1.)*(Y-1.))/4.
      DER(3,3)=-((X-Y+2.*Z+1.)*(X-1.)*(Y+1.))/8.
      DER(3,4)= ((X+1.)*(X-1.)*(Y+1.))/4.
      DER(3,5)=-((X+Y-2.*Z-1.)*(X+1.)*(Y+1.))/8.
      DER(3,6)= ((X+1.)*(Y+1.)*(Y-1.))/4.
      DER(3,7)= ((X-Y-2.*Z-1.)*(X+1.)*(Y-1.))/8.
      DER(3,8)=-((X+1.)*(X-1.)*(Y-1.))/4.
      DER(3,9)=-((X-1.)*(Y-1.)*Z)/2.
      DER(3,10)= ((X-1.)*(Y+1.)*Z)/2.
      DER(3,11)=-((X+1.)*(Y+1.)*Z)/2.
      DER(3,12)= ((X+1.)*(Y-1.)*Z)/2.
      DER(3,13)=-((X+Y-2.*Z+1.)*(X-1.)*(Y-1.))/8.
      DER(3,14)= ((X-1.)*(Y+1.)*(Y-1.))/4.
      DER(3,15)= ((X-Y-2.*Z+1.)*(X-1.)*(Y+1.))/8.
      DER(3,16)=-((X+1.)*(X-1.)*(Y+1.))/4.
      DER(3,17)= ((X+Y+2.*Z-1.)*(X+1.)*(Y+1.))/8.
      DER(3,18)=-((X+1.)*(Y+1.)*(Y-1.))/4.
      DER(3,19)=-((X-Y+2.*Z-1.)*(X+1.)*(Y-1.))/8.
      DER(3,20)= ((X+1.)*(X-1.)*(Y-1.))/4.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FM26NB (DER,IDER,FUN,SAMP,ISAMP,I)
C
C     This forms the shape functions and their
C     derivatives for 26-noded brick elements
C         Created using REDUCE by Dan Kidger,  March  1992
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      X = SAMP(I,1)
      Y = SAMP(I,2)
      Z = SAMP(I,3)
      FUN(1)=-((X*Y+X*Z+X+Y*Z+Y+Z+1.)*(X-1.)*(Y-1.)*(Z-1.))/8.
      FUN(2)= ((X+Z+1.)*(X-1.)*(Y+1.)*(Y-1.)*(Z-1.))/4.
      FUN(3)=-((X*Y-X*Z-X+Y*Z+Y-Z-1.)*(X-1.)*(Y+1.)*(Z-1.))/8.
      FUN(4)= ((X+1.)*(X-1.)*(Y-Z-1.)*(Y+1.)*(Z-1.))/4.
      FUN(5)=-((X*Y-X*Z-X-Y*Z-Y+Z+1.)*(X+1.)*(Y+1.)*(Z-1.))/8.
      FUN(6)= ((X-Z-1.)*(X+1.)*(Y+1.)*(Y-1.)*(Z-1.))/4.
      FUN(7)=-((X*Y+X*Z+X-Y*Z-Y-Z-1.)*(X+1.)*(Y-1.)*(Z-1.))/8.
      FUN(8)= ((X+1.)*(X-1.)*(Y+Z+1.)*(Y-1.)*(Z-1.))/4.
      FUN(9)= ((X+Y+1.)*(X-1.)*(Y-1.)*(Z+1.)*(Z-1.))/4.
      FUN(10)=-((X-Y+1.)*(X-1.)*(Y+1.)*(Z+1.)*(Z-1.))/4.
      FUN(11)=-((X+Y-1.)*(X+1.)*(Y+1.)*(Z+1.)*(Z-1.))/4.
      FUN(12)= ((X-Y-1.)*(X+1.)*(Y-1.)*(Z+1.)*(Z-1.))/4.
      FUN(13)= ((X*Y-X*Z+X-Y*Z+Y-Z+1.)*(X-1.)*(Y-1.)*(Z+1.))/8.
      FUN(14)=-((X-Z+1.)*(X-1.)*(Y+1.)*(Y-1.)*(Z+1.))/4.
      FUN(15)= ((X*Y+X*Z-X-Y*Z+Y+Z-1.)*(X-1.)*(Y+1.)*(Z+1.))/8.
      FUN(16)=-((X+1.)*(X-1.)*(Y+Z-1.)*(Y+1.)*(Z+1.))/4.
      FUN(17)= ((X*Y+X*Z-X+Y*Z-Y-Z+1.)*(X+1.)*(Y+1.)*(Z+1.))/8.
      FUN(18)=-((X+Z-1.)*(X+1.)*(Y+1.)*(Y-1.)*(Z+1.))/4.
      FUN(19)= ((X*Y-X*Z+X+Y*Z-Y+Z-1.)*(X+1.)*(Y-1.)*(Z+1.))/8.
      FUN(20)=-((X+1.)*(X-1.)*(Y-Z+1.)*(Y-1.)*(Z+1.))/4.
      FUN(21)=-((X+1.)*(X-1.)*(Y+1.)*(Y-1.)*(Z-1.))/2.
      FUN(22)= ((X+1.)*(X-1.)*(Y+1.)*(Y-1.)*(Z+1.))/2.
      FUN(23)=-((X+1.)*(X-1.)*(Y-1.)*(Z+1.)*(Z-1.))/2.
      FUN(24)= ((X+1.)*(X-1.)*(Y+1.)*(Z+1.)*(Z-1.))/2.
      FUN(25)=-((X-1.)*(Y+1.)*(Y-1.)*(Z+1.)*(Z-1.))/2.
      FUN(26)= ((X+1.)*(Y+1.)*(Y-1.)*(Z+1.)*(Z-1.))/2.
      DER(1,1)=-((2.*X*Y+2.*X*Z+2.*X+Y*Z)*(Y-1.)*(Z-1.))/8.
      DER(1,2)= ((2.*X+Z)*(Y+1.)*(Y-1.)*(Z-1.))/4.
      DER(1,3)=-((2.*X*Y-2.*X*Z-2.*X+Y*Z)*(Y+1.)*(Z-1.))/8.
      DER(1,4)= ((Y-Z-1.)*(Y+1.)*(Z-1.)*X)/2.
      DER(1,5)=-((2.*X*Y-2.*X*Z-2.*X-Y*Z)*(Y+1.)*(Z-1.))/8.
      DER(1,6)= ((2.*X-Z)*(Y+1.)*(Y-1.)*(Z-1.))/4.
      DER(1,7)=-((2.*X*Y+2.*X*Z+2.*X-Y*Z)*(Y-1.)*(Z-1.))/8.
      DER(1,8)= ((Y+Z+1.)*(Y-1.)*(Z-1.)*X)/2.
      DER(1,9)= ((2.*X+Y)*(Y-1.)*(Z+1.)*(Z-1.))/4.
      DER(1,10)=-((2.*X-Y)*(Y+1.)*(Z+1.)*(Z-1.))/4.
      DER(1,11)=-((2.*X+Y)*(Y+1.)*(Z+1.)*(Z-1.))/4.
      DER(1,12)= ((2.*X-Y)*(Y-1.)*(Z+1.)*(Z-1.))/4.
      DER(1,13)= ((2.*X*Y-2.*X*Z+2.*X-Y*Z)*(Y-1.)*(Z+1.))/8.
      DER(1,14)=-((2.*X-Z)*(Y+1.)*(Y-1.)*(Z+1.))/4.
      DER(1,15)= ((2.*X*Y+2.*X*Z-2.*X-Y*Z)*(Y+1.)*(Z+1.))/8.
      DER(1,16)=-((Y+Z-1.)*(Y+1.)*(Z+1.)*X)/2.
      DER(1,17)= ((2.*X*Y+2.*X*Z-2.*X+Y*Z)*(Y+1.)*(Z+1.))/8.
      DER(1,18)=-((2.*X+Z)*(Y+1.)*(Y-1.)*(Z+1.))/4.
      DER(1,19)= ((2.*X*Y-2.*X*Z+2.*X+Y*Z)*(Y-1.)*(Z+1.))/8.
      DER(1,20)=-((Y-Z+1.)*(Y-1.)*(Z+1.)*X)/2.
      DER(1,21)=-(Y+1.)*(Y-1.)*(Z-1.)*X
      DER(1,22)= (Y+1.)*(Y-1.)*(Z+1.)*X
      DER(1,23)=-(Y-1.)*(Z+1.)*(Z-1.)*X
      DER(1,24)= (Y+1.)*(Z+1.)*(Z-1.)*X
      DER(1,25)=-((Y+1.)*(Y-1.)*(Z+1.)*(Z-1.))/2.
      DER(1,26)= ((Y+1.)*(Y-1.)*(Z+1.)*(Z-1.))/2.
      DER(2,1)=-((2.*X*Y+X*Z+2.*Y*Z+2.*Y)*(X-1.)*(Z-1.))/8.
      DER(2,2)= ((X+Z+1.)*(X-1.)*(Z-1.)*Y)/2.
      DER(2,3)=-((2.*X*Y-X*Z+2.*Y*Z+2.*Y)*(X-1.)*(Z-1.))/8.
      DER(2,4)= ((X+1.)*(X-1.)*(2.*Y-Z)*(Z-1.))/4.
      DER(2,5)=-((2.*X*Y-X*Z-2.*Y*Z-2.*Y)*(X+1.)*(Z-1.))/8.
      DER(2,6)= ((X-Z-1.)*(X+1.)*(Z-1.)*Y)/2.
      DER(2,7)=-((2.*X*Y+X*Z-2.*Y*Z-2.*Y)*(X+1.)*(Z-1.))/8.
      DER(2,8)= ((X+1.)*(X-1.)*(2.*Y+Z)*(Z-1.))/4.
      DER(2,9)= ((X+2.*Y)*(X-1.)*(Z+1.)*(Z-1.))/4.
      DER(2,10)=-((X-2.*Y)*(X-1.)*(Z+1.)*(Z-1.))/4.
      DER(2,11)=-((X+2.*Y)*(X+1.)*(Z+1.)*(Z-1.))/4.
      DER(2,12)= ((X-2.*Y)*(X+1.)*(Z+1.)*(Z-1.))/4.
      DER(2,13)= ((2.*X*Y-X*Z-2.*Y*Z+2.*Y)*(X-1.)*(Z+1.))/8.
      DER(2,14)=-((X-Z+1.)*(X-1.)*(Z+1.)*Y)/2.
      DER(2,15)= ((2.*X*Y+X*Z-2.*Y*Z+2.*Y)*(X-1.)*(Z+1.))/8.
      DER(2,16)=-((X+1.)*(X-1.)*(2.*Y+Z)*(Z+1.))/4.
      DER(2,17)= ((2.*X*Y+X*Z+2.*Y*Z-2.*Y)*(X+1.)*(Z+1.))/8.
      DER(2,18)=-((X+Z-1.)*(X+1.)*(Z+1.)*Y)/2.
      DER(2,19)= ((2.*X*Y-X*Z+2.*Y*Z-2.*Y)*(X+1.)*(Z+1.))/8.
      DER(2,20)=-((X+1.)*(X-1.)*(2.*Y-Z)*(Z+1.))/4.
      DER(2,21)=-(X+1.)*(X-1.)*(Z-1.)*Y
      DER(2,22)= (X+1.)*(X-1.)*(Z+1.)*Y
      DER(2,23)=-((X+1.)*(X-1.)*(Z+1.)*(Z-1.))/2.
      DER(2,24)= ((X+1.)*(X-1.)*(Z+1.)*(Z-1.))/2.
      DER(2,25)=-(X-1.)*(Z+1.)*(Z-1.)*Y
      DER(2,26)= (X+1.)*(Z+1.)*(Z-1.)*Y
      DER(3,1)=-((X*Y+2.*X*Z+2.*Y*Z+2.*Z)*(X-1.)*(Y-1.))/8.
      DER(3,2)= ((X+2.*Z)*(X-1.)*(Y+1.)*(Y-1.))/4.
      DER(3,3)=-((X*Y-2.*X*Z+2.*Y*Z-2.*Z)*(X-1.)*(Y+1.))/8.
      DER(3,4)= ((X+1.)*(X-1.)*(Y-2.*Z)*(Y+1.))/4.
      DER(3,5)=-((X*Y-2.*X*Z-2.*Y*Z+2.*Z)*(X+1.)*(Y+1.))/8.
      DER(3,6)= ((X-2.*Z)*(X+1.)*(Y+1.)*(Y-1.))/4.
      DER(3,7)=-((X*Y+2.*X*Z-2.*Y*Z-2.*Z)*(X+1.)*(Y-1.))/8.
      DER(3,8)= ((X+1.)*(X-1.)*(Y+2.*Z)*(Y-1.))/4.
      DER(3,9)= ((X+Y+1.)*(X-1.)*(Y-1.)*Z)/2.
      DER(3,10)=-((X-Y+1.)*(X-1.)*(Y+1.)*Z)/2.
      DER(3,11)=-((X+Y-1.)*(X+1.)*(Y+1.)*Z)/2.
      DER(3,12)= ((X-Y-1.)*(X+1.)*(Y-1.)*Z)/2.
      DER(3,13)= ((X*Y-2.*X*Z-2.*Y*Z-2.*Z)*(X-1.)*(Y-1.))/8.
      DER(3,14)=-((X-2.*Z)*(X-1.)*(Y+1.)*(Y-1.))/4.
      DER(3,15)= ((X*Y+2.*X*Z-2.*Y*Z+2.*Z)*(X-1.)*(Y+1.))/8.
      DER(3,16)=-((X+1.)*(X-1.)*(Y+2.*Z)*(Y+1.))/4.
      DER(3,17)= ((X*Y+2.*X*Z+2.*Y*Z-2.*Z)*(X+1.)*(Y+1.))/8.
      DER(3,18)=-((X+2.*Z)*(X+1.)*(Y+1.)*(Y-1.))/4.
      DER(3,19)= ((X*Y-2.*X*Z+2.*Y*Z+2.*Z)*(X+1.)*(Y-1.))/8.
      DER(3,20)=-((X+1.)*(X-1.)*(Y-2.*Z)*(Y-1.))/4.
      DER(3,21)=-((X+1.)*(X-1.)*(Y+1.)*(Y-1.))/2.
      DER(3,22)= ((X+1.)*(X-1.)*(Y+1.)*(Y-1.))/2.
      DER(3,23)=-(X+1.)*(X-1.)*(Y-1.)*Z
      DER(3,24)= (X+1.)*(X-1.)*(Y+1.)*Z
      DER(3,25)=-(X-1.)*(Y+1.)*(Y-1.)*Z
      DER(3,26)= (X+1.)*(Y+1.)*(Y-1.)*Z
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FM27NB (DER,IDER,FUN,SAMP,ISAMP,I)
C
C     This forms the shape functions and their
C     derivatives for 27-noded brick elements
C         Created using REDUCE by Dan Kidger,  March  1992
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      X = SAMP(I,1)
      Y = SAMP(I,2)
      Z = SAMP(I,3)
      FUN(1)= ((X-1.)*(Y-1.)*(Z-1.)*X*Y*Z)/8.
      FUN(2)=-((X-1.)*(Y+1.)*(Y-1.)*(Z-1.)*X*Z)/4.
      FUN(3)= ((X-1.)*(Y+1.)*(Z-1.)*X*Y*Z)/8.
      FUN(4)=-((X+1.)*(X-1.)*(Y+1.)*(Z-1.)*Y*Z)/4.
      FUN(5)= ((X+1.)*(Y+1.)*(Z-1.)*X*Y*Z)/8.
      FUN(6)=-((X+1.)*(Y+1.)*(Y-1.)*(Z-1.)*X*Z)/4.
      FUN(7)= ((X+1.)*(Y-1.)*(Z-1.)*X*Y*Z)/8.
      FUN(8)=-((X+1.)*(X-1.)*(Y-1.)*(Z-1.)*Y*Z)/4.
      FUN(9)=-((X-1.)*(Y-1.)*(Z+1.)*(Z-1.)*X*Y)/4.
      FUN(10)=-((X-1.)*(Y+1.)*(Z+1.)*(Z-1.)*X*Y)/4.
      FUN(11)=-((X+1.)*(Y+1.)*(Z+1.)*(Z-1.)*X*Y)/4.
      FUN(12)=-((X+1.)*(Y-1.)*(Z+1.)*(Z-1.)*X*Y)/4.
      FUN(13)= ((X-1.)*(Y-1.)*(Z+1.)*X*Y*Z)/8.
      FUN(14)=-((X-1.)*(Y+1.)*(Y-1.)*(Z+1.)*X*Z)/4.
      FUN(15)= ((X-1.)*(Y+1.)*(Z+1.)*X*Y*Z)/8.
      FUN(16)=-((X+1.)*(X-1.)*(Y+1.)*(Z+1.)*Y*Z)/4.
      FUN(17)= ((X+1.)*(Y+1.)*(Z+1.)*X*Y*Z)/8.
      FUN(18)=-((X+1.)*(Y+1.)*(Y-1.)*(Z+1.)*X*Z)/4.
      FUN(19)= ((X+1.)*(Y-1.)*(Z+1.)*X*Y*Z)/8.
      FUN(20)=-((X+1.)*(X-1.)*(Y-1.)*(Z+1.)*Y*Z)/4.
      FUN(21)= ((X+1.)*(X-1.)*(Y+1.)*(Y-1.)*(Z-1.)*Z)/2.
      FUN(22)= ((X+1.)*(X-1.)*(Y+1.)*(Y-1.)*(Z+1.)*Z)/2.
      FUN(23)= ((X+1.)*(X-1.)*(Y-1.)*(Z+1.)*(Z-1.)*Y)/2.
      FUN(24)= ((X+1.)*(X-1.)*(Y+1.)*(Z+1.)*(Z-1.)*Y)/2.
      FUN(25)= ((X-1.)*(Y+1.)*(Y-1.)*(Z+1.)*(Z-1.)*X)/2.
      FUN(26)= ((X+1.)*(Y+1.)*(Y-1.)*(Z+1.)*(Z-1.)*X)/2.
      FUN(27)=-(X+1.)*(X-1.)*(Y+1.)*(Y-1.)*(Z+1.)*(Z-1.)
      DER(1,1)= (((X-1.)+X)*(Y-1.)*(Z-1.)*Y*Z)/8.
      DER(1,2)=-(((X-1.)+X)*(Y+1.)*(Y-1.)*(Z-1.)*Z)/4.
      DER(1,3)= (((X-1.)+X)*(Y+1.)*(Z-1.)*Y*Z)/8.
      DER(1,4)=-((Y+1.)*(Z-1.)*X*Y*Z)/2.
      DER(1,5)= (((X+1.)+X)*(Y+1.)*(Z-1.)*Y*Z)/8.
      DER(1,6)=-(((X+1.)+X)*(Y+1.)*(Y-1.)*(Z-1.)*Z)/4.
      DER(1,7)= (((X+1.)+X)*(Y-1.)*(Z-1.)*Y*Z)/8.
      DER(1,8)=-((Y-1.)*(Z-1.)*X*Y*Z)/2.
      DER(1,9)=-(((X-1.)+X)*(Y-1.)*(Z+1.)*(Z-1.)*Y)/4.
      DER(1,10)=-(((X-1.)+X)*(Y+1.)*(Z+1.)*(Z-1.)*Y)/4.
      DER(1,11)=-(((X+1.)+X)*(Y+1.)*(Z+1.)*(Z-1.)*Y)/4.
      DER(1,12)=-(((X+1.)+X)*(Y-1.)*(Z+1.)*(Z-1.)*Y)/4.
      DER(1,13)= (((X-1.)+X)*(Y-1.)*(Z+1.)*Y*Z)/8.
      DER(1,14)=-(((X-1.)+X)*(Y+1.)*(Y-1.)*(Z+1.)*Z)/4.
      DER(1,15)= (((X-1.)+X)*(Y+1.)*(Z+1.)*Y*Z)/8.
      DER(1,16)=-((Y+1.)*(Z+1.)*X*Y*Z)/2.
      DER(1,17)= (((X+1.)+X)*(Y+1.)*(Z+1.)*Y*Z)/8.
      DER(1,18)=-(((X+1.)+X)*(Y+1.)*(Y-1.)*(Z+1.)*Z)/4.
      DER(1,19)= (((X+1.)+X)*(Y-1.)*(Z+1.)*Y*Z)/8.
      DER(1,20)=-((Y-1.)*(Z+1.)*X*Y*Z)/2.
      DER(1,21)= (Y+1.)*(Y-1.)*(Z-1.)*X*Z
      DER(1,22)= (Y+1.)*(Y-1.)*(Z+1.)*X*Z
      DER(1,23)= (Y-1.)*(Z+1.)*(Z-1.)*X*Y
      DER(1,24)= (Y+1.)*(Z+1.)*(Z-1.)*X*Y
      DER(1,25)= (((X-1.)+X)*(Y+1.)*(Y-1.)*(Z+1.)*(Z-1.))/2.
      DER(1,26)= (((X+1.)+X)*(Y+1.)*(Y-1.)*(Z+1.)*(Z-1.))/2.
      DER(1,27)=-2.*(Y+1.)*(Y-1.)*(Z+1.)*(Z-1.)*X
      DER(2,1)= (((Y-1.)+Y)*(X-1.)*(Z-1.)*X*Z)/8.
      DER(2,2)=-((X-1.)*(Z-1.)*X*Y*Z)/2.
      DER(2,3)= (((Y+1.)+Y)*(X-1.)*(Z-1.)*X*Z)/8.
      DER(2,4)=-(((Y+1.)+Y)*(X+1.)*(X-1.)*(Z-1.)*Z)/4.
      DER(2,5)= (((Y+1.)+Y)*(X+1.)*(Z-1.)*X*Z)/8.
      DER(2,6)=-((X+1.)*(Z-1.)*X*Y*Z)/2.
      DER(2,7)= (((Y-1.)+Y)*(X+1.)*(Z-1.)*X*Z)/8.
      DER(2,8)=-(((Y-1.)+Y)*(X+1.)*(X-1.)*(Z-1.)*Z)/4.
      DER(2,9)=-(((Y-1.)+Y)*(X-1.)*(Z+1.)*(Z-1.)*X)/4.
      DER(2,10)=-(((Y+1.)+Y)*(X-1.)*(Z+1.)*(Z-1.)*X)/4.
      DER(2,11)=-(((Y+1.)+Y)*(X+1.)*(Z+1.)*(Z-1.)*X)/4.
      DER(2,12)=-(((Y-1.)+Y)*(X+1.)*(Z+1.)*(Z-1.)*X)/4.
      DER(2,13)= (((Y-1.)+Y)*(X-1.)*(Z+1.)*X*Z)/8.
      DER(2,14)=-((X-1.)*(Z+1.)*X*Y*Z)/2.
      DER(2,15)= (((Y+1.)+Y)*(X-1.)*(Z+1.)*X*Z)/8.
      DER(2,16)=-(((Y+1.)+Y)*(X+1.)*(X-1.)*(Z+1.)*Z)/4.
      DER(2,17)= (((Y+1.)+Y)*(X+1.)*(Z+1.)*X*Z)/8.
      DER(2,18)=-((X+1.)*(Z+1.)*X*Y*Z)/2.
      DER(2,19)= (((Y-1.)+Y)*(X+1.)*(Z+1.)*X*Z)/8.
      DER(2,20)=-(((Y-1.)+Y)*(X+1.)*(X-1.)*(Z+1.)*Z)/4.
      DER(2,21)= (X+1.)*(X-1.)*(Z-1.)*Y*Z
      DER(2,22)= (X+1.)*(X-1.)*(Z+1.)*Y*Z
      DER(2,23)= (((Y-1.)+Y)*(X+1.)*(X-1.)*(Z+1.)*(Z-1.))/2.
      DER(2,24)= (((Y+1.)+Y)*(X+1.)*(X-1.)*(Z+1.)*(Z-1.))/2.
      DER(2,25)= (X-1.)*(Z+1.)*(Z-1.)*X*Y
      DER(2,26)= (X+1.)*(Z+1.)*(Z-1.)*X*Y
      DER(2,27)=-2.*(X+1.)*(X-1.)*(Z+1.)*(Z-1.)*Y
      DER(3,1)= (((Z-1.)+Z)*(X-1.)*(Y-1.)*X*Y)/8.
      DER(3,2)=-(((Z-1.)+Z)*(X-1.)*(Y+1.)*(Y-1.)*X)/4.
      DER(3,3)= (((Z-1.)+Z)*(X-1.)*(Y+1.)*X*Y)/8.
      DER(3,4)=-(((Z-1.)+Z)*(X+1.)*(X-1.)*(Y+1.)*Y)/4.
      DER(3,5)= (((Z-1.)+Z)*(X+1.)*(Y+1.)*X*Y)/8.
      DER(3,6)=-(((Z-1.)+Z)*(X+1.)*(Y+1.)*(Y-1.)*X)/4.
      DER(3,7)= (((Z-1.)+Z)*(X+1.)*(Y-1.)*X*Y)/8.
      DER(3,8)=-(((Z-1.)+Z)*(X+1.)*(X-1.)*(Y-1.)*Y)/4.
      DER(3,9)=-((X-1.)*(Y-1.)*X*Y*Z)/2.
      DER(3,10)=-((X-1.)*(Y+1.)*X*Y*Z)/2.
      DER(3,11)=-((X+1.)*(Y+1.)*X*Y*Z)/2.
      DER(3,12)=-((X+1.)*(Y-1.)*X*Y*Z)/2.
      DER(3,13)= (((Z+1.)+Z)*(X-1.)*(Y-1.)*X*Y)/8.
      DER(3,14)=-(((Z+1.)+Z)*(X-1.)*(Y+1.)*(Y-1.)*X)/4.
      DER(3,15)= (((Z+1.)+Z)*(X-1.)*(Y+1.)*X*Y)/8.
      DER(3,16)=-(((Z+1.)+Z)*(X+1.)*(X-1.)*(Y+1.)*Y)/4.
      DER(3,17)= (((Z+1.)+Z)*(X+1.)*(Y+1.)*X*Y)/8.
      DER(3,18)=-(((Z+1.)+Z)*(X+1.)*(Y+1.)*(Y-1.)*X)/4.
      DER(3,19)= (((Z+1.)+Z)*(X+1.)*(Y-1.)*X*Y)/8.
      DER(3,20)=-(((Z+1.)+Z)*(X+1.)*(X-1.)*(Y-1.)*Y)/4.
      DER(3,21)= (((Z-1.)+Z)*(X+1.)*(X-1.)*(Y+1.)*(Y-1.))/2.
      DER(3,22)= (((Z+1.)+Z)*(X+1.)*(X-1.)*(Y+1.)*(Y-1.))/2.
      DER(3,23)= (X+1.)*(X-1.)*(Y-1.)*Y*Z
      DER(3,24)= (X+1.)*(X-1.)*(Y+1.)*Y*Z
      DER(3,25)= (X-1.)*(Y+1.)*(Y-1.)*X*Z
      DER(3,26)= (X+1.)*(Y+1.)*(Y-1.)*X*Z
      DER(3,27)=-2.*(X+1.)*(X-1.)*(Y+1.)*(Y-1.)*Z
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------


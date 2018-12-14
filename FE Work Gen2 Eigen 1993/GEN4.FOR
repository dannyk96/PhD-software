      OPTIONS  ( DREAL, INTL, FULLCHECK )
               P  R  O  G  R  A  M       G  E  N  2 

C-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C     D A N ' S    2 / 3 D  F E - P R O G R A M                       !
C-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C    This is a re-hash of my code for FOS of 3d slopes                !
C    ( Originally BOOK6.4 ) .. this version is heavily modified and   !
C    tidied-up .. MAY 1991                                            !
C                                                                     !
C-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C   once          (     PROGRAM 6.4 THREE-DIMENSIONAL ELASTO-PLASTIC  !
C    upon         (     ANALYSIS USING 20-NODE BRICK ELEMENTS         !
C     a time ...  (     MOHR-COULOMB FAILURE CRITERION                !
C-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
C --- latest changes ...
C       17-10-91  PARMETER's shuffled, FOS looping added, MUNG extended
C***********************************************************************
 
C     Max. element size parameters        Max. problem size params
C     ----------------------------        -------------------------
C        MN   = nodes per elem.           INF = nodes in mesh
C        MDOF = Dof. per mode             MEL = elements in mesh
C        MDF  = freedoms per elem         IKV = global stiffness matrix
C               
C        NS   = strains per GP            ILDS  = loaded nodes per mesh
C        IGRULE = GP's per elem           IQINC = load steps
C                                         IPRPTAB = material types
C***********************************************************************

C ----> change the next line to change problem size
c      PARAMETER (IKV=2 600 000, INF=3500, MEL=600
       PARAMETER (IKV=84 500, INF=6000, MEL=915
c      PARAMETER (IKV=9 000, INF=60, MEL=33
C      PARAMETER (IKV=150 100, INF=3000, MEL=600
     +                   ,ILDS=50,IQINC=50,IPRTAB=20)

C ---->    change the next line to change largest element size
C      PARAMETER ( MN = 8, MDF = 2 )
      PARAMETER ( MN = 20, MDF = 3 )

      PARAMETER ( IGRULE =27 )
      PARAMETER ( MDOF = MN*MDF, NS   = MDF*2 )
      PARAMETER ( ILOADS = INF*MDF  * 7/10 )

C ----> change the next line to =1 if running on the VP
      PARAMETER ( IVP=0 )
C ---- 'dependant' array sizes ----
      PARAMETER ( IEVPT =MEL*NS*IGRULE, ISMP = IGRULE,
     +     IDEE=NS, IBEE=NS+1, IJAC=MDF,IDER=MDF,IDERIV=MDF,
     +     ICOORD=MN, IKM=MDOF, IGC=INF, INUMS=MEL, IKDIAG=ILOADS)

C***********************************************************************
C ------- optional COMMON block for the largest arrays ---------------
c     COMMON/COM/KV, GC,NUMS,EVPT,EVSTRS
C    +        ,KM, LOADS,OLDIS,GRAVLO,TOTDIS,KDIAG
C ----------------- problem size depandant arrays --------------------
      REAL   KV  (IKV)
      REAL GC(IGC,MDF),
     +     EVSTRS(IEVPT),   EVPT  (IEVPT),
     +     LOADS (ILOADS),  OLDIS (ILOADS), TOTD(ILOADS),
     +     BDYLDS(ILOADS),  GRAVLO(ILOADS)
      INTEGER NF (INF,MDF), KDIAG (IKDIAG),  MATBM(10)
     +       ,NUMS(INUMS,MN+5),    NO (ILDS)            , QN  (IQINC)
      REAL    PRPTAB(IPRTAB,10), VAL(ILDS),  REACT(ILDS), QINC(IQINC)
C ---------------------- standard FE arrays --------------------------
      INTEGER     G(MDOF), NUM(MN), U0,U1,U2,U3,U4
      REAL DEE(IDEE,IDEE), FUN(MN),       COORD (ICOORD,MDF),
     +     BEE(IBEE,MDOF), DER(IDER,MN),  DERIV (IDERIV,MN),
     +     KM (IKM, IKM ), SMP(ISMP,MDF), WTS(ISMP),  JAC(IJAC,IJAC),
     +     ELD(MDOF),      EPS(NS),       SIGMA(NS),  STRESS(NS),
     +     GAMA(MDF),      XYZP(3), XYZBM(6), BMST(IGRULE)
C --------------------------------------------------------------------
      CHARACTER FILE*60, TIME*8
C .. push next to its' sub ?
      PARAMETER(RTD=180./3.14159265)
      DATA U0,U1,U2,U3,U4/10,11,12,13,14/, ITYPE/5/
C =====================================================================
C =========================== DATA INPUT ==============================
C =====================================================================
      IF(IVP.NE.1) THEN
        CALL FILES2(FILE,LEN)
        OPEN(U0  ,FILE='DAT/'//FILE//'.DAT')
        OPEN(U1  ,FILE=  'TEMP.DAT')
        OPEN(U2  ,FILE='RES/'//FILE//'.RES')
        OPEN(U3  ,FILE='PL/' //FILE//'.PL ')
        OPEN(U4  ,FILE='RES/'//FILE//'.BM1')
        OPEN(U4+1,FILE='RES/'//FILE//'.BM2')
        OPEN(U4+2,FILE='RES/'//FILE//'.BM3')
      ENDIF
      CALL CSTRIP(U0,U1)
      READ(U1,*) NODOF, NGP, IOPCYL
C -------------- read the basic mesh in ----------------------
      CALL READGM(U1,NODOF,NN,GC,IGC,NEL,NUMS,INUMS)
C ---------------- move nodes if required --------------------
      CALL NODMOV(U1,NN,GC,IGC,NNMOV,NODOF,IOPCYL)
C ----------- get material propertty groups ------------------
      CALL GETELT(U1,NUMS,INUMS,NEL,NODOF)

C ------ read in property types as E,v,c,phi,psi,cons,bw --------------
      READ(U1,*) NPROPS,NTYPES, ((PRPTAB(I,J),J=1,NTYPES),I=1,NPROPS)

C ------------ delete unwanted elements and close up the 'gaps' -------
      CALL ELDEL(NODOF,NN,GC,IGC,NEL,NUMS,INUMS,NF,INF)

C ------------------------ set up a few values --------------------
      IF(NODOF.EQ.2) IH = 4
      IF(NODOF.EQ.3) IH = 6
      NEVPT  = NEL * IH * NGP

C ---------------------- boundary conditions --------------------------
C ----- set ALL nodes to be freedoms !
      DO 95,I=1,NN
      DO 95,J=1,NODOF
   95 NF(I,J)=1
      CALL GETNF  (U1,NF,INF,GC,IGC,NODOF,NN)
      CALL SORTNF (NF,INF,NN,NODOF,N)

C ----------------- form the KDIAG pointer vector -------------------
      CALL FMKDAG(NEL,NUMS,INUMS,NF,INF,G,NODOF,N,IW,KDIAG,IKDIAG)

C ------------------------- read loadings ------------------------------
C --- first read in the X,Y,Z gravity factors (eg. 0.,-1.,0.)
C -- data as max # incs. , start load, X,Y axes lengths

      READ(U1,*) (GAMA(I),I=1,NODOF)
      READ(U1,*) INCS, XO, XS,YS, ITS
      CALL READLD (U1,NN,GC,IGC,NF,INF,NODOF,NO,VAL,NL)
      READ(U1,*) INCS,(QN(I),QINC(I),I=1,INCS)
      INCST=0
      DO 211,I=1,INCS
  211 INCST = INCST+QN(I)
C ------------------------- read printout options----------------------
C ... print frequency (always print-out first,last & end of each block)
      READ(U1,*)IPFREQ
C ... bending moments
      READ(U1,*) NBM,(MATBM(I),I=1,NBM)
C ========================== end of data input ========================
      CLOSE (U1)

C ----------------- transform the model coords ------------------------
      IF (IOPCYL.NE.0) THEN
        DO 201,I=1,NN
        DO 202,J=1,NODOF
          XYZP(1) = GC(I,1)
          XYZP(2) = GC(I,2)
          IF(NODOF.EQ.2) XYZP( 3 ) = 0.
  202     IF(NODOF.EQ.3) XYZP( 3 ) = GC(I,MDF)

          CALL MUNG (XYZP,IOPCYL)

        DO 201,J=1,NODOF
  201   GC(I,J) = XYZP(J)
      ENDIF
C ------------ output the geometry and steering ----------------------

      WRITE(U3,*)NODOF
      WRITE(U3,*)NN
      DO 101,I=1,NN
  101 WRITE(U3,'(I6,3F12.5)')I,(GC(I,J),J=1,NODOF)
      WRITE(U3,*)NEL
      DO 102,I=1,NEL
  102 WRITE(U3,'(25I6)')I,(NUMS(I,J),J=1,NUMS(I,1)+2)
C     STOP '** geometry written OK'
C --------------------- check array sizes ------------------------

      WRITE(*,1001)'nodes'    ,'NN  = ',NN      ,' INF   =',INF
      WRITE(*,1001)'elements' ,'NEL = ',NEL     ,' MEL   =',MEL
      WRITE(*,1001)'loadings' ,'NL  = ',NL      ,' ILDS  =',ILDS
      WRITE(*,1001)'freedoms' ,'N   = ',N       ,' ILOADS=',ILOADS
      WRITE(*,1001)'KV length','IR  = ',KDIAG(N),' IKV   =',IKV
      WRITE(*,1001)'G.points' ,'NEVPT=',NEVPT   ,' IEVPT =',IEVPT

 1001 FORMAT(A10,A10,I8,A10,I8)
      IF (NN.GT.INF .OR. KDIAG(N).GT.IKV .OR. NEL.GT.MEL .OR.
     +    NEVPT.GT.IEVPT .OR. N.GT.ILOADS) THEN
          PRINT*,'******* ARRAYS TOO SMALL ********'
          STOP
      ENDIF
C ----------------------- null some global matrices -------------------
      CALL NULVEC (KV,KDIAG(N))
      CALL NULVEC (TOTD,N)
      CALL NULVEC (OLDIS,N)
      CALL NULVEC (GRAVLO,N)
C ------------ form some matricess ----------------------
      CALL MYGAUS (SMP,ISMP,WTS,NODOF,NGP)

C =====================================================================
C ============ ELEMENT STIFFNESS INTEGRATION AND ASSEMBLY =============
C =====================================================================

C --------------- loop through the elements ---------------------------
      DO 3,IEL=1,NEL
      print*,'element',iel,' /',nel,char(27)//'[A'
      CALL GEOM1(IEL,NUMS,INUMS,GC,IGC,NOD,COORD,ICOORD,NUM,NODOF)
      CALL NUM2G (NF,INF,NOD,NUM,NODOF,G)

C ------- get material properties, hence form DEE ---------------------
      IDOF = NOD * NODOF
      IMAT = NUMS (IEL,NOD+2)
      E    = PRPTAB (IMAT   ,1)
      V    = PRPTAB (IMAT   ,2)
      RHO  = PRPTAB (IMAT   ,6)
      BW   = PRPTAB (IMAT   ,7)
      CALL FRMD (DEE,IDEE,E,V,BW,NODOF,1)
      CALL NULL (BEE,IBEE,IH,IDOF)

C ------- extract the local co-ordinates COORD and the steering G -----
      CALL NULL   (KM,IKM,IDOF,IDOF)
      CALL NULVEC (ELD,IDOF)

C ----------------------- loop the G.P.'s -----------------------------
      DO 4 IGP=1,NGP
        CALL GSF  (NODOF,NOD,ITYPE, DER,IDER,FUN,SMP,ISMP,IGP)
      CALL MATMUL (DER,IDER,COORD,ICOORD,JAC,IJAC,NODOF,NOD,NODOF)
      CALL INVERT (JAC,IJAC,DET,NODOF)
      CALL MATMUL (JAC,IJAC,DER,IDER,DERIV,IDERIV,NODOF,NODOF,NOD)
      SUM=1.
      IF(NODOF.EQ.2)THEN
        IF(IOPCYL.EQ.0)CALL FORMB  (BEE,IBEE,DERIV,IDERIV,NOD)
        IF(IOPCYL.EQ.1)CALL FMBRAD
     +        (BEE,IBEE,DERIV,IDERIV,FUN,COORD,ICOORD,SUM,NOD)
      ENDIF
      IF(NODOF.EQ.3)CALL FORMB3 (BEE,IBEE,DERIV,IDERIV,NOD)
      WEIGHT= SUM * DET * WTS(IGP)
      CALL FMBTDB (BEE,IBEE,DEE,IDEE,KM,IKM, WEIGHT, IH,IDOF)
      DO 4 L=1,IDOF
    4 ELD(L) = ELD(L) + WEIGHT * FUN((L-1)/NODOF+1)

      CALL FSPARV(KV,KM,IKM,G,KDIAG,IDOF)
c      write(*,'(a,i3,a,i3,a,a)')
c     +    'element',iel,' /',nel,' completed  '
      DO 3 J=1,NODOF
      DO 3 K=J,IDOF,NODOF
    3 IF(G(K).NE.0) GRAVLO(G(K)) = GRAVLO(G(K)) + ELD(K) * GAMA(J)*RHO
C =============== end of matrix assembly ==============================
      IR = KDIAG(N)
      DO 1963 I=1,IR
      IF(ABS(KV(I)).LE.0.0001) J = J +1
 1963 CONTINUE
      WRITE(*,'(I6,A,I6,A,F8.3,A)')  J ,' ZEROS +',IR-J,' NON-ZEROS ='
     + ,100.*REAL(IR-J)/REAL(IR),' %'

C ----------- set the 'big-springs' at displacemnt loaded nodes -------
      HUGE=0.
      DO 92,I=1,KDIAG(N)
   92 HUGE = MAX(HUGE,(KV(I)))
      HUGE = HUGE * 10000.
      DO 93,I=1,NL
   93 IF(NO(I).LT.0) KV(KDIAG(ABS( NO(I)) ))  = HUGE
      PRINT*, 'NOTE: big spring stiffness = HUGE = ', HUGE,'   '

C======================== REDUCE EQUATIONS ============================
       CALL SPARIN(KV,N,KDIAG)

C**********************************************************************

C======================================================================
C===================  FACTOR OF SAFETY LOOP  ==========================
C======================================================================
C >> need to toggle wether I am load-stepping or FOS stepping
C >> really this is only relavent in the AUTOFOS bit !
c      DO 777,IFOS=1,NFOS
       FOS = 1.

C ------------------ initial stresses and strains ---------------
C >> on the second pass these should be available on file.
      CALL NULVEC (SIGMA,IH)
      CALL NULVEC (EVPT,  NEVPT)
      IN = -IH
      DO 91,IEL=1,NEL
      CONS = PRPTAB(NUMS(IEL,NUMS(IEL,1)+2),6)
      SIGMA(1) = CONS
      SIGMA(2) = CONS
      SIGMA(3) = CONS
        IF(NODOF.EQ.2) SIGMA(4) = SIGMA(3)
        IF(NODOF.EQ.2) SIGMA(3) = 0.
      DO 91,I=1,NGP
      IN = IN+IH
      DO 91,J=1,IH
   91 EVSTRS( IN +J ) = SIGMA(J)

C =====================================================================
C ====================== LOAD INCREMENT LOOP ==========================
C =====================================================================

      IY1  = 1
      IY2  = 0
      PTOT = 0

      DO 700 IY=1,INCST
      PRINT*,'LOAD STEP',iy,' /',INCST,CHAR(27)//'[A'
      IY2 = IY2+1
      IF (IY2.GT.QN(IY1)) THEN
          IY1 = IY1+1
          IY2 = 1
      ENDIF

      STEP = QINC(IY1)
      PTOT = PTOT + STEP

      CALL NULVEC (BDYLDS,N)
      CALL NULVEC (EVPT,NEVPT)

C ====================== iteration loop ===============================
      DO 701,ITERS=1,ITS
      CALL NULVEC(LOADS,N)

C ------ get load/disp increment data from NO and VAL --------
      DO 42,I=1,NL
        WEIGHT = 1.
        IF(NO(I).LT.0) WEIGHT = HUGE
   42   LOADS (ABS(NO(I))) = STEP * VAL(I) * WEIGHT

C        --- ADD THE GRAVITY LOADS ON LOAD STEP #1
      IF(IY.EQ.1) CALL VECADD (LOADS,GRAVLO,LOADS,N)
      CALL VECADD (LOADS,BDYLDS,LOADS,N)
      CALL SPABAC (KV,LOADS,N,KDIAG)
C ---------------------- check convergence ----------------------------
      CALL CHECON(LOADS,OLDIS,N,.001,ICON)

C ... on the first pass force 'no-convergence' flag
      IF (ITERS.EQ.1) ICON=0
      IF (ICON.EQ.1.OR.ITERS.EQ.ITS) CALL NULVEC(BDYLDS,N)
      IF (ICON.EQ.1.OR.ITERS.EQ.ITS) CALL NULVEC(REACT ,NL)
      IOPPR=0
      IF(ICON.EQ.1.AND. (IY.EQ.1.OR.IY.EQ.INCST.OR. IY2.EQ.QN(IY1)
     +                   . OR.MOD(IY,IPFREQ).EQ.0))    IOPPR=1

C --------------- loop through the elements ---------------------------
      IN = -IH
      DO 777,IEL=1,NEL
      CALL GEOM1(IEL,NUMS,INUMS,GC,IGC,NOD,COORD,ICOORD,NUM,NODOF)
      CALL NUM2G (NF,INF,NOD,NUM,NODOF,G)

C ------- get material properties, hence form DEE ---------------------
      IDOF = NOD * NODOF
      IMAT = NUMS (IEL,NOD+2)
      E    = PRPTAB (IMAT   ,1)
      V    = PRPTAB (IMAT   ,2)
      C    = PRPTAB (IMAT   ,3)
      PHI  = PRPTAB (IMAT   ,4)
      PSI  = PRPTAB (IMAT   ,5)
      BW   = PRPTAB (IMAT   ,7)
      CALL FRMD (DEE,IDEE,E,V,BW,NODOF,1)
      CALL NULL (BEE,IBEE,IH,IDOF)
      CALL NULL (KM,IKM,IDOF,IDOF)
      IWBM=0
      IF(IOPPR.EQ.1.OR.IY.EQ.1)THEN
        DO 280,I=1,NBM
  280   IF(IMAT.EQ.MATBM(I)) IWBM=I
C       IF(IWBM.EQ.0)THEN
C         DO 281,I=1,NODOF+2
C 281     XYZBM(I)=0.
C       ENDIF
      ENDIF
C ---- get the local freedom forces -------
      DO 9, M=1,IDOF
        IF (G(M).EQ.0) ELD(M) = 0.
    9   IF (G(M).NE.0) ELD(M) = LOADS(G(M))
C ----------------------- loop the G.P.'s -----------------------------
C >>>> note that for 2d IH=4 for the following lines
      DO 778 IGP=1,NGP
      IN=IN+IH
        CALL GSF  (NODOF,NOD,ITYPE, DER,IDER,FUN,SMP,ISMP,IGP)
      CALL MATMUL (DER,IDER,COORD,ICOORD,JAC,IJAC,NODOF,NOD,NODOF)
      CALL INVERT (JAC,IJAC,DET,NODOF)
      CALL MATMUL (JAC,IJAC,DER,IDER,DERIV,IDERIV,NODOF,NODOF,NOD)
      SUM=1.
      IF(NODOF.EQ.2)THEN
        IF(IOPCYL.EQ.0)CALL FORMB  (BEE,IBEE,DERIV,IDERIV,NOD)
        IF(IOPCYL.EQ.1)CALL FMBRAD
     +        (BEE,IBEE,DERIV,IDERIV,FUN,COORD,ICOORD,SUM,NOD)
      ENDIF
      IF (NODOF.EQ.3) CALL FORMB3 (BEE,IBEE,DERIV,IDERIV,NOD)
      WEIGHT= SUM * DET * WTS(IGP)
C ----> ! what is the next for ??
C      PRINT*,IEL,IGP,SUM,DET,WTS(IGP)
      CALL MVMULT (BEE,IBEE,ELD,IH,IDOF,EPS)
      DO 120 L=1,IH
  120 EPS(L) = EPS(L) - EVPT(IN+L)
        CALL MVMULT(DEE,IDEE,EPS,IH,IH,SIGMA)
      DO 121,L=1,IH
  121 STRESS(L) = SIGMA(L) + EVSTRS(IN+L)
      CALL VECCOP(STRESS,SIGMA,IH)

      IF(NODOF.EQ.2) CALL INVAR  (STRESS,SIGM,DSBAR,THETA)
      IF(NODOF.EQ.3) CALL INVAR3 (STRESS,SIGM,DSBAR,THETA)

C ============= CHECK WHETHER YIELD IS VIOLATED ==================

C ---------- factor PHI and C if necessary --------------
      C    =           C           /FOS
      PHI  = ATAN (TAN(PHI/RTD)    /FOS) *RTD
      SNPH = SIN  (PHI/RTD)
      DT =     4.*(1.+V)*(1.-2.*V)      /
     +     (E*(1.-2.*V+SIN(PHI/RTD)**2 ))
      CALL MOCOUF(PHI,C, SIGM,DSBAR,THETA ,F)

C ---- if failed (F>0.) then put the 'plastic flow' stresses into SIGMA

C >> note: in general I CAN do plastic flow if ICON=1 if I want!
      IF (ICON.NE.1 .AND. ITERS.NE.ITS) THEN
      IF(F.LT.0.)GOTO 127
      CALL MCFLOW (DSBAR,THETA, PSI,DT,F,DEE,IDEE,
     +                       EVPT,IN,IH,NODOF, STRESS)
      ENDIF
        DO 122,L=1,IH
          T = STRESS(L) * WEIGHT
          DO 122 M=1,IDOF
C         PRINT*,IEL,IGP,IY,ITERS,M,L,   STRESS(L),WEIGHT
  122     IF(G(M).NE.0) BDYLDS(G(M)) = BDYLDS(G(M)) + BEE(L,M) * T

C --------- update stresses only if converged -------------------------
C ---- also sum the reaction forces into REACT  :-)    -----

  127 CONTINUE
      IF(ICON.EQ.1.OR.ITERS.EQ.ITS)THEN
        DO 123,L=1,IH
  123   EVSTRS(IN+L) = SIGMA(L)
C --- loop reaction nodes, find where(if) this Gp contributes
C --- if so, loop thru IH and sum into REACT, BEE*STRESS
        DO 831,I=1,NL
          DO 832,M=1,IDOF
            IF (G(M).EQ.ABS(NO(I))) THEN
              DO 833,L=1,IH
                REACT(I) = REACT(I) + BEE(L,M) * EVSTRS(IN+L) * WEIGHT
  833         CONTINUE
              GOTO 831
            ENDIF
  832     CONTINUE
  831   CONTINUE
      ENDIF
C ---------------- calculate bending moments if required --------------

      IF(IWBM.NE.0) BMST(IGP)=EVSTRS(IN+2)

      IF (IWBM.NE.0.AND.IY.EQ.1.AND.ITERS.EQ.1) THEN
        CALL GCOORD (FUN,COORD,ICOORD,NOD,NODOF,XYZP)
        WRITE(U4+IWBM-1,'(2I4,3F12.6)') IEL,IGP,(XYZP(J),J=1,NODOF)
      ENDIF

C ---------------------------------------------------------------------
  778 CONTINUE
      IF (IWBM.NE.0.AND.ITERS.NE.1)
     + WRITE(U4+IWBM-1,'(2I5,10E12.4,/,10X,10E12.4,/,10X,10E12.4)')
     +                       IY,IEL, (BMST(I),I=1,NGP)

  777 CONTINUE
C --------- end of GP and elem loops .. can do some monitoring --------
      IF (ICON.NE.0) GOTO 702
  701 CONTINUE
C --------- end of ITS loop - only here if no-convergence -----
      PRINT*, '*** NO-CONVERGENCE *** in load-step',IY
C --------------------------------------------------------------
  702 CONTINUE
C ------ we have convergence (or hit ITS limit) so print some values --
      CALL VECADD(TOTD,LOADS,TOTD,N)
      WRITE(U2,'(  E12.4,I6)') PTOT, ITERS
      WRITE(U2,'(12E12.4   )')(REACT(       L  ),L=1,NL)
      WRITE(U2,'(12E12.4   )')(TOTD (ABS(NO(L))),L=1,NL)
C--- print the total force and the mean displacemnt of the loaded nodes
      S1=0.
      S2=0.
      IF(IY.EQ.1)S2OLD=0.
      DO 841,L=1,NL
      S1 = S1 + TOTD (ABS(NO(L)))
  841 S2 = S2 + REACT(       L  )
      WRITE ( *,'(I3,3F18.9,I5,A9)') IY,S2,S1/NL,S2-S2OLD,ITERS,TIME()
C     WRITE (U4,'(I3,3F18.9,I5,A9)') IY,S2,S1/NL,S2-S2OLD,ITERS,TIME()
C      WRITE (U2,'( 2F20.10   )')
C     +    -S1/NL /.14,   -S2/67./(3.14159*.07**2) *4.
      S2OLD = S2
C--- output the nodal displacements (cf. nodal co-ords)
      IF (IOPPR.EQ.1) THEN
        WRITE(U3,'(A,I3,A,E12.5)')'#     Load-step No.',IY,' LOAD=',S2
        WRITE(U3,*)NN
        DO 103,I=1,NN
        DO 104,J=1,NODOF
          IF(NF(I,J).EQ.0)          XYZP(J) = 0.
  104     IF(NF(I,J).NE.0)          XYZP(J) = TOTD(NF(I,J))
  103     WRITE(U3,'(I5,3E12.4)')I,(XYZP(J),J=1,NODOF)
C------- old LOADS output (only 2d with no pile is actually usefull !)
C      WRITE(19,'(/,I12)') IY
C      WRITE(19,'(1X,6E12.4)')(LOADS(I),I=1,N)
      ENDIF
  700 CONTINUE
C --- end of the load steeping loop.. print final state ? ------
c  778 CONTINUE
C --- end of the FOS  steeping loop.. print final state ? ------
      STOP
      END
C =====================================================================
C ================= end of program ====================================
C =====================================================================

C      WRITE(8,'(A,F7.3)')'FOS=',FOS
C     IF(NOD.EQ.20) CALL LIN(ITERS,LOADS,NF,INF,ICR,2*NZE)
C     IF(NOD.EQ.14) CALL LIN(ITERS,LOADS,NF,INF,ICR,  NZE)
C     WRITE(8,*)
C     CALL PRNTV (8,LOADS,N)
C     CALL MAXPRN(LOADS,NF,INF,ICR,ITERS)
C     CALL MAXDS(LOADS,ILOADS,NF,INF,NN,NXE,NZE)
C      CALL TABLE2(FOS,ITERS,.00, LOADS(1),1,
C    +     LOADS(NF(ICR+1,1)),     LOADS(NF(ICR+1,2)),ICR)
C ----------------------------------------------------------------------
      SUBROUTINE MAXDS2 (LOADS,ILOADS,NF,INF,NN)
C
C     TO CALC. MAX DISPLACEMENTS - this will be the rewrite in 1/2/3D
C     also add the a final column of the maximum nodal magnitude
C
      REAL LOADS(0:ILOADS-1)
      INTEGER I,J,NN,NF(INF,*)
      DO 4,I=ILOADS-1,1,-1
    4 LOADS(I)=LOADS(I-1)
      LOADS(0)=0.
      XM=0.
      YM=0.
      ZM=0.
      DO 10,J=1,NN
      IF(ABS(LOADS(NF(J,1))).GT.ABS(XM))THEN
        IX=J
        XM=LOADS(NF(J,1))
      ENDIF
      IF(ABS(LOADS(NF(J,2))).GT.ABS(YM))THEN
        IY=J
        YM=LOADS(NF(J,2))
      ENDIF
      IF(ABS(LOADS(NF(J,3))).GT.ABS(ZM))THEN
        IZ=J
        ZM=LOADS(NF(J,3))
      ENDIF
   10 CONTINUE
        WRITE(*,'(3(I5,F18.11))')IX,XM,IY,YM,IZ,ZM
      DO 5,I=1,ILOADS-1
    5 LOADS(I-1)=LOADS(I)
      RETURN
      END
C***********************************************************************
      SUBROUTINE READLD (IO,NN,GC,IGC,NF,INF,NODOF,NO,VAL,NL)
C
C     this reads in the load at a node based on its coordinate
C
C  -- if IDIR  -ve then a presc. disp.
      REAL GC(IGC,*),POINT(5),VAL(*)
      INTEGER NO(*), NF(INF,*), NUMLST(5)
      READ(IO,*) NL
      IC = 0
      DO 1, I=1,NL
      READ(IO,*)  (POINT(J),J=1,NODOF),  IDIR,VALUE
      CALL FNDNOD (NN,GC,IGC,POINT,NODOF,NUMLST,K)
      IF(K.EQ. 0) THEN
        PRINT*,'** WARNING : no nodes found for load #',I
        GOTO 1
      ENDIF
      IF(K.GT. 1) PRINT*,'** WARNING : >1 node  found for load #',I
      VAL(I) = VALUE
      NO (I) = SIGN(NF(NUMLST(1),ABS(IDIR)), IDIR  )
      IF (NO(I).EQ.0)PRINT*,'*** WARNING : attempt to load a fixity',I
    1 CONTINUE
      RETURN
      END
C***********************************************************************
        SUBROUTINE MUNG (XYZP,IOP)
C
C       this MUNGS a co-ord from one co-ord system to another
C
C       eg. IOP = 1 cartesian to cylindricl
C           IOP = 2 cartesian to spherical
C           IOP = 3 'spiral' a cartesian in the z
C           IOP = 4  ... insert your own rule ??
C
        REAL XYZP(*)

        XO = XYZP(1)
        YO = XYZP(2)
        ZO = XYZP(3)
        RTD = 180. / 3.14159265
      IF (IOP.EQ.1) THEN

C ---- apply cylindical transformation
        YO = -YO
        XYZP(1) =  YO * COS(ZO/RTD)
        XYZP(3) =  YO * SIN(ZO/RTD)
        XYZP(2) =  XO
      ELSEIF (IOP.EQ.2) THEN

C ----- apply a 'twist' in the y-direction
        ANG =  YO*180/16. /RTD
        XYZP(1) =  XO * COS(ANG) - ZO * SIN(ANG)
        XYZP(3) =  XO * SIN(ANG) + ZO * COS(ANG)
        XYZP(2) =  YO
C       XYZP(2) = -YO
      ELSEIF (IOP.EQ.3) THEN
C ----- transform to a sphere radius = 4.-ish
      RADBX = MAX(ABS(XO),ABS(YO),ABS(ZO))
      RADSP = SQRT (XO**2 + YO**2 + ZO**2 + 1.E-10)
      XYZP(1) = XO * RADBX / RADSP * 2.
      XYZP(2) = YO * RADBX / RADSP * 2.
      XYZP(3) = ZO * RADBX / RADSP * 2.
      ENDIF
      RETURN
      END
C***********************************************************************
      SUBROUTINE NODMOV(IO,NN,GC,IGC,NNMOV,NODOF,IOPCYL)
C
C     this moves sets of nodes to specified positions
C     (note  230964.e10 = wildcard)
C
      REAL GC(IGC,*), COLD(3), CNEW(3)
      INTEGER NUMLST(40)

      WILD = 230964.E10
      READ(IO,*) NNMOV
      DO 1,K=1 , NNMOV
        READ(IO,*) (COLD(I),I=1,NODOF) ,(CNEW(I),I=1,NODOF)
        CALL FNDNOD (NN,GC,IGC,COLD,NODOF,NUMLST,NNOD)
        IF(NNOD.EQ. 0) THEN
          PRINT*,'** WARNING : no nodes found for node move #',NNOD
          GOTO 1
        ENDIF
        DO 2,I=1,NNOD
        DO 2,J=1,NODOF
          IF (ABS(COLD(J)-WILD)/WILD.GT.1.E-4)
     +                         GC(NUMLST(I),J) = CNEW(J)
   2  CONTINUE
   1  CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE FILES2 (FILE,L)
C
C     this returns the root name of the data files FILE
C
      CHARACTER FILE*(*)
      OPEN(10,FILE='FILE')
      READ(10,'(A)')FILE
      L = INDEX (FILE,' ') -1 
      PRINT*,'file="',FILE(1:L),'"'
      RETURN
      END
C***********************************************************************
      SUBROUTINE CSTRIP (U0,U1)
C
C     this copies a data file from unit U0 to unit U1
C     stripping out all characters after a '#' comment mark
C
      INTEGER U0,U1
      CHARACTER LINE*(80)
      REWIND U0
      REWIND U1
      DO 1,I=1,99999
      READ(U0,'(A)', END=99) LINE
      IP = INDEX(LINE,'#')-1
      IF (IP.EQ.0) GOTO 1
      IF (IP.EQ.-1) IP=LEN(LINE)
        DO 2,J=IP,1,-1
          IF (LINE(J:J).NE.' ')GOTO 3
    2   CONTINUE
    3   IP=J
      IF(IP.GT.0)WRITE(U1,'(A)')LINE(1:IP)
    1 CONTINUE
   99 REWIND U0
      REWIND U1
      END
C***********************************************************************
      SUBROUTINE ELDEL(NODOF,NN,GC,IGC,NEL,NUMS,INUMS,NF,INF)
C
C     this will delete any elements of material type=0 in a mesh
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*),NF(INF,*)

      NN1  = NN
      NEL1 = NEL
C ----- use column 1 of NF as a working space to index the nodes -------
      DO 1,I=1,NN
    1 NF(I,1)=0

      IC=0
      DO 2,IEL=1,NEL
        NOD  = NUMS(IEL,1)
        IMAT = NUMS(IEL,NOD+2)
C --- if a 'good' element .. activate all its nodes and store it at IC
        IF(IMAT.NE.0)THEN
          IC = IC + 1
          DO 3,I=2,NOD+1
            IF (IMAT.NE.0) NF(NUMS(IEL,I),1)=1
    3     CONTINUE
          DO 4,I=1,NOD+5
            NUMS(IC,I)=NUMS(IEL,I)
    4     CONTINUE
        ENDIF
    2 CONTINUE
      NEL = IC

C --- now delete the unreferenced nodes ------
C ----- label each valid node and update the co-ordinate list..
      IC = 0
      DO 5,I=1,NN
        IF(NF(I,1).NE.0)THEN
          IC = IC + 1
          NF(I,1) = IC
            DO 6,J=1,NODOF
    6       GC(IC,J) = GC (I,J)
        ENDIF
    5 CONTINUE
      NN = IC
C --- update the nodal pointers in NUMS
      DO 7,IEL=1,NEL
        NOD = NUMS(IEL,1)
          DO 7,J=2,NOD+1
          NUMS(IEL,J) = NF(NUMS(IEL,J),1)
    7 CONTINUE
      WRITE(*,'(A,I4,A,I4,A)')
     +     'NOTE :',NEL1-NEL,' elements (',NN1-NN,' ) nodes deleted.'
      RETURN
      END
C***********************************************************************
      CHARACTER*8 FUNCTION TIME()
C
C     returns the 'time' as a character string
C
      character time@*8
C      TIME='hh:mm'
      TIME=time@()
      RETURN
      END
C***********************************************************************
C***********************************************************************



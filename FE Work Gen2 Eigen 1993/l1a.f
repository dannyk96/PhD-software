      OPTIONS (fullcheck,undef,dreal,intl)
C***********************************************************************
C    this is my library L1A .. full of 'mainstream' routines          **
C    see L1B for 'secondary' routines                                 **
C    THESE ARE THE LIBRARY OF FINITE ELEMENT SUBROUTINES              **
C    FOR 'SLOPE' TYPE MODELS .. ie. beams, plates piles, etc. etc.    **
C    WRITTEN BY DANIEL KIDGER BETWEEN 1987 and 1989 FOR MY            **
C     PH.D. WORK ... YOU MAY COPY THEM FREELY BUT LIKE GIVE           **
C     ME THE CREDIT Y'KNOW.                                           **
C                                                                     **
C***********************************************************************

C************************************************************************
C*****
C***** routines using generalised geometry and steering GCOORD/NUMS
C*****
C************************************************************************

      SUBROUTINE READGM(IO,NODOF,NN,GCOORD,IGCRD,NEL,NUMS,INUMS)
C
C     this routine reads the element geometry puts it into 'long' form
C
      REAL GCOORD(IGCRD,*), COORD(20,3)
      REAL    TOP(90), BOT(90), DEPTH(90), BRDTH(90)
      INTEGER NUM(20), NUMS(INUMS,*)

      READ(IO,*) NOD, NXE,NXS, NYE,NYS 
      NZE=1
      IF(NODOF.EQ.3) READ(IO,*) NZE
      NEL   = NZE*(NXE*NYE-NYS*(NXE-NXS))
      IF(NEL.GT.INUMS) PRINT*,'** WARNING: NUMS too small; NEL=',NEL
      READ(IO,*) (TOP(I),I=1,NXS+1)
      READ(IO,*) (BOT(I),I=1,NXE+1)
      READ(IO,*) (DEPTH(I),I=1,NYE+1)
      IF(NODOF.EQ.3)READ(IO,*) (BRDTH(I),I=1,NZE+1)
C --------------- form the global co-ord array GCOORD ------------

      NN  = 0 
      IEL = 0
      DO 6,IP=1,NXE
      DO 6,IQ=1+NYS*MIN(MAX(IP-NXS,0),1),NYE
      DO 6,IS=1,NZE
      IEL = IEL + 1
      NUMS(IEL,     1) = NOD
      NUMS(IEL, NOD+2) = 1
      NUMS(IEL, NOD+3) = IP
      NUMS(IEL, NOD+4) = IQ
      NUMS(IEL, NOD+5) = IS

      IF(NODOF.EQ.2)CALL SLOGE2(IP,IQ,NXS,NYS,NYE,TOP,BOT,DEPTH
     +       ,COORD,20,NUM)
      IF(NODOF.EQ.3)
     +CALL SLOGET(IP,IQ,IS,NXS,NYS,NYE,NZE,TOP,BOT,DEPTH,BRDTH
     +       ,COORD,20,NUM,NOD)
      DO 7,I=1,NOD
        NUMS(IEL,I+1) = NUM(I)
        NN =      MAX (NN,NUM(I))
        IF (NN.GT.IGCRD) PRINT*,'** WARNING: GCOORD too small NN>',NN
        DO 7,J=1,NODOF
    7   GCOORD(NUM(I),J)=COORD(I,J)
    6 CONTINUE
C     WRITE(*,'(I5,3F12.5)')(I,(GCOORD(I,J),J=1,3),I=1,NN)
      RETURN
      END
C ------------------------------------------------------------------
      SUBROUTINE GETELT (IO,NUMS,INUMS,NEL,NODOF) 
C                
C     this reads in the element material types in x,y,z from/to
C     format - storing the material property type in NUMS after
C     the number of nodes and node numbers 
C     -- data in the format /mat type/,/xfrom,xto/, etc

      INTEGER NUMS(INUMS,*), XF(5),XT(5)
      READ(IO,*)NMAT
      DO 1,L=1,NMAT
      READ(IO,*) MATYPE, (XF(I),XT(I),I=1,NODOF)
      IC = 0
      DO 2,I=1,NEL
        NOD = NUMS (I,1)
        DO 3,IDIR=1,NODOF
          IP = NUMS(I,NOD+2+IDIR)
          IF (XF(IDIR).GT. XT(IDIR))
     +          PRINT*,'** WARNING  XF>XT (GETELT)'
          IF (IP.LT.XF(IDIR).OR.IP.GT.XT(IDIR)) GOTO 2 
    3   CONTINUE
        NUMS(I,NOD+2) = MATYPE 
      IC=IC+1
    2 CONTINUE
      IF(IC.EQ.0) PRINT*,'*** WARNING: no MATERIAL #', NMAT,' found'  
    1 CONTINUE
      END
C ----------------------------------------------------------------
      SUBROUTINE GETNF (IO,NF,INF,GCOORD,IGCRD,NODOF,NN)
C
C     this reads in a set of boundary conditions based on  
C     'X-from, X-to, Y-from, Y-to, Z-from, Z-to' co-ords
C      where GCOORD contains all the nodal co-ordinates
C   >> note set NF to all non-zero first and call SORTNF afterwards
C
      REAL GCOORD(IGCRD,*), XMIN(5), XMAX(5)
      INTEGER  NF(INF,   *), BC(5)
      TOL=  1.E-20
      READ(IO,*)NR
      DO 1,I=1,NR
        READ(IO,*) (BC(J),J=1,NODOF),(XMIN(J),XMAX(J),J=1,NODOF)
        IC = 0
        DO 2,J=1,NN
          DO 4,K=1,NODOF
            IF (GCOORD(J,K).LT.XMIN(K)-TOL) GOTO 2
            IF (GCOORD(J,K).GT.XMAX(K)+TOL) GOTO 2
    4     CONTINUE
          DO 3,K=1,NODOF
    3     NF(J,K) = BC(K) * NF(J,K)
          IC=IC+1
    2   CONTINUE
        PRINT*,'FOR BC #',I,IC,'points were found'
        IF(IC.EQ.0)PRINT*,'*** WARNING no nodes found for BC #',I
    1 CONTINUE
      RETURN
      END
C ------------------------------------------------------------------
      SUBROUTINE FMKDAG(NEL,NUMS,INUMS,NF,INF,G,NODOF,N,IW,KDIAG,IKDIAG)
C 
C     this calculates the max. (and mean) bandwidth 
C     optionally KDIAG may be formed ( only formed if N<IKDIAG )
C
      INTEGER KDIAG(*), G(*),NUMS(INUMS,*), NF(INF,*), NUM(20)
      IWMAX = 0
      IWSUM = 0
      IF(N.LE.IKDIAG)THEN
        DO 1 I=1,N
    1   KDIAG(I)=0
      ELSE
        PRINT*,'warning: KDIAG not formed as N>IKIAG'
      ENDIF
      DO 3 IEL=1,NEL
        NOD  = NUMS (IEL,1)
        IDOF = NOD * NODOF
        DO 2,I=1,NOD
    2   NUM(I)=NUMS(IEL,1+I)
        CALL NUM2G (NF,INF,NOD,NUM,NODOF,G)
        IF(N.LT.IKDIAG) CALL FKDIAG(KDIAG,G,IDOF)
        MING = N
        MAXG = 0
        DO 4,I=1,IDOF
                      MAXG = MAX (MAXG,G(I))
    4   IF(G(I).NE.0) MING = MIN (MING,G(I))
        IW    = MAXG - MING
        IWMAX = MAX(IWMAX,(IW))
    3   IWSUM = IWSUM + IW

      IF(N.LT.IKDIAG)THEN
        KDIAG(1)=1
        DO 5 I=2,N
    5   KDIAG(I)=KDIAG(I)+KDIAG(I-1)
      ENDIF
      IW = IWMAX
      PRINT*,'total IW=',IW,' average IW=',IWSUM/NEL
      END
C************************************************************************
      SUBROUTINE GEOM1(IEL,NUMS,INUMS,GCOORD,IGCRD,NOD,COORD,ICOORD
     +                                                     ,NUM,NODOF)
C
C     general geometry routine - returns COORD from the global GCOORD 
C                                also G and NOD
C                                   
      REAL GCOORD(IGCRD,*), COORD(ICOORD,*)
      INTEGER NUM(*), NUMS(INUMS,*)
      NOD  = NUMS (IEL,1)
      IC=0
      DO 1, I =1,NOD
      DO 2, IFREE = 1,NODOF
      IC = IC + 1
      INODE = NUMS (IEL,1+I)
    2 COORD(I,IFREE) = GCOORD(INODE,IFREE)
    1 NUM(I) = NUMS(IEL,1+I)
      RETURN
      END
C************************************************************************
      SUBROUTINE FNDNOD(NN,GCOORD,IGCRD,POINT,NODOF,NUMLST,IC)
C
C     given a co-ordinate, this returns its node number
C
      REAL GCOORD(IGCRD,*),POINT(*)
      INTEGER NUMLST(*)
      TOL  = 1.E-12
      WILD = 230964.E10
      IC=0
      DO 1,I=1,NN
      DO 2,J=1,NODOF
      IF (ABS(POINT(J)-  WILD     ).LT.TOL) GOTO 2
      IF (ABS(POINT(J)-GCOORD(I,J)).GT.TOL) GOTO 1
    2 CONTINUE
      IC=IC+1
      NUMLST(IC)=I
    1 CONTINUE
      RETURN
      END
C************************************************************************
C      ____________________________________________________________      
C     /                                                          / \     
C    /                                                          /   \   
C   /  u s e f u l l   e x t e n s i o n s   t o   F E 5 L I B /     \ 
C  /                                                          /_______\
C /__________________________________________________________/         

      SUBROUTINE INVERT (JAC,IJAC,DET,NODOF)
C
C     this calls TWOBY2,TREEX3 as needed, putting JAC1 back into JAC
C
      REAL JAC(IJAC,*), JAC1(3,3)
      IF (NODOF.EQ.1) THEN
        DET = JAC(1,1)
        JAC(1,1) = 1./DET
      ELSE IF (NODOF.EQ.2) THEN
        CALL TWOBY2 (JAC,IJAC,JAC1,3,DET)
        CALL MATCOP (JAC1,3,JAC,IJAC,NODOF,NODOF)
      ELSE IF (NODOF.EQ.3) THEN
        CALL TREEX3 (JAC,IJAC,JAC1,3,DET)
        CALL MATCOP (JAC1,3,JAC,IJAC,NODOF,NODOF)
      ELSE
        CALL MATINV(JAC,IJAC,NODOF)
        PRINT*,'** WARNING, DET not calculated (INVERT)' 
      ENDIF
        IF(DET.LT.0.)PRINT*,'** WARNING negative Jacobian (INVERT)' 
      END
C ------------------------------------------------------------------
      SUBROUTINE SORTNF(NF,INF,NN,NODOF,N)
C
C     this resequences all the non-zero terms in NF to an ascending 
C     sequence and returns the total number of freedoms in N
C
      INTEGER NF(INF,*)
      N = 0
      DO 1,I=1,NN
      DO 1,J=1,NODOF
      IF (NF(I,J).NE.0) THEN
        N = N + 1
        NF(I,J) = N
      ENDIF
    1 CONTINUE
      RETURN
      END
C ------------------------------------------------------------------
      SUBROUTINE NUM2G(NF,INF,NOD,NUM,NODOF,G)
C
C     this short routine copies the NF data into the steering G
C
      INTEGER NF(INF,*),NUM(*),G(*)
      IC=0
      DO 1,I=1,NOD
      DO 1,J=1,NODOF
      IC=IC+1
    1 G(IC) = NF(NUM(I),J)
      RETURN
      END
C************************************************************************

      SUBROUTINE MCFLOW (DSBAR,THETA,PSI,DT,F, DEE,IDEE,EVPT,IN, 
     +                        IH,NODOF, SIGMA)
C
C     ... GIVEN(PSI,SIGMA,DT  ) .. GIVES 'DEVP' = PLASTIC STRAIN BIT
C .... could really pass it all the material parameters 
C  .. and updates evpt ...

      REAL M1(6,6), M2(6,6), M3(6,6), FLOW(6,6), EVPT(*), EVP(6)
      REAL T1(4,4), T2(4,4), T3(4,4), ERATE(6) , SIGMA(*),DEE(IDEE,*)
      IF (NODOF.EQ.2) THEN
        CALL FORMM  (SIGMA,T1,T2,T3)
        CALL MATCOP (T1,4,M1,6,4,4)
        CALL MATCOP (T2,4,M2,6,4,4)
        CALL MATCOP (T3,4,M3,6,4,4)
      ELSEIF (NODOF.EQ.3) THEN
        CALL FORMM3 (SIGMA,M1,M2,M3)
      ELSE
        STOP 'MCFLOW - NOT 2D OR 3D'
      ENDIF
      CALL MOCOUQ(PSI,DSBAR,THETA,DQ1,DQ2,DQ3)
      DO 1 L=1,IH   
      DO 1 M=1,IH
    1 FLOW(L,M)=F*(M1(L,M)*DQ1+M2(L,M)*DQ2+M3(L,M)*DQ3)
      CALL MVMULT (FLOW,6,SIGMA,IH,IH,ERATE)
      DO 2 L=1,IH
      EVP(L) = ERATE(L) * DT
    2 EVPT(IN+L) = EVPT(IN+L)+EVP(L)
      CALL MVMULT (DEE,IDEE,EVP,IH,IH,SIGMA)
      RETURN
      END
C********************************************************************
      SUBROUTINE FMBTDB(BEE,IBEE,DEE,IDEE,KM,IKM,FACT,IH,IDOF)
C
C     ** this is a mega-simplification of the long BT * DEE * BEE
C        and MSMULT, MATADD proceedures. Simply pass it the BEE and
C        DEE matrices and this willl do the rest.
C     It has condensed looping, it avoids unnecessary work when
C     zero terms are found in DEE ( hence not much slow-down when doing
C     SRI). It also exploits the symmetry of KM thus halving the effort.
C
C     * Note it needs BEE to be one row longer (as a workspace)
C
C                                         Dan KIdger    May 1991
C
      REAL  BEE(IBEE,*), DEE(IDEE,*), KM(IKM,*)
      LOGICAL LSYM
      DATA LSYM/.TRUE./
      IWORK = IH+1
      DO 6,I=1,IH
      DO 6,J=1,IH
    6 LSYM = LSYM.AND. (ABS(DEE(I,J)-DEE(J,I)).LT.1.E-12)
      DO 1,K=1,IH

C --- first form one column of BT * DEE in the workspace
        DO 2,J=1,IDOF
    2   BEE(IWORK,J) = 0.0
        DO 3,I=1,IH
          F = DEE(I,K) * FACT
          IF (ABS(F).GT. 1.E-12) THEN
            DO 4,J=1,IDOF
    4       BEE(IWORK,J) =  BEE(IWORK,J) + F * BEE(I,J)
          ENDIF
    3   CONTINUE

C ---  now product this column with BEE and sum into KM
        DO 1,I=1,IDOF
        IF(LSYM)THEN
        DO 5,J=I,IDOF
          KM(I,J) = KM(I,J) + BEE(IWORK,I) * BEE(K,J)
    5     KM(J,I) = KM(I,J)
        ELSE
        DO 7,J=1,IDOF
    7     KM(I,J) = KM(I,J) + BEE(IWORK,I) * BEE(K,J)
        ENDIF
    1 CONTINUE
      RETURN
      END

C********************************************************************
      SUBROUTINE AUTO2 (X,Y,XO,YO,XS,YS,STEPL)
C
C       ---- this is the SEQUEL to AUTOFOS -----
C      X,Y are the old positions
C      STEPL is the non-dimensional step length  (eg 1. = 5% diagonal)
C      XO,YO are the old X,Y positions (auto-updated)
C      X,Y are the expected X/Y 'graph' ranges
C      --- typically one is used, the other compared
C
       DX   =  X - XO  /XS
       DY   = (Y - YO) /YS
       VLEN = SQRT (DX*DX+DY*DY) * 0.05
       XN = X + STEPL * DX / VLEN * XS
       YN = Y + STEPL * DY / VLEN * YS
       XO = X
       YO = Y
      RETURN
      END
C************************************************************************
      SUBROUTINE FRMD (DEE,IDEE,E,V,BW,NODOF,IOP)
C
C    THIS FORMS MOST STRESS/STRAIN MATRICES
C
C    SET IOP = (1)NORM, (2)K, (3)G, 4(LAME), 5(G),6(PLANE-STRESS)
C    BW = modulus of the water ( =0.0 if none)
C    
      REAL DEE(IDEE,*),FACT(10)
      INTEGER M(3,6)
      DATA M/5,4,2, 3,3,1, 7,8,2, 4,4,1, 6,1,2,  9,10,2/
      CALL NULL(DEE,IDEE,IDEE,IDEE)

C      fact= (1)0.,(2)G,(3)K,(4)LAME,(5)L2,(6)2G,(7)4G/3,(8)=2G/3

      FACT(1) = 0.0
      FACT(2) = E/2./(1.+V)
      FACT(3) = E/3./(1.-2.*V)              + BW
      FACT(4) =     V *E/(1.+V)/(1.-2.*V)   + BW
      FACT(5) = (1.-V)*E/(1.+V)/(1.-2.*V)   + BW
      FACT(6) =   2.*FACT(2)                  
      FACT(7) =   4.*FACT(2)/3.
      FACT(8) = - 2.*FACT(2)/3.
      FACT(9) =     E/(1.-V*V)              
      FACT(10)=   V*E/(1.-V*V)

      DO 2,I=1,NODOF
        DEE(I,I)=FACT(M(1,IOP))
        DO 2,J=1,NODOF
    2   IF(I.NE.J)DEE(I,J)=FACT(M(2,IOP))
      DO 1,I=NODOF+1, NODOF*(NODOF+1)/2
    1 DEE(I,I)=FACT(M(3,IOP))
C ... if 2d and DEE big enough add the 4th row
      IF(NODOF.EQ.2.AND.IDEE.GE.4)THEN
        DO 3,I=1,4
        DEE(4,I) = DEE(1,I)
    3   DEE(I,4) = DEE(4,I)
      ENDIF
      RETURN
      END
C*********************************************************************
C*****
C*****  s h a p e   f u n c t i o n   r o u t i n e s
C*****
C*********************************************************************

      SUBROUTINE GSF(NODOF,NOD,ITYPE,DER,IDER,FUN,SMP,ISMP,I)
C
C     a macro header routine to call up the shape functions
C     and derivatives of ANY known element.
C     ITYPE = optional element type
C                                       Dan Kidger   April '91
C
C .. samp order reversed for 4nq,8nq,8nb,11nb etc. 
C    so that the 'first' coord is ALWAYS local-X !         8-4-92
C
      REAL SMP(ISMP,*), DER(IDER,*), FUN(*), SAMP(3,2)

C ... quads etc. need to mung SMP into old SAMP format ...

      DO 1,J=1,NODOF
    1 SAMP(J,1) = SMP(I,J)

C ... need to add in 1d elements such as the beams etc.
      FUN(1) = 230964.e10
      IF (NODOF.EQ.1)THEN
        IF(NOD.EQ. 2) CALL FM02N1 (DER,IDER,FUN,SMP, ISMP,I)
        IF(NOD.EQ. 3) CALL FM03N1 (DER,IDER,FUN,SMP, ISMP,I)
      ELSEIF (NODOF.EQ.2)THEN
        IF(NOD.EQ. 3) CALL FMTRI3 (DER,IDER,FUN,SMP, ISMP,I)
        IF(NOD.EQ. 6) CALL FMTRI6 (DER,IDER,FUN,SMP, ISMP,I)
        IF(NOD.EQ.10) CALL FMTR10 (DER,IDER,FUN,SMP, ISMP,I)
        IF(NOD.EQ.15) CALL FMTR15 (DER,IDER,FUN,SMP, ISMP,I)
        IF(NOD.EQ. 4) CALL FORMLN (DER,IDER,FUN,SAMP, 3,2,1)
        IF(NOD.EQ. 5) CALL FM05N2 (DER,IDER,FUN,SAMP, 3,1,2)
        IF(NOD.EQ. 8) CALL FMQUAD (DER,IDER,FUN,SAMP, 3,2,1)
        IF(NOD.EQ. 9) CALL FMLAG9 (DER,IDER,FUN,SAMP, 3,2,1)
        IF(NOD.EQ.12) CALL FM12N2 (DER,IDER,FUN,SAMP, 3,1,2) 
        IF(NOD.EQ.17) CALL FM17N2 (DER,IDER,FUN,SAMP, 3,1,2) 
 
      ELSEIF (NODOF.EQ.3)THEN
        IF(NOD.EQ. 4) CALL FMTET4 (DER,IDER,FUN,SMP, ISMP,I)
        IF(NOD.EQ. 8) CALL FMLIN3 (DER,IDER,FUN,SAMP, 3,2,1,3)
        IF(NOD.EQ.11) CALL FMQM6  (DER,IDER,FUN,SAMP, 3,2,1,3)
        IF(NOD.EQ.14) CALL FM14   (DER,IDER,FUN,SAMP, 3,1,2,3,ITYPE)
        IF(NOD.EQ.20) CALL FMQUA3 (DER,IDER,FUN,SAMP, 3,1,2,3)
        IF(NOD.EQ.26) CALL FM26N3 (DER,IDER,FUN,SAMP, 3,1,2,3)
        IF(NOD.EQ.27) CALL FM27N3 (DER,IDER,FUN,SAMP, 3,1,2,3)
      ENDIF
C .... possibility of 4d elements etc. also.
      IF (FUN(1) .eq. 230964.e10) WRITE(*,'(A,I2,A,I3,A,I3,a,i2)')
     + '** WARNING: element not found: nodof=',nodof,' nod=',nod
     +,' type=',itype
      END
C***********************************************************************
c-----------------------------------------------------------------------
      SUBROUTINE WTHATN (NEN,NDIME,ITYPE,LCOORD,ILCOORD)
c
c     This recursively searches coords in the range -1 to +1
c     searching for the local coords of an element usng GSF.
c
      parameter( max_nn=30,max_level=12, idim=5,ider=idim,ismp=1,
     + tol = 1.e-5 )
      real lcoord(ilcoord,*)
      real fun(max_nn), der(idim,max_nn), smp(ismp,idim)
      integer point(5)
      do i=1,nen       ! null FUN incase GSF can't find it :-)
        fun(i) = 0.
      enddo
      inode = 0
      do n=1,max_level    ! loop the 'recursion' level 
c       print*,'n=',n
        do j=1,ndime
          point(j) = 0        ! set-up the 'initial values'
        enddo
C--------- loop through the local co-ords ------------------------------
      do level=1,99999
        if (level.gt.1) then
C-----------------------------------------------------------------------
c         ........ 'toggle' up the POINT pointers .............
    2     ipt = ndime
    1     point(ipt) = point(ipt) + 1
          if (point(ipt).gt.n)then
            point(ipt) = 0
            ipt = ipt - 1
            if (ipt.gt.0) goto 1    ! push left
            goto 99    ! finished all possiblilities
          endif
        endif
C .. Sieve of Erasmthus (sp?) search of the prime factors
      if (n.ge.2) then
        do j=2,max_level/2
          if (mod(n,j).eq.0)then
            iflag=0
            do ipt=1,ndime
              if (mod(point(ipt),j).ne.0) iflag=1
            enddo  
            if (iflag.eq.0) goto 2
          endif
        enddo
      endif
c     write(*,'(a,3i2)') 'point=',(point(j),j=1,ndime)
C-----------------------------------------------------------------------
c..... form the sampling point into SMP
        do j=1,ndime
          smp(1,j) = -1. + 2. * real(point(j)) / real(n) 
        enddo
        call gsf(ndime,nen,itype,der,ider,fun,smp,ismp,1)

c.......... loop the nodes to find which is = 1.0  ...................
        maxpos = 0
        nzero = 0
        do i=1,nen
          if ( abs(fun(i)-1.).lt.tol)then  ! 'own node'
            maxpos = i
          elseif (abs(fun(i)).lt.tol)then  ! 'other node'
            nzero = nzero +1
          endif
        enddo

        if (maxpos.ne.0.and.nzero+1.eq.nen) then      ! a 'hit'
          inode = inode + 1
          do j=1,ndime
            lcoord(maxpos,j) = smp(1,j)
          enddo
        endif
        if (inode.eq.nen)  return   ! all done

      enddo     ! end of the 'pseudo' local co-ord permutation loop
   99 continue
      enddo     ! end of the recusion 'level'
      print*,'*** ERROR: coords not found after',max_level,' iterations'
      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


C***********************************************************************
      SUBROUTINE MYGAUS(SMP,ISMP,WTS,NODOF,NGP)
C
C     a header routine to GAUSS which returns the quadrature
C     points as a single list (no implied looping) ie. in the same
C     form as NUMINT, NUMIN3 and QTURE
C     
C                                            Dan Kidger  April '91
C
C  .. need to add a flag to choose triangle/quad rules ...
C  .. and find a new name for 'mygauss'
C
      REAL SAMP(8,2), SMP(ISMP,*), WTS(*)
      INTEGER POINT(5)
      NGPI = NINT(REAL(NGP)**(1./REAL(NODOF)) )

      IF((NODOF.EQ.3).AND.
     + (NGP.EQ.6.OR.NGP.EQ.13.OR.NGP.EQ.14.OR.NGP.EQ.15)) THEN
        CALL QTURE(SMP,ISMP,WTS,NGP)
      ELSE IF (NGP.NE.NGP**NODOF)THEN

      DO 4,I=1,NODOF
    4 POINT(I)=1
      CALL GAUSS(SAMP,8,NGPI)
      DO 2,IC=1,NGP

      WTS(IC)=1.
      DO 3,I=1,NODOF
      SMP(IC,I) = SAMP(POINT(I),1)
    3 WTS(IC)   = SAMP(POINT(I),2) * WTS(IC)

      I=NODOF
    1 POINT(I) = POINT(I) + 1
      IF (POINT(I).GT.NGPI)THEN
       POINT(I)=1
       I = I - 1
       IF (I.GT.0) GOTO 1
      ENDIF
    2 CONTINUE
      NGP = IC -1
      ELSE
      PRINT*,' *** WARNING ',NGP,' Point Integration not understood'
      ENDIF
      END
      SUBROUTINE QTURE(SAMPL,ISAMPL,WT,NQP)
C
C         WEIGHTS AND SAMPLING POINTS FOR IRONS 3 D
C  13 point rule added 3/5/90    (so now 6,13,14,15 rules)
C
      REAL SAMPL(ISAMPL,*),WT(*)
      CALL NULL(SAMPL,ISAMPL,NQP,3)
      IF(NQP.EQ.6) THEN
      DO 1 I=1,NQP
    1 WT(I)=4./3.
      SAMPL(1,1)=-1.
      SAMPL(2,1)=1.
      SAMPL(3,2)=-1.
      SAMPL(4,2)=1.
      SAMPL(5,3)=-1.
      SAMPL(6,3)=1.
      ENDIF

      IF(NQP.EQ.13)THEN
      B=-0.49584802
      C= 0.88030430
      D= 0.79562143
      E= 0.025293237
      DO 21,I=1,6
   21 WT( I)= 0.54498736
      DO 22,I=7,12
   22 WT( I)= 0.507644216
      WT(13)= 1.68421056
      DO 23,I=1,6
      DO 23,J=1,3
      SAMPL(I,J)=B
   23 SAMPL(I,MOD(I-1,3)+1) =C
      DO 24,I=7,12
      DO 24,J=1,3
      SAMPL(I,J)= D
   24 SAMPL(I,MOD(I-1,3)+1) =E
      DO 25,J=1,3
      SAMPL(13,J)=0.
      DO 25,I=4,9
   25 SAMPL(I,J)=-SAMPL(I,J)
      ENDIF

      IF(NQP.EQ.14) THEN
      B=0.795822426
      C=0.758786911
      DO 2 I=1,6
    2 WT(I)=0.886426593
      DO 3 I=7,NQP
    3 WT(I)=0.335180055
      SAMPL(1,1)=-B
      SAMPL(2,1)=B
      SAMPL(3,2)=-B
      SAMPL(4,2)=B
      SAMPL(5,3)=-B
      SAMPL(6,3)=B
      DO 4 I=7,NQP
      DO 4 J=1,3
    4 SAMPL(I,J)=C
      SAMPL(7,1)=-C
      SAMPL(7,2)=-C
      SAMPL(7,3)=-C
      SAMPL(8,2)=-C
      SAMPL(8,3)=-C
      SAMPL(9,1)=-C
      SAMPL(9,3)=-C
      SAMPL(10,3)=-C
      SAMPL(11,1)=-C
      SAMPL(11,2)=-C
      SAMPL(12,2)=-C
      SAMPL(13,1)=-C
      ENDIF
      IF(NQP.EQ.15) THEN
      B=1.
      C=0.674199862
      WT(1)=1.564444444
      DO 7 I=2,7
    7 WT(I)=0.355555556
      DO 8 I=8,15
    8 WT(I)=0.537777778
      SAMPL(2,1)=-B
      SAMPL(3,1)=B
      SAMPL(4,2)=-B
      SAMPL(5,2)=B
      SAMPL(6,3)=-B
      SAMPL(7,3)=B
      DO 9 I=8,NQP
      DO 9 J=1,3
    9 SAMPL(I,J)=C
      SAMPL(8,1)=-C
      SAMPL(8,2)=-C
      SAMPL(8,3)=-C
      SAMPL(9,2)=-C
      SAMPL(9,3)=-C
      SAMPL(10,1)=-C
      SAMPL(10,3)=-C
      SAMPL(11,3)=-C
      SAMPL(12,1)=-C
      SAMPL(12,2)=-C
      SAMPL(13,2)=-C
      SAMPL(14,1)=-C
      ENDIF
      RETURN
      END
C************************************************************************
      SUBROUTINE FMTR10(DER,IDER,FUN,SAMP,ISAMP,I)
C
C     SHAPE FUNCS AND DERIVS - 10 NODE TRIANGLE
C       ... created using REDUCE as this wasn't in FE5LIB
C            DJK MAY 1990
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

C************************************************************************
      SUBROUTINE FMQM6(DER,IDER,FUN,SAMP,ISAMP,I,J,K)
C
C      this forms the SF's and Derivs for the 11 node brick
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)

      CALL FMLIN3 (DER,IDER,FUN,SAMP, ISAMP,I,J,K)

      ETA = SAMP(I,1)
      XI  = SAMP(J,1)
      ZETA= SAMP(K,1)
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
C************************************************************************
      SUBROUTINE FMBQ6(BEE,IBEE,DERIV,IDERIV,NOD)
C
C      THIS FORMS THE 3-D STRAIN-DISPLACEMENT MATRIX
C       ..... for an 8 node brick with 3 'nodeless' functions
C       .....  from IMS   July '90 -- 'the 11 node brick'
C       .....  this is the same as FORMB3 with 9 lines added
C
C     .... can't I sub-call FORMB3 or pass IDOF to FORM_BEE  :-)
C     --> isn't this FORMB3 with NOD=11 ???
C
      REAL BEE(IBEE,*),DERIV(IDERIV,*)
      DO 1 M=1,NOD
      N=3*M
      K=N-1
      L=K-1
      X=DERIV(1,M)
      BEE(1,L)=X
      BEE(4,K)=X
      BEE(6,N)=X
      Y=DERIV(2,M)
      BEE(2,K)=Y
      BEE(4,L)=Y
      BEE(5,N)=Y
      Z=DERIV(3,M)
      BEE(3,N)=Z
      BEE(5,K)=Z
      BEE(6,L)=Z
    1 CONTINUE
      BEE(1,25)=DERIV(1,9)
      BEE(2,29)=DERIV(2,10)
      BEE(3,33)=DERIV(3,11)
      BEE(4,26)=DERIV(2,10)
      BEE(4,28)=DERIV(1,9)
      BEE(5,30)=DERIV(3,11)
      BEE(5,32)=DERIV(2,10)
      BEE(6,27)=DERIV(3,11)
      BEE(6,31)=DERIV(1,9)
      RETURN
      END
C********************************************************************
C*****
C*****   g e o m e t r y   g e n e r a t i n g   r o u t i n e s 
C*****
C********************************************************************

      SUBROUTINE SLOGET(IP,IQ,IS,NXS,NYS,NYE,NZE,TOP,BOT,DEPTH,BRDTH
     +       ,COORD,ICOORD,NUM,NEN)
C
C    header to SLOGE3 (20nb), SLOG14 (14nb) and SLOG08 (8nb/11nb)
C
C    .. needs fixing so that it can call SLOGE2 also 
C
      REAL TOP(*),BOT(*),DEPTH(*),BRDTH(*),COORD(ICOORD,*)
      INTEGER NUM(*)
      IF(NEN.EQ.20)THEN
      CALL SLOGE3(IP,IQ,IS,NXS,NYS,NYE,NZE,TOP,BOT,DEPTH,BRDTH
     +       ,COORD,ICOORD,NUM)
      ELSE IF(NEN.EQ.14)THEN
      CALL SLOG14(IP,IQ,IS,NXS,NYS,NYE,NZE,TOP,BOT,DEPTH,BRDTH
     +       ,COORD,ICOORD,NUM)
      ELSE IF(NEN.EQ. 8 .OR. NEN.EQ.11)THEN
      CALL SLOG08(IP,IQ,IS,NXS,NYS,NYE,NZE,TOP,BOT,DEPTH,BRDTH
     +       ,COORD,ICOORD,NUM)
      IF(NEN.EQ.11)THEN
        DO 2,I=9,11
    2   NUM(I)=I
        ENDIF
      ELSE
      PRINT*,'  OOPS NOT A 8,11,14 OR 20 NODE BRICK'
      STOP
      ENDIF
      RETURN
      END
C********************************************************************

      SUBROUTINE SLOGE2(IP,IQ,NXS,NYS,NYE,TOP,BOT,DEPTH,COORD,ICOORD
     +,NUM)
C
C      THIS  FORMS THE COORDINATES AND STEERING VECTOR
C      FOR 8-NODE QUADRILATERALS IN A 'SLOPE' GEOMETRY
C      (NUMBERING IN THE Y-DIRECTION)
C      ---> corrected 23/6/88 by DJK to get the X-values right!   8-)
C      ...... ammended 7-9-88 by DJK
C
      REAL TOP(*),BOT(*),COORD(ICOORD,*),DEPTH(*)
      INTEGER NUM(*)
      NUM(1)=(IP-1)*(3*NYE+2)+2*IQ+1
      NUM(4)=(IP-1)*(3*NYE+2)+2*NYE+IQ+1
      NUM(5)=IP*(3*NYE+2)+2*IQ-1
      IF(IP.GT.NXS)THEN
      NUM(1)=NUM(1)-(((IP-1-NXS)*3)+0)*NYS
      NUM(4)=NUM(4)-(((IP-1-NXS)*3)+1)*NYS
      NUM(5)=NUM(5)-(((IP-1-NXS)*3)+3)*NYS
      ENDIF
      NUM(2)=NUM(1)-1
      NUM(3)=NUM(1)-2
      NUM(6)=NUM(5)+1
      NUM(7)=NUM(5)+2
      NUM(8)=NUM(4)+1
c ... lines to form 'G' deleteed
      IF(IQ.LE.NYS)THEN
      FAC1=(BOT(IP  )-TOP(IP  ))/ (DEPTH(NYS+1)-DEPTH(1))
      FAC2=(BOT(IP+1)-TOP(IP+1))/ (DEPTH(NYS+1)-DEPTH(1))
      COORD(1,1)=TOP(IP  )+(DEPTH(IQ+1)-DEPTH(1))*FAC1
      COORD(3,1)=TOP(IP  )+(DEPTH(IQ  )-DEPTH(1))*FAC1
      COORD(5,1)=TOP(IP+1)+(DEPTH(IQ  )-DEPTH(1))*FAC2
      COORD(7,1)=TOP(IP+1)+(DEPTH(IQ+1)-DEPTH(1))*FAC2
      ELSE
      COORD(1,1)=BOT(IP)
      COORD(3,1)=BOT(IP)
      COORD(5,1)=BOT(IP+1)
      COORD(7,1)=BOT(IP+1)
      ENDIF
      COORD(2,1)=.5*(COORD(1,1)+COORD(3,1))
      COORD(6,1)=.5*(COORD(5,1)+COORD(7,1))
      COORD(4,1)=.5*(COORD(3,1)+COORD(5,1))
      COORD(8,1)=.5*(COORD(7,1)+COORD(1,1))
      COORD(1,2)=DEPTH(IQ+1)
      COORD(8,2)=DEPTH(IQ+1)
      COORD(7,2)=DEPTH(IQ+1)
      COORD(3,2)=DEPTH(IQ)
      COORD(4,2)=DEPTH(IQ)
      COORD(5,2)=DEPTH(IQ)
      COORD(2,2)=.5*(COORD(1,2)+COORD(3,2))
      COORD(6,2)=.5*(COORD(5,2)+COORD(7,2))
      RETURN
      END
C************************************************************************
      SUBROUTINE SLOGE3(IP,IQ,IS,NXS,NYS,NYE,NZE,TOP,BOT,DEPTH,BRDTH
     +       ,COORD,ICOORD,NUM)
C
C      THIS FORMS THE COORDINATES AND STEERING VECTOR
C      FOR 20-NODE BRICKS IN A 'SLOPE' GEOMETRY
C      (NUMBERING IN THE Z-Y-X -DIRECTIONS)
C      ---> corrected 23/6/88 by DJK to get the X-values right!   8-)
C      ...... ammended 7-9-88 by DJK for slopes on a base 8->
C      -----> extension of SLOGE2 by DJK ....  10-2-89
C
      REAL TOP(*),BOT(*),COORD(ICOORD,*),DEPTH(*),BRDTH(*)
      INTEGER NUM(*),IFROM(12),ITO(12),IE(3,12)
      DATA IFROM/15,16,17,10,11, 3, 4, 5,14,18, 2, 6/
      DATA ITO  /14,20,18, 9,12, 2, 8, 6,13,19, 1, 7/
      DATA IE   / 1, 2, 3, 3, 4, 5, 5, 6, 7, 7, 8, 1,
     +           13,14,15,15,16,17,17,18,19,19,20,13,
     +            1, 9,13, 3,10,15, 5,11,17, 7,12,19/
      IF(IQ.LE.NYS.AND.IP.GT.NXS)THEN
      PRINT*,'OUT OF MESH RANGE..... ie this element aint ''ere'
      STOP
      ENDIF
      IF(IP.LE.NXS)THEN
      BASE=(IP-1)*(4*NYE*NZE+3*(NYE+NZE)+2)
      NYEH=NYE
      IQ1=IQ
      ELSE
      NYEH=NYE-NYS
      IQ1=IQ-NYS
      BASE = 4*NXS*NYE*NZE+3*(NXS*(NYE+NZE)+NYE*NZE)+2*(NXS+NYE+NZE)+1
     + - 3*NZE*NYEH -2*(NZE+NYEH) -1
     + + (IP-NXS-1)*(4*NYEH*NZE+3*(NYEH+NZE)+2)
      ENDIF
      N1=3*NZE*NYEH+2*(NZE+NYEH)+1
      N2=N1+(NZE+1)*(NYEH+1)
      NUM(15)=BASE+ (3*NZE+2)*(IQ1-1)+ 2*(IS-1)+ 1
      NUM(10)=NUM(15)+ (2*NZE+1)- (IS-1)
      NUM( 3)=NUM(15)+ (3*NZE+2)
      NUM(16)=BASE+N1+ (NZE+1)*(IQ1-1)+  (IS-1)+ 1
      NUM( 4)=NUM(16)+ (NZE+1)
      NUM(17)=NUM(15)+ N2
      NUM(11)=NUM(10)+ N2
      NUM( 5)=NUM( 3)+ N2
      DO 5,I=1,12
    5 NUM(ITO(I))=NUM(IFROM(I)) + 1
C     WRITE(*,'(20I4)')(NUM(I),I=1,20)
      IF(IQ.LE.NYS)THEN
      FAC1=(BOT(IP  )-TOP(IP  ))/ (DEPTH(NYS+1)-DEPTH(1))
      FAC2=(BOT(IP+1)-TOP(IP+1))/ (DEPTH(NYS+1)-DEPTH(1))
      COORD( 3,1)=TOP(IP  )+(DEPTH(IQ+1)-DEPTH(1))*FAC1
      COORD(15,1)=TOP(IP  )+(DEPTH(IQ  )-DEPTH(1))*FAC1
      COORD(17,1)=TOP(IP+1)+(DEPTH(IQ  )-DEPTH(1))*FAC2
      COORD( 5,1)=TOP(IP+1)+(DEPTH(IQ+1)-DEPTH(1))*FAC2
      ELSE
      COORD( 3,1)=BOT(IP)
      COORD(15,1)=BOT(IP)
      COORD(17,1)=BOT(IP+1)
      COORD( 5,1)=BOT(IP+1)
      ENDIF
      COORD( 3,2)=DEPTH(IQ+1)
      COORD(15,2)=DEPTH(IQ)
      COORD(15,3)=BRDTH(IS  )
      COORD(13,3)=BRDTH(IS+1)
      COORD(13,1)=COORD(15,1)
      COORD(19,1)=COORD(17,1)
      COORD( 1,1)=COORD( 3,1)
      COORD( 7,1)=COORD( 5,1)
      COORD(17,2)=COORD(15,2)
      COORD(13,2)=COORD(15,2)
      COORD(19,2)=COORD(15,2)
      COORD( 5,2)=COORD( 3,2)
      COORD( 1,2)=COORD( 3,2)
      COORD( 7,2)=COORD( 3,2)
      COORD(17,3)=COORD(15,3)
      COORD( 3,3)=COORD(15,3)
      COORD( 5,3)=COORD(15,3)
      COORD(19,3)=COORD(13,3)
      COORD( 1,3)=COORD(13,3)
      COORD( 7,3)=COORD(13,3)
      DO 8,I=1,12
      DO 8,J=1,3
    8 COORD(IE(2,I),J)=.5*(COORD(IE(1,I),J)+COORD(IE(3,I),J))
      RETURN
      END
C************************************************************************
      SUBROUTINE SLOG08 (IP,IQ,IS,NXS,NYS,NYE,NZE,TOP,BOT,DEPTH,BRDTH
     +       ,COORD,ICOORD,NUM)
C
C      THIS FORMS THE COORDINATES AND STEERING VECTOR
C      FOR  8-NODE BRICKS IN A 'SLOPE' GEOMETRY
C      (NUMBERING IN THE Z-Y-X -DIRECTIONS)
C
C     --> based on SLOG14 (cf)   . . . . . . .   11-Aug-90  by DJK
C
      REAL TOP(*),BOT(*),COORD(ICOORD,*),DEPTH(*),BRDTH(*)
      INTEGER NUM(*)
     + ,BASE,F1,F2
      IF(IQ.LE.NYS.AND.IP.GT.NXS)THEN
        PRINT*,'OUT OF MESH RANGE..... ie this element aint ''ere'
        STOP
      ENDIF
      IF(IP.LE.NXS)NYEH=NYE
      IF(IP.GT.NXS)NYEH=NYE-NYS
      F1 = (NYEH+1) * (NZE+1)
      F2= (NYE +1) * (NZE+1)
      IF(IP.LE.NXS)THEN
        IQ1=IQ
      BASE= (IP-1) * F2
      ELSE
        IQ1=IQ-NYS
        BASE = NXS* F2  + (IP-1-NXS) * F1
      ENDIF
C   ------- simply step around the element,to number it-------
C     PRINT*,BASE,F1,F2
      NUM( 6)=BASE+ (NZE+1) * (IQ-1) + IS
      NUM( 5)=NUM( 6)+ 1
      NUM( 2)=NUM( 5)+ NZE
      NUM( 1)=NUM( 2)+ 1
      NUM( 7)=NUM( 6)+ F1
      NUM( 8)=NUM( 7)+ 1
      NUM( 3)=NUM( 8)+ NZE
      NUM( 4)=NUM( 3)+ 1
C     WRITE(*,'(8I4)')(NUM(I),I=1,8)
C -------- set the 4 x-values of the 'trapezoidal' faces -----
      IF(IQ.LE.NYS)THEN
        FAC1=(BOT(IP  )-TOP(IP  ))/ (DEPTH(NYS+1)-DEPTH(1))
        FAC2=(BOT(IP+1)-TOP(IP+1))/ (DEPTH(NYS+1)-DEPTH(1))
        COORD( 1,1)=TOP(IP  )+(DEPTH(IQ+1)-DEPTH(1))*FAC1
        COORD( 5,1)=TOP(IP  )+(DEPTH(IQ  )-DEPTH(1))*FAC1
        COORD( 8,1)=TOP(IP+1)+(DEPTH(IQ  )-DEPTH(1))*FAC2
        COORD( 4,1)=TOP(IP+1)+(DEPTH(IQ+1)-DEPTH(1))*FAC2
      ELSE
        COORD( 1,1)=BOT(IP)
        COORD( 5,1)=BOT(IP)
        COORD( 8,1)=BOT(IP+1)
        COORD( 4,1)=BOT(IP+1)
      ENDIF
      COORD( 1,2)=DEPTH(IQ+1)
      COORD( 5,2)=DEPTH(IQ  )
      COORD( 1,3)=BRDTH(IS+1)
      COORD( 2,3)=BRDTH(IS  )
C --------  copy  to nodes of the same x,y,z ------
      DO 4,I=1,5,4
      COORD( 8-I,1)=COORD(9-I,1)
    4 COORD( 1+I,1)=COORD(I  ,1)
      DO 5,I=2,4
      COORD(I  ,2)=COORD(1,2)
    5 COORD(I+4,2)=COORD(5,2)
      COORD(4,3)=COORD(1,3)
      COORD(5,3)=COORD(1,3)
      COORD(8,3)=COORD(1,3)
      COORD(3,3)=COORD(2,3)
      COORD(6,3)=COORD(2,3)
      COORD(7,3)=COORD(2,3)
      RETURN
      END
C************************************************************************
      SUBROUTINE SLOG14 (IP,IQ,IS,NXS,NYS,NYE,NZE,TOP,BOT,DEPTH,BRDTH
     +       ,COORD,ICOORD,NUM)
C
C      THIS FORMS THE COORDINATES AND STEERING VECTOR
C      FOR 14-NODE BRICKS IN A 'SLOPE' GEOMETRY
C      (NUMBERING IN THE Z-Y-X -DIRECTIONS)
C
C     --> based on SLOGE3 (cf)   . . . . . . .   31-July-89  by DJK
C
      REAL TOP(*),BOT(*),COORD(ICOORD,*),DEPTH(*),BRDTH(*)
      INTEGER NUM(*) ,IFN(6),IF(6,4)
     + ,BASE,F1,F2
      DATA IFN /  9, 10, 11, 12, 13, 14/
     +    ,IF  /  1,  5,  1,  2,  1,  3,
     +            2,  6,  3,  4,  2,  4,
     +            3,  7,  5,  6,  5,  7,
     +            4,  8,  7,  8,  6,  8/
      IF(IQ.LE.NYS.AND.IP.GT.NXS)THEN
        PRINT*,'OUT OF MESH RANGE..... ie this element aint ''ere'
        STOP
      ENDIF
      IF(IP.LE.NXS)NYEH=NYE
      IF(IP.GT.NXS)NYEH=NYE-NYS
C     PRINT*,NYE,NZE
      F1=2*NYEH*NZE + NYEH+ NZE
      F2=2*NYE *NZE + NYE + NZE
      IF(IP.LE.NXS)THEN
        IQ1=IQ
        BASE=(IP-1)*(2*F2+1)
      ELSE
        IQ1=IQ-NYS
        BASE = NXS* (2*F2+1) +NYS* (2*NZE+1) + (IP-1-NXS) * (2*F1+1)
      ENDIF
C   ------- simply step around the element,to number it-------
C     PRINT*,BASE,F1,F2
      NUM( 6)=BASE+ (2*NZE+1)*(IQ1-1)+ IS
      NUM( 2)=NUM( 6)+ 1
      NUM(13)=NUM( 2)+ NZE
      NUM( 5)=NUM(13)+ NZE
      NUM( 1)=NUM( 5)+ 1
      NUM(12)=NUM( 2)+ F1
      NUM(10)=NUM(12)+ NZE
      NUM( 9)=NUM(10)+ 1
      NUM(11)=NUM( 9)+ NZE
      NUM( 8)=NUM(12)+ F1
      NUM( 4)=NUM( 8)+ 1
      NUM(14)=NUM( 4)+ NZE
      NUM( 7)=NUM(14)+ NZE
      NUM( 3)=NUM( 7)+ 1
C -------- set the 4 x-values of the 'trapezoidal' faces -----
      IF(IQ.LE.NYS)THEN
        FAC1=(BOT(IP  )-TOP(IP  ))/ (DEPTH(NYS+1)-DEPTH(1))
        FAC2=(BOT(IP+1)-TOP(IP+1))/ (DEPTH(NYS+1)-DEPTH(1))
        COORD( 2,1)=TOP(IP  )+(DEPTH(IQ+1)-DEPTH(1))*FAC1
        COORD( 1,1)=TOP(IP  )+(DEPTH(IQ  )-DEPTH(1))*FAC1
        COORD( 4,1)=TOP(IP+1)+(DEPTH(IQ  )-DEPTH(1))*FAC2
        COORD( 3,1)=TOP(IP+1)+(DEPTH(IQ+1)-DEPTH(1))*FAC2
      ELSE
        COORD( 2,1)=BOT(IP)
        COORD( 1,1)=BOT(IP)
        COORD( 4,1)=BOT(IP+1)
        COORD( 3,1)=BOT(IP+1)
      ENDIF
      COORD( 1,2)=DEPTH(IQ+1)
      COORD( 2,2)=DEPTH(IQ  )
      COORD( 1,3)=BRDTH(IS+1)
      COORD( 5,3)=BRDTH(IS  )
C --------  copy  to nodes of the same x,y,z ------
      DO 4,I=5,8
    4 COORD( I,1)=COORD(I-4,1)
      COORD( 3,2)=COORD( 1,2)
      COORD( 5,2)=COORD( 1,2)
      COORD( 7,2)=COORD( 1,2)
      COORD( 4,2)=COORD( 2,2)
      COORD( 6,2)=COORD( 2,2)
      COORD( 8,2)=COORD( 2,2)
      COORD( 2,3)=COORD( 1,3)
      COORD( 3,3)=COORD( 1,3)
      COORD( 4,3)=COORD( 1,3)
      COORD( 6,3)=COORD( 5,3)
      COORD( 7,3)=COORD( 5,3)
      COORD( 8,3)=COORD( 5,3)
C  ------------- put face nodes at the 'mid-face' ----------
      DO 8,I=1,6
        DO 8,J=1,3
        SUM=0
          DO 9,K=1,4
    9     SUM=SUM+COORD(IF(I,K),J)
    8   COORD(IFN(I),J)=SUM/4.
C     DO 10,I=1,14
C  10 WRITE(*,'(3F10.2)')(COORD(I,J),J=1,3)
      RETURN
      END
C************************************************************************
C*****
C*****  m y    b i t s   a n d   p i e c e s 
C*****
C********************************************************************

      SUBROUTINE FILES
C
C           to set up I/O files from file 'FILE'
C
      CHARACTER FILE*15
      OPEN(11,FILE='FILE')
      READ(11,'(A)')FILE
      CLOSE (11)
      LEN=INDEX(FILE,' ')-1
      OPEN(11,FILE=FILE)
      OPEN(11,FILE=FILE(1:LEN)//'.DAT')
      OPEN(12,FILE=FILE(1:LEN)//'.RES')
      OPEN(13,FILE=FILE(1:LEN)//'.PL')
      RETURN
      END
C************************************************************************
      SUBROUTINE AUTOFOS(IY,FOS,ITERS,SCALEITERS,STEPL)
C
C     at last .. an auto-fos tracing routine
C
      REAL FOSL(3,60)
      SAVE
      FOSL(1,IY)=FOS
      FOSL(2,IY)=ITERS
      IF(IY.EQ.1)THEN
        FOS=FOS+.001
        FOSL(3,1)=0
        FOSL(3,2)=0
      ELSE
         WRITE(*,'(F10.6,F10.0,F12.4)')(FOSL(I,IY),I=1,3)
         DX= FOS -FOSL(1,IY-1)
         DY=REAL(ITERS)-FOSL(2,IY-1)
         DY=DY/SCALEITERS
         VLEN=SQRT(DX*DX+DY*DY)
C        WRITE(8,*)FOS,STEPL,DX,VLEN
         FOS=FOS + STEPL*DX / VLEN
        FOSL(3,IY+1)= ITERS + STEPL*DY / VLEN * SCALEITERS
C       WRITE(8,*)DX,DY,FOS,ITERS
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FM20N3(DER,IDER,FUN,SAMP,ISAMP,I,J,K)
C
C     SHAPE FUNCS AND DERIVS - 20 NODE BRICK
C            DJK MARCH  1992
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      X = SAMP(I,1)
      Y = SAMP(J,1)
      Z = SAMP(K,1)
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
      SUBROUTINE FM26N3(DER,IDER,FUN,SAMP,ISAMP,I,J,K)
C
C     SHAPE FUNCS AND DERIVS - 26 NODE BRICK
C            DJK MARCH  1992
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      X = SAMP(I,1)
      Y = SAMP(J,1)
      Z = SAMP(K,1)
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
      SUBROUTINE FM27N3(DER,IDER,FUN,SAMP,ISAMP,I,J,K)
C
C     SHAPE FUNCS AND DERIVS - 27 NODE BRICK
C            DJK MARCH  1992
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      X = SAMP(I,1)
      Y = SAMP(J,1)
      Z = SAMP(K,1)
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
      SUBROUTINE FM05N2(DER,IDER,FUN,SAMP,ISAMP,I,J)
C
C     SHAPE FUNCS AND DERIVS - 5 NODE quad
C         by  Dan Kidger         13th July 1992
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      X = SAMP(I,1)
      Y = SAMP(J,1)
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
      SUBROUTINE FM12N2(DER,IDER,FUN,SAMP,ISAMP,I,J)
C
C     SHAPE FUNCTIONS AND DERIVS - 12 node quad
C         by  Dan Kidger         8th April 1992
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      X = SAMP(I,1)
      Y = SAMP(J,1)
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
      SUBROUTINE FM17N2(DER,IDER,FUN,SAMP,ISAMP,I,J)
C
C     SHAPE FUNCS AND DERIVS - 26 NODE BRICK
C            DJK APRIL  1992
C
      REAL DER(IDER,*),FUN(*),SAMP(ISAMP,*)
      X = SAMP(I,1)
      Y = SAMP(J,1)
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
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------


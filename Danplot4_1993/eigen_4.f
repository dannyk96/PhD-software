      OPTIONS ( FULLCHECK,UNDEF,DREAL, INTL)
      PROGRAM EIGEN
C---------------------------------------------------------------------------
C  general purpose EIGENVALUE/VECTOR of F.E. stiffness matrrices program
C---------------------------------------------------------------------------
C
C       c. Dan.Kidger   April 1992
C
C  extensively developed during my Ph.D work (1986-1990)
C  this is the tidied up version (splitting PLOTMS off)
C  ( the subroutines may be found in L.1.F77, l.3.f77 and FE5LIB )
C
C  latest changes --> pushing out GSF and MYGAUSS subs.
C ---> LINPACk Eigenvalue routines to replace the NAg. F02ABF
C ---> 7/10/91 'column' output restored
C ---> 9/04/92 WTHATN interface to directy get the local co-ords
C

C    ... Max. problem size parameters : MN=nodes, MDOF=Dof., MDF=space,
C   ... NS=strains, ISAMP=G.P.rule order
      PARAMETER ( MN = 27, MDOF = 3*MN )
      PARAMETER ( MDF = 3, NS   = 6  )
      PARAMETER (ISAMP= 6, ISMP = 64 )

C   put 'BIG' arrays in COMMON blocks
c      COMMON/COM/KM

C -------------------- standard FE arrays ---------------------
      REAL DEE(NS,NS), BEE(NS+1,MDOF), KM(MDOF,MDOF),
     +     FUN(MN), DER(MDF,MN), DERIV(MDF,MN), COORD(MN,MDF),
     +     JAC(MDF,MDF), JAC1(MDF,MDF),
     +     SMP(ISMP,3),  WTS(ISMP),
     +     V(MDOF,MDOF), R(MDOF)
C .. next is obsolete !! ???
C     REAL   BT(MDOF,NS), DBEE(NS,MDOF), BTDB(MDOF,MDOF)
      INTEGER NGPR(2), U1,U2,U3
      LOGICAL L1D,L2D,L3D, LTRI, LRTRI
      CHARACTER FN*15, RESP*5, FRMAT1*12

C ------------- double precision for NAg routines -----------------
c .... could just RUN the prog. in double precision ??
      REAL*8    KMD(MDOF,MDOF),RD(MDOF),VD(MDOF,MDOF) ,WORK(MDOF)
      INTEGER*4 IFAIL, MDOFD, IDOFD
C ----------- 'standard' FEPROGS parameters -----------------------
      DATA IKM/MDOF/, IDEE/NS/,IBEE/NS+1/
      DATA IJAC,IJAC1,IDER,IDERIV/4*MDF/,ICOORD/MN/

      DATA Q/1.00/ ,RESIZE/1.00/
      DATA U1,U2,U3 /11,12,13/

C ============= M A I N    P R O G R A M =============================
C ... better to write as a table !

      WRITE (  *,'( 79(''*''),//,A,//,A,//,79(''*'') //)'    )
     +'    Welcome to Dan''s Wonderful Eigenmode Program',
     +'      >> just follow the instructions >> {2} = default option'
C----------------- new interface -------------------------
  992  PRINT*,'Enter the element code as NEN,NODOF,ITYPE'
     +      //'.... (type HELP for a list)'
      READ(*,*,iostat=iostat) NEN,NODOF,ITYPE
      IF (iostat.ne.0)then
      WRITE(*,'(T12,A)') 
     +'   Elements           No. of Nodes         Integration rule    ',
     +'   Available              NEN                   NGP            ',
     +' =============       ==============        ================    ',
     +'                                                               ',
     +'  1D: lines          1...2...3...4           1  2  3  4        ',
     +'                                                               ',
     +'  2D: triangles:     3...6..10..15           1  3  4 ..        ',
     +'      quads:         4...8...9..12..17       1  4  9 16 25     ',
     +'                                                               ',
     +'  3D; tetrahedra:    4                       1  4              ',
     +'      bricks:        8..11..14..20..26..27   1  6  8 14 15 27  ',
     +' '                                
        GOTO 992
      ENDIF
      L1D = (NODOF.EQ.1)
      L2D = (NODOF.EQ.2)
      L3D = (NODOF.EQ.3)
        IF(L1D) IH = 1
        IF(L2D) IH = 3
        IF(L3D) IH = 6
      IDOF = NODOF*NEN

      LTRI = (NEN.LT.0)
      LRTRI= (L2D.AND.(NEN.EQ.3.OR.NEN.EQ.6.OR.NEN.EQ.10.OR.NEN.EQ.15))
C .. ( LRTRI if a real triangle element , LTRI if a distorted quad.)

C -------------------- set integration rules -------------------------
      IRULE = 1
        PRINT*,'Enter Integration rule type:           {1}',
     + '1=Full Int,  2 =S.R.I (K,G), 4=S.R.I (Lame,G)'//
     + '5 = axisym, 6= plane-stress'
        READ*,IRULE
        NGPASS=1
        IF(IRULE.EQ.2.OR.IRULE.EQ.4) NGPASS=2 
      PRINT*,'Enter Number of Quadrature Points  {4}'
      IF (L1D)                PRINT*,'= 1, 2, 3, 4, ... '
      IF (L2D.AND..NOT.LRTRI) PRINT*,'= 1, 4, 9,16, ...'
      IF (L2D.AND.     LRTRI) PRINT*,'= 1, 3, 4, 6, 7, 12, 16'
      IF (L3D.AND.NEN.NE.4)   PRINT*,'= 1, 8, 27  / 6, 13, 14, 15'
      IF (L3D.AND.NEN.EQ.4)   PRINT*,'= 1, 4, 5'

      READ*,(NGPR(I) ,I=1,NGPASS)
      IF (IRULE.EQ.1) NGPR(2) = NGPR(1)

C --------------- Input Poisson's Ratio etc. -------------------------  
      PRINT*,'ENTER Poisson''s ratio {0.3}      (-ve for looped )'
      READ*,PRL
      IF (PRL.LT.0) IOP = 2
      IF (PRL.GE.0) IOP = 1
      IF (IOP.EQ.1) THEN
        PRH = PRL + .00001
      ELSE
        PRL=0.
        PRH=.5-.00001
      ENDIF

C  .. can modify next to allow aspect ratio effect of ANY element
C      IF (L2D.AND.NEN.EQ.8)THEN
C        PRINT*,'ENTER ELEMENT DISTORTION FACTOR'
C       READ*,Q
C     ENDIF
C ---------------------- Open output files -----------------------------
      WRITE (FN,'(I2.2,A,I2.2,A,A)')  NEN,'nb',NGPR(1),'gp','.p'
        IF (NODOF.EQ.2) FN( 4: 4) = 'q'
        IF (IRULE.EQ.2) FN( 8: 9) = 'sa'
        IF (IRULE.EQ.4) FN( 8: 9) = 'sb'
      OPEN  (7,FILE='FILE')
      WRITE (7,'(A)')FN
      CLOSE (7)
      OPEN  (U2,FILE=FN)
      FN(11:11)='l'
      OPEN  (U3,file=fn)
C ********************* loop Poisson's ratio **********************
C ....... nice to be able to adjust the step length
C .... even to be able to input a new value each time 
C ..should really have a INTEGER loop counter
      DO 322,PR=PRL,PRH, .02
      WRITE(U2,'(A20,F8.4)')'POISSON''S RATIO=',PR

  432 CONTINUE
      CALL NULL(KM,IKM,IDOF,IDOF)

C --------------- LOOP THE FULL/SRI CASES ----------------------------
      DO 21 ,IGPASS=1, NGPASS
      CALL FRMD (DEE,IDEE,1.,PR,0.,NODOF,IRULE+IGPASS-1)

C ------- get Gauss point as a long list of ALL g.p's ----------------
C .. nice if MYGAUS provided the appropriate steering !
      NGP = NGPR(IGPASS)
      IF     (LRTRI) THEN
          CALL NUMINT (SMP,ISMP,WTS,NGP)
      ELSEIF (L3D .AND. NEN.EQ.4) THEN
          CALL NUMIN3 (SMP,ISMP,WTS,NGP)
      ELSE
        CALL MYGAUS (SMP,ISMP,WTS,NODOF,NGP)
C ... my header to GAUSS .. gives SMP like QTURE and NUMINT :-)
      ENDIF
C ------------------- get the element geometry -----------------------
C --> Note that I would rather use the local geometry directly
C --> rather than having to rely on a set of geometry routines
C --> which might not 'match' the right configuration.
C
c      CALL NULL (COORD,ICOORD,ICOORD,MDF)
c      IF (LRTRI)   CALL GEOTRI (                COORD,ICOORD,NEN)
c      IF (LTRI)    CALL TRQUAD (                COORD,ICOORD,NEN)
c      ENDIF
c      DO 5,I=1,IDOF
c    5 G(I)=I
c

      call wthatn (nen,nodof,itype, coord,icoord)

C  ------- normalise coords to centre at (0,0,0) edge length = 1 -----

C ... first find mean x (y,z) coord, then push to +/- 1/2
      DO J=1,NODOF
        SUM = 0.
        DO I=1,NEN
          SUM = SUM + COORD(I,J)
        ENDDO
        SUM = SUM / REAL(NEN)
        DO I=1,NEN
          COORD(I,J) = (COORD(I,J)-SUM) * RESIZE

c ... apply the 'q' fudge factor to 'lengthen' x and 'shorten 'y' !
c ... should have an IF for being not 1D and being the FIRST pass 
c.. NOTE HERE that I really ought to use modal paricipation factors
C             to get the 'match' not just assume they are 'close'

          IF (j.eq.1) coord(i,j) = coord(i,j) * q
          IF (j.eq.2) coord(i,j) = coord(i,j) / q
        ENDDO
      ENDDO

C .. ie. following only writes if 'pure' coords (pass 1 ?)
C ... therefore I should TEST that it is pass #1 !
      IF (ABS(Q-1.).LT. .00001 ) THEN
        WRITE (U2,'(A)')'NODAL CO-ORDS ARE :-'
C .. forget formats ! put a DO-loop around it !
C       WRITE (FRMAT1,'(A,I2,A)') '(',NEN,'F9.5)'
C       WRITE (U2,FRMAT1) ((COORD(I,J),I=1,NEN),J=1,NODOF)
        WRITE (FRMAT1,'(A,I2,A)') '(I3,',NODOF,'F9.5)'
        WRITE (U2,FRMAT1) (I,(COORD(I,J),J=1,NODOF),I=1,NEN)
      ENDIF

C --------------------- LOOP THE GAUSS POINTS ------------------------
      DO 20 I=1,NGP
        CALL GSF(NODOF,NEN,ITYPE,DER,IDER,FUN,SMP,ISMP,I)
C .... 'standard bit' .....
      IF(NEN.EQ.11)THEN
        CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,NODOF,8,NODOF)
      ELSE
        CALL MATMUL(DER,IDER,COORD,ICOORD,JAC,IJAC,NODOF,NEN,NODOF)
      ENDIF

C ---> note: I can now use 'invert'
      IF(L3D)CALL TREEX3(JAC,IJAC,JAC1,IJAC1,DET)
      IF(L2D)CALL TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)

      CALL MATMUL(JAC1,IJAC1,DER,IDER,DERIV,IDERIV,NODOF,NODOF,NEN)
C .. only need to null at the top !
      CALL NULL(BEE,IBEE,NS+1,IDOF)

C .. write and include my OWN FMB ! , FMBQ6 is the same !!??
      IF (L2D)       CALL FORMB  (BEE,IBEE,DERIV,IDERIV,NEN)
      IF (L3D)       CALL FORMB3 (BEE,IBEE,DERIV,IDERIV,NEN)
      IF (NEN.EQ.11) CALL FMBQ6  (BEE,IBEE,DERIV,IDERIV,NEN)

C <><><><><><><><><>   Patch to do BTDB directly   <><><><><><><><><>

   20 CALL FMBTDB (BEE,IBEE,DEE,IDEE,KM,IKM,DET*WTS(I),IH,IDOF)

C <><><><><><>><><><><><><><><><><><><><><><><><><><><><><>><><><>><>

   21 CONTINUE

C ------------------------ NAG  BIT ----------------------------------
            DO 40,I=1,IDOF
            DO 40 J=1,IDOF
   40       KMD(J,I)=KM(I,J)
            IFAIL = -1
            MDOFD = MDOF
            IDOFD = IDOF
C ...... NAG eigenvalues and vectors  
C            CALL F02ABF(KMD,MDOFD,IDOFD,RD,VD,MDOFD,WORK,IFAIL)
C ... pull double precision back to single prec.
C---> E I S P A C K   eigenvectors and eigenvalues 
C      DO 801,I=1,IDOFD
C      DO 802,J=1,IDOFD
C  802 VD(I,J)=0.
C  801 VD(I,I)=1.

      CALL TRED2 (MDOFD,IDOFD,KMD,RD,WORK, VD)
      CALL TQL2  (MDOFD,IDOFD,    RD,WORK, VD,IFAIL)

C---> E I S P A C K   eigenvectors and eigenvalues 

            DO 41,I = 1,IDOF
            IF (IOP.NE.-1) R(I) = RD(I)
            DO 41,J = 1,IDOF
   41       V(I,J)  = VD(I,J)
C -------------- tidy up the eigenvectors ----------------------------
        IF ( IOP.EQ.1)THEN
        IOP = -1
C *** CHANGE NEXT ***
C ... better to properly loop the Q factor 
C         ( NOTE not done for looped Poisson's Ratio ? )

        IF (ABS(Q-1.).LT..1) Q=1.01
        GOTO 432
        ENDIF
      r(idof+1) = 0
      CALL WRVAL (1,R,IDOF)

C  ------------ redefine unit translations & rotations ---------------
C .. better check is to test as a f(NODOF) !
      IF(.NOT.(L2D.AND. (ABS(R(4)).LT..001)   .OR.
     +         L3D.AND. (ABS(R(7)).LT..001))) THEN
        CALL REDEF(V,MDOF,COORD,MN,NEN,NODOF)

C ---- if 11NB then we need to UNSET the last three freedoms ! -------
c ... but only for the last 9 modes ! why ??
        IF (NEN.EQ.11)THEN
          DO 127,I=1,6
            DO 127,J=25,33
  127     V(J,I)=0.0
        ENDIF
      ENDIF
C ------------------- normalise Eigenvectors, so max = 1. ------------
      CALL NRMLIZ (V,MDOF,NEN,NODOF)
C -------------------- output Eigenvectors ---------------------------

      WRITE (U2,'(A)') 'Eigenvalues are......'
      WRITE (U2,'(8F10.5)')(R(I),I=1,IDOF)

      WRITE (U2,'(A)') 'Eigenvectors are......'
C      WRITE (FRMAT1,'(A,I2,A)') '(',NEN,  'F9.5)'
       WRITE (FRMAT1,'(A,I2,A)') '(',NODOF,'F9.5)'

      DO 51,K=1,IDOF
      WRITE (U2,  '(I9)') K
      WRITE (U2,  '(I9)') NEN
      WRITE (U2,FRMAT1)  ( V(J,K), J = 1,IDOF)
C     WRITE (U2,FRMAT1) 
C    +  ((JJ+nodof-1)/nodof,( V(J,K), J=JJ,JJ+NODOF-1), JJ=1,IDOF,NODOF)

   51 CONTINUE

C -------------- output .PL format file ------------------------------
      WRITE(U3,'(I5)') NODOF,NEN
      WRITE (FRMAT1,'(A,I2,A)') '(I3,',NODOF,'F9.5)'
      WRITE (U3,FRMAT1) (I,(COORD(I,J),J=1,NODOF),I=1,NEN)
      WRITE (U3,'(  I5)') 1
      WRITE (U3,'(30I3)') 1,NEN,(I,I=1,NEN),1
      DO K=1,IDOF
        WRITE(U3,'(A,I3,A,F10.5)')
     +          '# Eigen mode No.',K,' Eigenvalue=',R(K),' ',NEN
        WRITE (FRMAT1,'(A,I2,A)') '( I3,',NODOF,'F9.5)'
        WRITE (U3,FRMAT1) 
     +  ((JJ+nodof-1)/nodof,( V(J,K), J=JJ,JJ+NODOF-1), JJ=1,IDOF,NODOF)
      ENDDO

C --------------- optional eigenvector diagnostics -------------------
C --> such as check Radial/ Tangential components, moments, weights
C --> 
C      IF (L3D)  CALL RADTAN  (U2,COORD,ICOORD,V,IV,IDOF,NEN)
C      IF (L3D)  CALL MOMENT  (U2,COORD,ICOORD,V,R,IDOF,G)
C      IF (L3D)  CALL DIAGVEC (U2,V,MDOF,R,IDOF)

c       IF (IOP.EQ.-1) CLOSE(U2)
c       IF (IOP.EQ.-1) CLOSE(U3)

C --------------- end of Poisson's ratio loop ------------------------
  322 CONTINUE
      END

C*********************************************************************


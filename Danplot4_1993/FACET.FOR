C-----------------------------------------------------------------------
c...... These Subroutines are for:
c          x    Handling the FACETS.. 
c          x    and SAMPLE-ing to get interpolated values
c
c       Dan Kidger 15-6-93
c
C    20-8-93 : extended FACETS (3->7 cols) to include edge info
C-----------------------------------------------------------------------
      SUBROUTINE SAMPLE
     +  (SC,ISC,DISP,IDISP,NEN,NDIME,NDIM,VAL,VEC, POINTS,IPTS,IPT, IOP)
C
C     IN: COORD, DISP = nodal co-ords and disps of the element(face)
C         PNT = sampling point co-ord
C     this returns the co_ord, dIsp, strain,stress,etc.
C     of a sampled point within an element
C     local x,y in CO(),
C     returned single value in  VAL, {vector in VEC if appropriate}
C
C     IOP = 20   x,y,z coordinate
C     IOP = 14   .. normal to the face
C     IOP = 30+  x,y,z displacements
C     IOP = 40+  x,y,z strains
C     IOP = 50+  x,y,z stresses            etc.   etc.
C
C-------------- arguements -------------------
      REAL VAL, VEC(*), SC(ISC,*), DISP(IDISP,*),POINTS(IPTS,*)
      INTEGER NEN,NDIME,NDIM

C---------------------- FE5LIB interface ----------------------------
      PARAMETER (IDER=2,ISMP=2,IDERIV=3,IJAC=3,ICOORD=20)
      REAL DER(IDER,icoord), FUN(icoord), SMP(ISMP,3)
     +  ,COORD(ICOORD,3), JAC(IJAC,IJAC), DERIV(IDERIV,icoord)
     +  ,EPS(9), STRESS(9),DET , JAC1(IJAC,IJAC), DEE (3,3)
C---------------------------------------------
      IOP_2 = mod(iop,10)
C---------------------------------------------

      DO J=1,NDIME
        SMP(1,J) = POINTS(IPT,J)
      ENDDO
      DO J=1,NDIM
        DO I=1,NEN
          COORD(I,J) = SC(I,J)
        ENDDO
      ENDDO
C---------------------------------------------
c... OK first get the shape functions (and ders)
      CALL GSF (NDIME,NEN,1, DER,IDER,FUN,SMP,ISMP,1)

c-----------------------------------------------------------------------
      IF (IOP.EQ.10) THEN                 !-- the co-ords of the point
        DO J=1,NDIM
          VEC(J) = 0
          DO I=1,NEN
            VEC(J) = VEC(J) + FUN(I) * COORD(I,J)
          ENDDO
        ENDDO
        RETURN

c-----------------------------------------------------------------------
c.... normal vector to the face (X-product of 2 in-plane vectors)
      ELSEIF (IOP.EQ.14.or.iop.eq.15) THEN

        DO I=1,NDIME
          DO J=1,NDIM
            X=0.0
            DO K=1,NEN
              X=X+DER(I,K)*COORD(K,J)
            ENDDO
            JAC(I,J)=X
          ENDDO
        ENDDO
        VEC(1) = JAC(1,2)*JAC(2,3) - JAC(2,2)*JAC(1,3)
        VEC(2) = JAC(1,3)*JAC(2,1) - JAC(2,3)*JAC(1,1)
        VEC(3) = JAC(1,1)*JAC(2,2) - JAC(2,1)*JAC(1,2)
        RETURN

c-----------------------------------------------------------------------
      ELSEIF (IOP.GE.20.AND.IOP.LE.39) THEN      !-----coords/disps-----
        DO J=1,NDIM
          VEC(J) = 0
          DO I=1,NEN
            IF (IOP.le.29) THEN
              VEC(J) = VEC(J) + FUN(I) * SC  (I,J) ! coord
            ELSEIF (IOP.le.39) THEN
              VEC(J) = VEC(J) + FUN(I) * DISP(I,J) ! disp
            ENDIF
          ENDDO
        ENDDO
        IF(IOP_2.EQ.0) val = sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
        IF(IOP_2.EQ.1) val = vec(1)
        IF(IOP_2.EQ.2) val = vec(2)
        IF(IOP_2.EQ.3) val = vec(3)

c-----------------------------------------------------------------------
      ELSEIF (IOP.ge.40.and.iop.le.59) THEN  

      CALL MAT_MULT (DER,ider,COORD,icoord,JAC,ijac,2,NEN,2)
c     CALL INVERT (JAC,IJAC,DET,2L)
      DET = JAC(1,1)*JAC(2,2) - JAC(2,1) * JAC(1,2)
      IF (ABS(DET).LT. 1.E-10) THEN
        VAL =0.
        RETURN
      eNDIF

      CALL TWOBY2 (JAC,IJAC,JAC1,IJAC,DET)
      CALL MAT_MULT (JAC1,IJAC,DER,IDER,DERIV,IDERIV,2,2,NEN)

C      CALL NULL   (BEE,IBEE,4,2*NEN)
C      CALL FORMB  (BEE,IBEE, DERIV,IDERIV,NEN)
C      CALL MVMULT (BEE,IBEE, ELD,3,2*NEN,EPS)
      DO I=1,4
        EPS(I)=0.      !   OK form the 2D strains in EPS form DERIV
      ENDDO
      DO I=1,NEN 
        EPS(1) = EPS(1) + DERIV(1,I) * DISP(I,1)
        EPS(2) = EPS(2) + DERIV(2,I) * DISP(I,2)
        EPS(3) = EPS(3) + DERIV(2,I) * DISP(I,1)
     +                  + DERIV(1,I) * DISP(I,2)
      ENDDO

C--------------- Now get the 'elastic' stresses ------------------------
      IF (IOP.GE.50) THEN 
c       CALL FMDEPS (DEE,3,1.E6,0.3)    !-- is this plane-stress ?
        CALL FORM_DEE (DEE,3,1.e6,.3,0.,  2,1)

        CALL MVMULT (DEE,3, EPS,3,3, STRESS)
        DO J=1,3
          EPS(J) = STRESS(J)
        ENDDO
      ENDIF
C----------------------- get the principal values ----------------------
      
      IF (IOP_2.GE.6 .and. IOP_2.LE.8) THEN   ! if INVAR: put shear in 3
        IF (IOP.LT.50) EPS(3) = EPS(3) /2.          !-- divide by 2
        CALL INVARIENTS (EPS,2, EPS(6),EPS(7),EPS(8))
c----   CALL INVAR (EPS, sigm,  DSBAR ,THETA ) ------
      ELSEIF (IOP_2.EQ.9) THEN
        EPS(9) = DET
      ENDIF
      VAL = EPS (IOP_2)

c-----------------------------------------------------------------------
      ENDIF   ! end of IOP option_choice
      END

C-----------------------------------------------------------------------
      SUBROUTINE FSTRIP2 (NEL,NUMS,INUMS,FACETS,PL2,IFCETS,NFCETS,NDIME)
C
C     THIS turns NUMS into FACETS-a list of element face pointers
C     .. the sequel to FSTRIP.. but all FACETS are recorded
C
C     Data as:  Element #, Facet # (1-6), 'touching' facet #
C
C 19-11-93 hmm. could make 'touching' # -ve if of a different IMAT ??
C----------------------------------
C     A similar method should be possible for 'edges' to give a list
C     of element edges with pointers to their neighbours.
C     ( in 3D edges have more than one neighbour so this is invalid)
C  -- so <elem #, edge #, touching edge #>.
C
C     The clever 'bit' is to daisy-chain these pieces up into strings
C     Hence one polygon for each piece ?
C----------------------------------
C
      INTEGER
     +     NUMS (INUMS,*)        !- The element steerings
     +  ,FACETS (7,IFCETS)       !- 3->7 so can store edge info too?
     +     ,PL2 (IFCETS)         !- workspace only
     +     ,NUM (27)  ,NUM2(27)  !- steering for one element
     +     ,FS1 (19)  ,FS2 (19)  !- nodes on the first/second facet

c... this is used to make 'OFF' style files appear as simply polygons
c... (I could allow 3nt and 4nq's in .OFF ??)
c     ITYPE= 1   !--should be an arguement (eg. ITYPE=9)

C------------- patch for speed-up if NDIME=2D --------------------------
c.. NDIME is a property of an element NOT the whole mesh
      IF (NDIME.LT.3) THEN
        NFCETS = NEL
        DO I=1,NFCETS
          FACETS(1,I) = I   !--> = the element number
          FACETS(2,I) = 1   !--> only one face (so face # == 1)
          FACETS(3,I) = 0   !--> all facets on the 'outside' (no pointer)
        ENDDO
        RETURN
      ENDIF
C------------------- find the MIN node # of each facet -----------------
      IFC = 0
      DO IEL1=1,NEL
        CALL GET_ELEMENT 
     +  (NUMS,INUMS,IEL1, NUM, NOD1,NDIME,ITYPE, IMAT1,IU1,IU2,IU3)
c        NOD1  = NUMS (1,IEL1)
c        IMAT1 = NUMS (NOD1+2,IEL1)
C.. in the next *only* NN_F is returned for: option=1
        CALL GFACE (NOD1,NDIME,ITYPE,NFACES,FS1,NN_F1,NN_FT1, 1 )

        DO IFACE1=1,NFACES         ! loop this element's faces
c          print*,'iel1=',iel1,' imat1=',imat1,
c     +       ' iface1=', iface1,' ndime=',ndime
          CALL GFACE (NOD1,NDIME,ITYPE,IFACE1,FS1,NN_F1,NN_FT1, 2 ) 
          N_MIN= NUM (FS1(1))
          N_MAX= NUM (FS1(1))
          DO I=2,NN_F1
            N_MIN= MIN (N_MIN,NUM (FS1(I)))
            N_MAX= MAX (N_MAX,NUM (FS1(I)))
          ENDDO
          IFC = IFC + 1
             PL2(  IFC) = -N_MIN      ! ie. this facets' MAX and MIN
          FACETS(3,IFC) = -N_MAX      ! . node numbers
          FACETS(2,IFC) = IFACE1
          FACETS(1,IFC) = IEL1
        ENDDO
      ENDDO
      NFCETS = IFC
C---------- now loop facets and compare with the other facets ----------
      DO IFC1 = 1,NFCETS
        DO IFC2 = IFC1+1,NFCETS
          IF (PL2(IFC2).NE.PL2(IFC1)) GOTO 1        !- skip cos Mins not ==
          IF (FACETS(3,IFC2).NE.FACETS(3,IFC1)) GOTO 1   !- or Max's not ==
          IF (FACETS(3,IFC2).GT.0) GOTO 1         !- skip (already done)

C-------------- a possible 'hit' so test all interface nodes -----------
            IFACE1= FACETS(2,IFC1)
            IFACE2= FACETS(2,IFC2)

            IEL1  = FACETS(1,IFC1)
            IEL2  = FACETS(1,IFC2)

            CALL GET_ELEMENT 
     +     (NUMS,INUMS,IEL1, NUM, NOD1,NDIME,ITYPE, IMAT1,IU1,IU2,IU3)
            CALL GET_ELEMENT 
     +     (NUMS,INUMS,IEL2, NUM2,NOD2,NDIME,ITYPE, IMAT2,IU1,IU2,IU3)

            CALL GFACE (NOD1,NDIME,ITYPE,IFACE1,FS1,NN_F1,NN_FT1, 2 )
            CALL GFACE (NOD2,NDIME,ITYPE,IFACE2,FS2,NN_F2,NN_FT2, 2 )

            IF (NN_FT1.NE.NN_FT2) GOTO 1   ! skip if incompatible
            IBOT = NUM(FS1(1))
            DO I2=1,NN_F1        ! loop to try and find the 'base' node
              IF (NUM2(FS2(I2)) .EQ. IBOT ) GOTO 3
            ENDDO
            GOTO 1    ! no match so skip straight to the next facet

    3     IBASE = I2      ! so node I2==IBOT.. the 'base node'

          DO I1=2,NN_F1
            I2 = 1+MOD (NN_F1+IBASE -I1, NN_F1)     
            IF (NUM(FS1(I1)).ne.NUM2(FS2(I2)) ) GOTO 1
          ENDDO
C         ----- OK, a 'hit' so mark each of the 2 facets ------
          FACETS(3,IFC1) = IFC2
          FACETS(3,IFC2) = IFC1
          
   1      CONTINUE   ! 'skip-to' pointer :-)
        ENDDO
      ENDDO
C---------------- 'zap' all 'orphans' to zero :-) ----------------------
      DO I=1,NFCETS
        IF (FACETS(3,I).LT.0) FACETS(3,I) = 0
      ENDDO

      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GET_FAC (LC_SF,ILC_SF, I_SF, N_SF, NN_SF, NN_FT)
C
C     this returns the I_SF'th sub-facet of NN_SF points in LC_SF
C     of a N_SF*N_SF patch
C

C.... really the looping of the sub-facet nodes should be in the main
c.... program.. a sub-routine could return them (eg. 4 node or 8 node
C....                                             or even sub-triangles) 
C.... ¨ or have 2 interfaces to SAMPLE .. one will do all nodes
C....  the other just a given (by #?) point

      REAL LC_SF (ILC_SF,*)
      INTEGER LC4 (2,4)
      DATA LC4/-1,-1, -1,0,  0,0, 0,-1/  ! changed to clockwise 5-6-92 !
c--> 4nq (nice to have an 8nq also)
      K = NN_FT

C------------ quadrilateral subdivision --------------------------------
      IF (K.eq. 4.or.K.eq. 5.or.K.eq. 8.or.
     +    K.eq. 9.or.K.eq.12.or.K.eq.17     ) THEN 
        NN_SF = 4     ! 4-node quadrilateral sub-division
        IY = 1 +    (I_SF-1)/    N_SF
        IX = 1+ MOD (I_SF-1+N_SF,N_SF)
        DO I=1,NN_SF
          LC_SF(I,1) = -1. + 2. * REAL (IX+LC4(1,I)) / REAL (N_SF) 
          LC_SF(I,2) = -1. + 2. * REAL (IY+LC4(2,I)) / REAL (N_SF) 
        ENDDO

C--------------- triangle subdivision ----------------------------------
C --> wow ! what a nice piece of code (even if I do say so myself)!
      ELSEIF (K.eq.3.or.K.eq.6.or.K.eq.10.or.K.eq.15) THEN
        NN_SF = 3
        delta = 1./real(n_sf)
        ix = sqrt(real(i_sf-1)+ .01) +1
        iy = n_sf - (ix*ix-i_sf)/2
        LC_SF(1,1) = 1.- delta * (ix-1)
        LC_SF(1,2) = 1.- delta *  iy
        IF (MOD(I_SF,2).EQ.MOD(IX,2))THEN    ! an 'up' pointing triangle
          LC_SF(2,1) = LC_SF(1,1) - delta
          LC_SF(2,2) = LC_SF(1,2) + delta
          LC_SF(3,1) = LC_SF(1,1) - delta
          LC_SF(3,2) = LC_SF(1,2)
        ELSE
          LC_SF(2,1) = LC_SF(1,1)
          LC_SF(2,2) = LC_SF(1,2) + delta
          LC_SF(3,1) = LC_SF(1,1) - delta
          LC_SF(3,2) = LC_SF(1,2) + delta
        ENDIF
      ELSE
        PRINT*,'*** WARNING: NEN=',NN_FT,' not known in GET_FAC!'
      ENDIF
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE GFACE (NEN,NDIME,ITYPE,IFACE,FACE,NN_F,NN_FT,IOP)
C
C     this returns a list of nodes making up a face of an element
C       ( of NOD nodes in NDIME dimensions, type = ITPYE )
C     in FS(), for face no. IFACE of NN_F nodes on the element boundary,
C       (NN_FT is the total number of nodes, including 'face' nodes). 
C     if ITYPE = 9 then they are just polygons (eg .G format)
C  
C     IOP = 1 --> return just the number of facets per element
C                ( this is only really needed fot facet-stripping )
C     IOP = 2 --> return the actual face info
C     IOP = 3 --> return only the 'corner node' info
C
C                                                  Dan Kidger 21-4-92
C     --> all this 'could' be done automatically using WTHATN and SHAPE !
C         .. so maybe generate the FORTRAN using WTHATN ?
C
C   * need to add NN_E = # of nodes on each edge ..
C

C--- increment the next line whenever a new element is added :-)
      PARAMETER    (MAX_INFO=227+35+10)  !14nb as 5nq's
C     PARAMETER    (MAX_INFO=227+35+10-6)  !14nb as 4nq's
      INTEGER INFO (MAX_INFO)  , FACE(*)
      SAVE NEN_OLD, NDIME_OLD, IFACE_OLD, IBASE
      DATA NEN_OLD, NDIME_OLD, IFACE_OLD /0,0,0/

c----------------------------------------------------------------------
C------>  = NEN,NDIME,#faces,#nodes per face,#mid-face nodes,
C          /nodes (-ve if 'not-a-corner')/ 

c     DATA (INFO (I),I=i,200) / 
      DATA INFO /
     +  2,1, 1, 2,0,  1,2       ! 2n-line
     +, 3,1, 1, 4,0,  1,-2,3,-2 ! 3n-line

     +, 4,2, 1, 4,0, 1,2,3,4                               ! 4nq
     +, 5,2, 1, 5,1, 1,2,3,4,5                             ! 5nq
     +, 8,2, 1, 8,0, 1,-2,3,-4,5,-6,7,-8                   ! 8nq
     +, 9,2, 1, 9,1, 1,-2,3,-4,5,-6,7,-8,9                 ! 9nq
     +,12,2, 1,12,0, 1,-2,-3,4,-5,-6, 7,-8,-9,10,-11,-12   ! 12nq
     +,17,2, 1,17,1, 1,-2,-3,-4,5,-6,-7,-8,9,-10,-11,-12   ! 17nq
     +                                ,13,-14,-15,-16,17 

     +, 3,2, 1, 3,0, 1,2,3                         ! 3n-tri
     +, 6,2, 1, 6,0, 1,-2,3,-4,5,-6                ! 6n-tri
     +,10,2, 1,10,1, 1,-2,-3,4,-5,-6,7,-8,-9, 10   !10n-tri 
     +,15,2, 1,15,3, 1,-2,-3,-4,5,-6,-7,-8,9,-10,-11,-12,13,14,15 !15n-tri

     +, 8,3, 6,4,0, 1,4,8,5,1,5,6,2,1,2,3,4, 7,3,2,6,7,8,4,3,7,6,5,8 ! 8nb

     +,14,3, 6,5,1, 3,4,2,1,-9,  1,2,6,5,-13, 1,5,7,3,-11         ! 14nb
     +,             5,6,8,7,-10, 7,8,4,3,-14, 2,4,8,6,-12 

c     +,14,3, 6,4,0, 3,4,2,1,  1,2,6,5, 1,5,7,3         ! 14nb
c     +,             5,6,8,7, 7,8,4,3, 2,4,8,6

     +,20,3, 6,8,0, 1,-8,7,-12,19,-20,13,-9, 1,-9,13,-14,15,-10,3,-2 
     +,      1,-2,3,-4,5,-6,7,-8, 17,-11,5,-4,3,-10,15,-16         !20nb
     +,     17,-18,19,-12,7,-6,5,-11,  17,-16,15,-14,13,-20,19,-18
     +/
C--------------------------- end_of_list--------------------------------

C---- rescan the database if we haven't currently 'got' this element
      IF (ITYPE.eq.9) THEN    !   'OFF/NFF' format
        if (iop.eq.1) then
          iface = 1           ! ie. return that there is only one face
          return
        elseif (iop.eq.2.or.iop.eq.3) then 
          nn_f  = nen
          nn_ft = nen
          if (iface.gt.0) then
            do i=1,nn_ft
              face(i) = i           ! return a simple list 
            enddo                   ! of consecutive numbers
          else
            do i=1,nn_ft
              face(i) = nn_ft+1-i   ! a reversed list
            enddo
          endif
          return
        endif
      ENDIF
c..................... handle 'Danplot' style .........................

      IF (NEN_OLD.NE.NEN.OR.NDIME_OLD.NE.NDIME) THEN     !- if not 'current'
c * * * * *  Previous was .AND. until 23-6-93  
c -> so it was S.L.O.W and wrong! (eg.if 3nt & 4nq)

        NEN_OLD   = NEN
        NDIME_OLD = NDIME
        IBASE = 1

c.... hmm really I would like to be able to speed this up by storing
c     pointers directly to the required element ? ( or at least put the
C     most 'common' elements first )
c
    2   CONTINUE                 !- loop back point
c       PRINT*,'trying,nen=,ndime=',INFO(IBASE),INFO(IBASE+1)
        IF (NEN.NE.INFO(IBASE).OR.NDIME.NE.INFO(IBASE+1) ) THEN  
          IBASE = IBASE + INFO(IBASE+2) * INFO(IBASE+3) + 5     !* 5 *
          IF (IBASE.LT.MAX_INFO) GOTO 2
        ELSE
          GOTO 1                 !- found this element
        ENDIF
C         ...... unfound so warn user and exit .........
          PRINT*,'*** WARNING, this element not found'//
     +         ' NEN=',NEN,' ndime=',ndime
     +     ,'ITYPE=', itype,' IFACE=', iface
        RETURN
      ENDIF

C-----------------------------------------------------------------------
    1 CONTINUE  ! OK .. extract this facets' info   
      IF (IOP.EQ.1) THEN
        IFACE = INFO(IBASE+2) 
        RETURN
      ENDIF
c.... hmm I 'may' modify FACE() so must always create it !
c.. ? why do I modify FACE() ?
c********** try uncommenting the next line ! ***************
c     IF (IFACE.EQ.IFACE_OLD) RETURN     ! no work at all to do !

c**** Also as we tend to have the same element (eg.20nb) throughout
c**** It must be faster to KEEP a copy of the 6 facets for v.fast access

      IFACE_OLD = ABS (IFACE)
      IBASE_F = IBASE + 4 + (ABS(IFACE)-1) * INFO(IBASE+3)
      NN_FT =         INFO(IBASE+3)
      NN_F  = NN_FT - INFO(IBASE+4)
      IF (IFACE.gt.0) THEN           ! a 'clockwise' face
        DO I=1,NN_FT
          FACE(I) = ABS(INFO(IBASE_F+I))   ! .. just use ABS for now 
        ENDDO
      ELSE                          ! an 'anticlockwise' face (mirrored)
        FACE(1) = ABS(INFO(IBASE_F+1))  !- I would rather not have 
        DO I=2,NN_F                     !- Mirrored facets.. 
          FACE(NN_F+2-I) = ABS(INFO(IBASE_F+I)) 
        ENDDO
        DO I=NN_F+1,NN_FT
          FACE(I) = ABS(INFO(IBASE_F+I))
        ENDDO
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C----------------------------------------------------------------------


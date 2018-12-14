C-----------------------------------------------------------------------       
      SUBROUTINE GET_FILE_TYPE (file_dat,idt)
C
C     This simply gives a file a 'data-type' IDT depending upon its
C     filename extension  
C
      CHARACTER  FILE_DAT*(*)
                                       IDT = 0    ! default = 'unknown'
      IF (INDEX(FiLE_DAT,'.DIS').ne.0) IDT = 1   !- my .PL format(IAN's) 
      IF (INDEX(FiLE_DAT,'.MES').ne.0) IDT = 1   !- my .PL format(IAN's) 
      IF (INDEX(FiLE_DAT,'.MSH').ne.0) IDT = 1   !- my .PL format(MESH5) 
      IF (INDEX(FILE_DAT,'.PL' ).ne.0) IDT = 1   !- my .PL format 
      IF (INDEX(FILE_DAT,'.PL2').ne.0) IDT = 1   !- my .PL + xtra lines
      IF (INDEX(FILE_DAT,'.G'  ).ne.0) IDT = 2   !- 'Objext File Format'
      IF (INDEX(FILE_DAT,'.GEO').ne.0) IDT = 2   !- 'Objext File Format'
      IF (INDEX(FILE_DAT,'.HO' ).ne.0) IDT = 3   !- David Ho's format ?
      IF (INDEX(FILE_DAT,'.NFF').ne.0) IDT = 4   !- NFF format
      IF (IDT.eq.0) THEN
        PRINT*,'file=',FILE_DAT,'is of unknown type !'
      ENDIF
      RETURN
      END

C-----------------------------------------------------------------------       
      SUBROUTINE READ_NFF (U3,NN,NDIM,
     +     GCOORD,IGCRD,INF,NEL,NUMS,INUMS,MEL,PAL)   

c..  don't really need NDIME,ITYPE etc., NLDS too
c
c     this reads in data from the NFF= 'Neutral FILE Format'
c     -- the 'p', and 'pp' read in the polygons
c     'f' the material number and RGB definition
c     -- however NFF also contains the viewing params, colours, etc.
c

      INTEGER*4  INF,MEL
      REAL  GCOORD (IGCRD,INF)
      INTEGER NUMS (INUMS,MEL)       !- the global element steering
     +        ,NUM (27)              !- an element's steering
     +        ,PAL (3,0:255)         !-
     +         ,U3                   !- hmm
      CHARACTER code*2,t*1, line*80

      NDIM = 3

      NN   = 0
      NEL  = 0
      IMAT = 0        !- (was = 2 before)
 1111 CONTINUE
      READ (U3,'(a)',ERR=999) CODE
      IF (CODE.eq.'f ') THEN    ! a material colour
       BACKSPACE (U3)
         IMAT = IMAT + 1
         ICOL = IMAT + 3        ! (Hard coded palette offset !)
       READ (U3,'(A)') LINE
       LINE(1:2) ='  '   !-- zap the 'key-letter'
       READ (LINE,*) R,G,B
       PAL(1,ICOL) = R *256 
       PAL(2,ICOL) = G *256 
       PAL(3,ICOL) = B *256 
       DO J=1,3
         IF (PAL(J,ICOL).EQ.256) PAL(J,ICOL) = 255   !- clip
       ENDDO

C-----------  a polygon / polygon-with-normals ('p','pp') --------------
C .. if 'pp' should really read-in the normals
C .. if 'p'  I can produce the 'flat' normals directly

      ELSEIF (CODE.eq.'p '.or.CODE.eq.'pp') THEN 
      BACKSPACE (U3)
          IF (CODE.eq.'p ') READ (U3,'( a,i2)')   t,NOD
          IF (CODE.eq.'pp') READ (U3,'(2a,i2)') t,t,NOD
        NEL = NEL+1

c       NUMS(NEL,1) = NOD           !(rev. subscripts)
        DO I=1,NOD
          NN = NN + 1
          READ (U3,*) (GCOORD(J,NN),J=1,3)
          NUM (I+1) = NN
        ENDDO
        NDIME = 2
        ITYPE = 9
        CALL PUT_ELEMENT 
     +    (NUMS,INUMS,NEL, NUM, NOD,NDIME,ITYPE, IMAT, 0,0,0)

      ENDIF            !- end_of_keyword_types
      GOTO 1111      !- loop back
  999 CONTINUE
      CALL SET_ALL_DACS ()       !- hmm 
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE READ_NODES (IO,IDT,NN,NDIM,GCOORD,IGCRD,INF,NEL)
c
c  this reads in the nodal co-ordinates for various file formats
c     IDT=1 (my .PL), =2 (OFF) , =3 (D.Ho's)

      INTEGER*4  INF
      REAL GCOORD(IGCRD,INF)
      IF (IDT.eq.1) THEN
C------------------- Dan's .PL format ----------------------------------
        READ (IO,*) NDIM
        READ (IO,*) NN_
        NN = 0
        DO I=1,NN_  ! nicer to have a 'loop-until-error' structure 
          READ(IO,*,ERR=99) II,(GCOORD(J,II),J=1,NDIM)
          NN =  MAX (II,NN)
        ENDDO
   99   CONTINUE
C------------------ 'Object File Format' -------------------------------
      ELSEIF (IDT.eq.2) THEN
        NDIM  = 3       !---  ie. 3D
        READ (IO,*) NN, NEL, NEDGES          ! NEL needs storing !
        READ (IO,*) ((GCOORD(J,I),J=1,NDIM), I=1,NN)
C---------------------- 'David Ho's Format -----------------------------
      ELSEIF (IDT.eq.3) THEN
c       NDIM  = 3
        READ (IO,*)  NDIM
        READ (IO,*)  NN
        READ (IO,*) (II,(GCOORD(J,II),J=1,NDIM),I=1,NN)
      ELSE
        PRINT*,'*** WARNING: data_type=',IDT,' unknown (READ_NODES)'
        NN = 0
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE READ_ELEMS (IO,IDT,NEL,NUMS,INUMS,MEL,NDIM)
c
c     this reads in the 'elements' into NUMS for different file formats
c     IDT=1 (my .PL), =2 (OFF) , =3 (D.Ho's)
c     ** fixed a IMAT bug 9-11-93 !
c     YUK! better to handle each file format seperately
c
      INTEGER*4  MEL
      INTEGER NUMS(INUMS,MEL), NUM(27)
C------------------- Dan's .PL format ----------------------------------
      IF (IDT.eq.1) THEN
        READ (IO,*) NEL_
        NEL = 0
        DO I = 1,NEL_
          READ(IO,*,ERR=99) II,NOD, (NUM(J),J=1,NOD), IMAT
        NDIME = NDIM         ! ie. same as the geometry :-(
        ITYPE = 1                ! default = 1
C                               ... (careful with NDIME)
          CALL PUT_ELEMENT 
     +    (NUMS,INUMS,II, NUM, NOD,NDIME,ITYPE, IMAT,0,0,0)
          NEL = MAX (II,NEL)
        ENDDO
C------------------ 'Object File Format' -------------------------------
      ELSEIF (IDT.eq.2) THEN 
        NDIME = 2       !      with 2D 'polygon' elements
        ITYPE = 9       !      a 'polygonal' element
        DO II=1,NEL
          READ(IO,*) NOD, (NUM(J),J=1,NOD)
          IMAT = NOD       !- set IMAT == # nodes per elem
          ITYPE = 9        !(ie. just a polygon)  (but 3nt & 4nq's are OK)
          CALL PUT_ELEMENT 
     +    (NUMS,INUMS,II, NUM, NOD,NDIME,ITYPE, IMAT,0,0,0)
        ENDDO
C---------------------- 'David Ho's Format -----------------------------
c  *obsolete* surely
c      ELSEIF (IDT.eq.3) THEN
c        ITYPE = 1            ! default = 1
c        READ (IO,*) NEL
c        DO I = 1,NEL
c          NEN = 14
c          READ(IO,*) II, (NUMS(J+1,II),J=1,NEN)
c          NUMS(1,II)     = NEN   ! #nodes per element
c        ENDDO
c        READ(IO,*) (NUMS(2+NEN,II),II=1,NEL)        ! material types
      ELSE
        PRINT*,'*** WARNING: data_type=',IDT,' unknown (READ_ELEMS)'
        NN = 0
      ENDIF
      RETURN
  99  STOP 'Error in reading Elements in DANPLOT format'
      END

C-----------------------------------------------------------------------
      SUBROUTINE READ_LOADS (IO,IDT,NLDS,GDISPS,IGD,IGDISPS,NODOF,NN)
c
c     this reads in the 'displacements' for different file formats
c     IDT=1 (my .PL), =2 (OFF) , =3 (D.Ho's)
c
      INTEGER*4 IGDISPS,IBASE
      REAL GDISPS(IGD,IGDISPS)
      CHARACTER LINE*80

C------------------------- Dan's .PL format ----------------------------
      IF (IDT.eq.1) THEN    ! 
      DO NLDS=1,9999                 !- *NOT* just 1->99 !
        IBASE = NN * (NLDS-1) 
        IF (IBASE+NN.GT.IGDISPS) THEN
          PRINT*,'*** TOO Many Load steps.. only',NLDS-1,' used <CR>'
          read*
          GOTO 123
        ENDIF
        READ(IO,'(A)',end=123) LINE
        WRITE(*,'(a,i3,a,a)') 'load step >',nlds,' <', LINE(1:50)
        DO I=1,NN
          DO J=1,3
            GDISPS(J,IBASE+I) = 123.e-30  ! set the disps to zero
          ENDDO
        ENDDO

        READ (IO,*,ERR=123) NN_
        DO I=1,NN_
c... next line change 14-07-92 !
c         READ (IO,*,ERR=123) II,dummy,(GDISPS(J,IBASE+II),J=1,NODOF)
          READ (IO,*,ERR=123) II      ,(GDISPS(J,IBASE+II),J=1,NODOF)
        ENDDO
      ENDDO
  123 NLDS=NLDS-1
      print*, ' Total No. of Load steps = ',NLDS,'           '
C------------------ 'Object File Format' -------------------------------
c... ie. no-need to ever call this
      ELSEIF (IDT.eq.2.or.IDT.eq.4) THEN
        print*,'** ignoring disps for ''OFF'' format'
        NLDS = 0
C---------------------- 'David Ho's Format -----------------------------
c   *obsolete* surely
c      ELSEIF (IDT.eq.3) THEN
c        NODOF = 3      !--- always in 3D
c        DO NLDS=1,99
c        IBASE = NN * (NLDS-1) 
c        IF (IBASE+NN.GT.IGDISPS)THEN
c          PRINT*,'*** TOO Many Load steps.. only',NLDS-1,' used'
c          GOTO 99
c        ENDIF
c        READ(IO,'(A)',ERR=99) LINE
c        WRITE(*,'(a,i3,a,a)') 'load step >',nlds,' <', LINE(1:50)
c        DO I=1,NN
c          DO J=1,3
c            GDISPS(J,IBASE+I) = 123.e-30  ! set the disps to zero
c          ENDDO
c        ENDDO
c        DO I=1,NN
c          READ (IO,*) II,(GDISPS(J,IBASE+II),J=1,NODOF)
c        ENDDO
c      ENDDO
   99 NLDS=NLDS-1
      print*, ' Total No. of Load steps = ',NLDS,'           '
      LD = NLDS
      ELSE
        PRINT*,'*** WARNING: data_type=',IDT,' unknown (READ_ELEMS)'
        NLDS = 0
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE GET_MESH_RANGE (GCOORD,IGCRD, NN,NDIM, COA,DIAG,diag2)
C
C     this calculates range of the data in GCOORD
C      also the diagonal of the bounding box DIAG
C     (NN+1)  x,y,z minima,    
C     (NN+2)  x,y,z lengths,   
C     (NN+3)  x,y,z (centroid)
C     COA = 'middle' of the range of co-ordinats in GCOORD
C
C     3-11-93 DIAG2 = the Max. of the x,y and z widths

      REAL GCOORD (IGCRD,*), COA(*)

C---------------- null the z-coords if necesary ------------------
      DO J=NDIM+1,3
        DO I=1,NN
          GCOORD(J,I) = 0.
        ENDDO
      ENDDO
C------------  set maxima/minima to first node --------
      DO J=1,3    
        DO I= NN+1,NN+3
          GCOORD(J,I) = GCOORD(J,1)
        ENDDO
      ENDDO
      DIAG=0.
C------------- loop and find the MAX,MIN & CENTROID --------------
      DO I=1,NN
        DO J=1,3
          GCOORD(J,NN+1) = MIN(GCOORD(J,NN+1), GCOORD(J,I) )
          GCOORD(J,NN+2) = MAX(GCOORD(J,NN+2), GCOORD(J,I) )
          GCOORD(J,NN+3) =    (GCOORD(J,NN+3) +GCOORD(J,I) /REAL(NN))
        ENDDO
      ENDDO
C------- find the 'chararistic model length' = major daigonal ---------
      DIAG = SQRT (( GCOORD(1,NN+2)-GCOORD(1,NN+1) )**2
     +        +    ( GCOORD(2,NN+2)-GCOORD(2,NN+1) )**2
     +        +    ( GCOORD(3,NN+2)-GCOORD(3,NN+1) )**2 )

      DIAG2 = MAX ( GCOORD(1,NN+2)-GCOORD(1,NN+1) ,
     +              GCOORD(2,NN+2)-GCOORD(2,NN+1) ,
     +              GCOORD(3,NN+2)-GCOORD(3,NN+1)   )

C--------- set the CENTRE_OF_ATTENTION to the 'middle' of the object ---
      DO I=1,3
        COA(I)= (GCOORD(I,NN+1) + GCOORD(I,NN+2)) /2.
      ENDDO
      RETURN
      END
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


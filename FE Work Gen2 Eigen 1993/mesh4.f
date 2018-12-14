      options (dreal, intl, fullcheck, undef)

      PROGRAM DANS_MESH_GENERATOR

      PARAMETER ( MELS =8000 )   !-- Max # of elements (change as needed)
c------------------------------------------------------------------------c
c           D A N   K I D G E R 'S   M E S H    G E N E R A T I O N      c
c------------------------------------------------------------------------c
      CHARACTER  VERSION*15
      PARAMETER (VERSION =             'Version 1.0'                   )
c------------------------------------------------------------------------c
c   ** DESCRIPTION **                                                    c
c                                                                        c
c   This programs is designed to create 2d/3d finite element meshes      c
c   from a 'Rule-based' data file using *Keywords                        c
c------------------------------------------------------------------------c
c   ** REVISION HISTORY **
c
c
c  8-01-93 Added co-incident node deletion
c 24-01-93 built-up .. first working copy ?     
c 19-01-93 MESH1 .. first version (abstracted from GEN4.FOR)
c-----------------------------------------------------------------------

      PARAMETER ( MNN = MELS     ! = 2 for 2d (=3 for 2d/3d)
     +           ,INF = 3        ! = 2 for 2d (=3 for 2d/3d)
     +           ,IGC = INF      !
     +           ,INUMS=28   )   ! = max nodes per element

c----------- problem sized arrays --------------------------------------
      REAL  GC  (IGC,MELS)       !- co-ordiantes of all the nodes
      INTEGER NUMS (INUMS,MNN)   !- the elements' nodes etc.
     +          ,P (MNN)         !- a 'workspace' of pointers

c----------------- 'small' arrays --------------------------------------
c     INTEGER  IOPS (30)      !  program control options
 
c--------------- character strings -------------------------------------
      CHARACTER KEYWORD*40   ! data input keyword 'token'
     +            ,FILE*100  ! input data file name
     +          ,CMNAM@*100  ! command line function
     +              ,COL*7   !-- ANSI colours
     +              ,COL1*7   !-- ANSI colours
     +              ,COL2*7   !-- ANSI colours
c    +          ,FDATE@*40   ! Date function
c    +           ,TIME@*8    ! Time function
c    +      ,TITLES(10)*80   ! The main titles
c    +            ,LINE*255  ! 'generic' line of text
      INTEGER U0             ! .INP : the original data file
     + ,U1                   !  --- : the current data file being read
c    + ,U6                   ! .DP2 : The Output 'keyword' mesh file
     + ,U7                   ! .PL  : The Danplot output file
      EXTERNAL COL
c------------------ DATA stataments ------------------------------------
      DATA    U0/40/,  U7/10/
      DATA NEL,NN/2*0/            !-- initially NO elements etc.
      DATA NDIM,NODOF/3,3/        !-- default is 3D

c      DATA (IOPS(I),I=1,10)/     !-------- default options ---------
c     + 0,   !1!    = 1  to inhibit any output for this load-step

C =====================================================================
      WRITE (*,'(A)') '-------- Dan''s Mesh Generator -------'
     +               ,'Version: '//VERSION

      FILE = CMNAM@()    !-- the command-line
      IF (FILE.EQ.' ') CALL SELECT_FILE@('*.d*',FILE,*999)
      CALL UPCASE (FILE)

      OPEN (U0,FILE=FILE,STATUS='OLD') 
C----------------- Parse the data file (s) -----------------------------
c.... get the next keyword: may skip out  at EOF or return keyword
      U1 = U0                         !-- start from the base file = U0
      DO IKEYWORD = 1, 999999
        CALL GET_KEYWORD (U1,U0,KEYWORD)
        col1 = col(33)
        col2 = col(37)
c       WRITE(*,*) KEYWORD//COL(36)
        WRITE(*,'(I3,A,A,A,A,A)') IKEYWORD, ' ',COL1,KEYWORD,COL2
c       WRITE(*,'(I3,A,A,A,A,A)') IKEYWORD, ' ',COL(33)//KEYWORD//COL2

      IF (KEYWORD.EQ.'*EOF') THEN
        GOTO 888                          !-- exit the program nicely
      ELSEIF (KEYWORD.EQ.'*TWO_DIMENSIONAL') THEN
        NODOF = 2
      ELSEIF (KEYWORD.EQ.'*THREE_DIMENSIONAL') THEN
        NODOF = 3

      ELSEIF (KEYWORD.EQ.'*IMPORT_OFF') THEN
        CALL R_OFF (U1,NODOF,NN,GC,IGC,NEL,NUMS,INUMS)
      ELSEIF (KEYWORD.EQ.'*IMPORT_PL') THEN
        CALL R_PL (U1,NODOF,NN,GC,IGC,NDIM,NEL,NUMS,INUMS)
      ELSEIF (KEYWORD.EQ.'*BASIC_SLOPE_MESH') THEN
        CALL R_READGM (U1,NODOF,NN,GC,IGC,NEL,NUMS,INUMS)

      ELSEIF (KEYWORD.EQ.'*MOVE_NODES') THEN
        CALL R_NODMOV(U1,GC,IGC,NDIM,NN)

      ELSEIF (KEYWORD.EQ.'*ELEMENT_TYPES') THEN
        CALL R_GETELT(U1,NUMS,INUMS,NEL,NODOF)

      ELSEIF (KEYWORD.EQ.'*DELETE_MATERIALS') THEN
        CALL DEL_MATS    (U1,GC,IGC,NN,NUMS,INUMS,NEL,P,NDIM)
      ELSEIF (KEYWORD.EQ.'*DELETE_ORPHAN_NODES') THEN
        CALL DEL_O_NODES (U1,GC,IGC,NN,NUMS,INUMS,NEL,P,NDIM)
      ELSEIF (KEYWORD.EQ.'*DELETE_COINCIDENT_NODES') THEN
        CALL DEL_COIN_NODES (U1,GC,IGC,NN,NUMS,INUMS,NEL,P,NDIM)

      ELSEIF (KEYWORD.EQ.'*4NQS_FROM_4_3NTS') THEN   !-- 'teapot.g'
      CALL MUNG_4NQ  (U1,GC,IGC,NUMS,INUMS,P,NDIM,NN,NEL)

      ELSEIF (KEYWORD.EQ.'*MIRROR_EVEN_ELEMENTS') THEN
      CALL MIRROR_EVENS (U1,GC,IGC,NUMS,INUMS,P,NDIM,NN,NEL)

      ELSEIF (KEYWORD.EQ.'*SHIFT') THEN
        CALL SHIFT (U1,GC,IGC,NN,NDIM)
      ELSEIF (KEYWORD.EQ.'*SCALE') THEN
        CALL SCALE (U1,GC,IGC,NN,NDIM)
      ELSEIF (KEYWORD.EQ.'*X_TO_Y_TO_Z') THEN
        CALL X2Y2Z  (GC,IGC,NN,NDIM)
      ELSEIF (KEYWORD.EQ.'*X_TO_Y') THEN
        CALL X2Y    (GC,IGC,NN,NDIM)
      ELSEIF (KEYWORD.EQ.'*WRAP_AROUND_Y') THEN
        CALL WRAP_Y (GC,IGC,NN,NDIM)
      ELSEIF (KEYWORD.EQ.'*WRAP_CIRCLE_Y') THEN
        CALL WRAP_CIRCLE_Y (GC,IGC,NN,NDIM)
      ELSEIF (KEYWORD.EQ.'*WRAP_SPHERE') THEN
        CALL WRAP_SPHERE (GC,IGC,NN,NDIM)
      ELSEIF (KEYWORD.EQ.'*WRAP_TWIST_Y') THEN
        CALL WRAP_TWIST_Y (U1,GC,IGC,NN,NDIM)

      ELSEIF (KEYWORD.EQ.'*WRITE_DANPLOT') THEN
        OPEN(U7,FILE=FILE(1:INDEX(FILE,'.D'))//'PL')
        CALL WR_DANPLOT(U7,FILE,GC,IGC,NUMS,INUMS,NDIM,NN,NEL)

      ELSEIF (KEYWORD.EQ.'*WRITE_OFF') THEN
        OPEN(U7,FILE=FILE(1:INDEX(FILE,'.D'))//'GEO')
        CALL WR_OFF (U7,FILE,GC,IGC,NUMS,INUMS,NDIM,NN,NEL)

C----------------------------------------------------------------------
        ELSE
          CALL MYERROR(1,'unknown keyword:'//KEYWORD)
        ENDIF      !--- end of possible keywords

      WRITE(*,'(A,I5,A,I5,A)') 'Now',NN,' Nodes and',NEL,' Elements'
      ENDDO
C------------------------- end of program ----------------------------
  888 STOP '> > >       Program completed   :-)       < < <' 
  999 STOP  '***** CRASHED OUT ******'
      END



C-----------------------------------------------------------------------
      FUNCTION COL(icode)
C
C     function to output ANSI colour control strings
C
      CHARACTER COL*7
      WRITE (COL,'(A,I2.2,A)')    CHAR(27)//'[1;',ICODE,'m'
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GET_KEYWORD (U1,U0,KEYWORD)
C
C     This gets the next keyword from the data file
C     '*EOF' is returned when all data is exhausted
C
      CHARACTER KEYWORD*(*), LINE*80
      INTEGER U1,U0
    1 READ (U1,'(A)',IOSTAT=IOS) LINE
      CALL IN_TEST(U1,IOS,*1,*999)  
      IF (U1.LT.U0) GOTO 999
      IF (LINE(1:1).NE.'*') GOTO 1
      KEYWORD  = LINE(1:INDEX(LINE,' '))   !-- strip off any comments
      CALL UPCASE(KEYWORD)                 !-- Uppercase is easier to handle
      RETURN
 999  KEYWORD ='*EOF'         !- end of the all the data
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE IN_TEST(IO,IOSTAT,*,*)
C
C     this tests the nature of a data read error
C     return to *1 if 'soft', and *2 if fatal
C
      CHARACTER LINE*80
      LOGICAL FILE_OPEN

      IF (IOSTAT.EQ.0) RETURN          !---- No error  !

      BACKSPACE (IO,ERR=99)
      READ (IO,'(A)',END=99) LINE

      IF (LINE(1:1).EQ.'='      !--- allowable comment lines 
     +  .OR.  LINE(1:1).EQ.'#'
     +  .OR.  LINE(1:1).EQ.'c'
     +  .OR.  LINE(1:1).EQ.'C'
     +  .OR.  LINE(1:1).EQ.' '
     +  .OR.  LINE(1:1).EQ.'/') THEN
        RETURN 1
      ELSEIF (LINE(1:1).EQ.'!') THEN !----An 'echoed' comment line -----
        WRITE (*,'(A)') LINE
        RETURN 1
      ELSEIF (LINE(1:1).EQ.'%') THEN !----An 'echoed' comment line -----
        WRITE (90,'(A)') LINE        !       to the '.PLT' file
        RETURN 1
      ELSEIF (LINE(1:1).eq.'*') THEN !------- A new Keyword -----------
        BACKSPACE (IO)  ! point at the '*'
        RETURN 2
      ELSEIF (LINE(1:1).eq.'<') THEN !------- read_from_another_file ---
        IO = IO + 1
        OPEN (IO,FILE = LINE(2:INDEX(LINE,'  ')),IOSTAT=IOSTAT)
        IF (IOSTAT.NE.0) THEN ! error = 'new data file not found'
          CALL MYERROR (2,'missing data file-'//LINE)
          RETURN 2  !-- crash out 
        ENDIF
        RETURN 1   !--- this line was missing ! 14-9-92
c- - - - - - - - - - - - - - - - - - - 
      ELSE                    !--------- 'rubbish' found in file -------
        CALL MYERROR(2,'error in data file format !')
        RETURN 2
      ENDIF
c- - - - - - - - - - - - - - - - - - - 
C... should always jump past this line ...
c- - - - - - - - - - - - - - - - - - - 
   99 CONTINUE  
      IO = IO - 1           !-- revert to the last file
      INQUIRE (IO,EXIST=FILE_OPEN)
      IF (FILE_OPEN) THEN   !------------ 'end-of-file'----------------
        IF (IO.ge.40) RETURN 1           !-- was 2 !! 1-2-93
        RETURN 2
      ELSE                  !---------- 'file was not open'------------
        RETURN 2
      ENDIF
      END
C-----------------------------------------------------------------------
      SUBROUTINE MYERROR (ILEVEL,STRING)
C
C     This outputs a warning/error message to stdout 
C
C.... maybe extend with further options..
C
C------------ exit the program via the FTN77 error hander --------
      CHARACTER STRING*(*)
      IF (ILEVEL.EQ.1) THEN
        WRITE (*,'(A,2X,A)') '<WARNING> :',STRING
      ELSEIF (ILEVEL.EQ.2) THEN
        WRITE (*,'(A,2X,A)') '**ERROR** :',STRING
        CALL ERROR@ (STRING)
      ELSEIF (ILEVEL.EQ.3) THEN
        WRITE (*,'(A,2X,A)') '** PROGRAM BUG  ** :',STRING
        CALL ERROR@ (STRING)
      ELSE
        WRITE (*,'(A,I2)') 'hmm.. unknown error level =',ILEVEL
      ENDIF
      END
c-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE R_OFF (IO,NODOF,NN,GC,IGC,NEL,NUMS,INUMS)
C
C     This imports a mesh in 'Object File Format'
C     ( These are 2d polygons (mostly 3 node triangles) in a 3d mesh
C     .. no need to test for comments as they never have any!
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*)

      NN_OLD = NN
      NEL_OLD = NEL

    1 READ (IO,*,IOSTAT=IOS)  NN,NEL,NEDGES
      CALL IN_TEST(IO,IOS,*1,*999)
      READ (IO,*)  ((GC(J,I+NN_OLD),J=1,3),I=1,NN)
      WRITE(*,'(A,I4,A)')  '>> ',NN,' nodes read'
      DO IEL=1,NEL
        READ (IO,*) NEN,(NUMS(J,IEL+NEL_OLD),J=1,NEN)
        DO J=1,NEN
          NUMS(J,IEL) = NUMS(J,IEL) + NN_OLD   !-- offsets
        ENDDO
        DO J=NEN+1,INUMS
          NUMS(J,IEL) = 0    !- blank-out
        ENDDO
        NUMS (INUMS  ,IEL) = NEN      !- NEN   (mostly = 3)
        NUMS (INUMS-1,IEL) = 2        !- NDIME = 2d
        NUMS (INUMS-2,IEL) = 9        !- ITYPE = 9
        NUMS (INUMS-3,IEL) = NEN      !- mat #  (use NEN)
        NUMS (INUMS-4,IEL) = IEL      !- IP -x
        NUMS (INUMS-5,IEL) = IEL      !- IQ -y
        NUMS (INUMS-6,IEL) = IEL      !- IS -z
      ENDDO
      IO = IO - 1     !-- exit from this file  :-)
      WRITE(*,'(A,I4,A)')  '>> ',NEL,' elements read'
  999 CONTINUE
      NN = NN + NN_OLD
      NEL= NEL+NEL_OLD
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE R_PL (IO,NODOF,NN,GC,IGC,NDIM,NEL,NUMS,INUMS)
C
C     This imports a mesh in 'Danplot format'
C     .. again there are never any comments to handle
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*)

      NN_OLD = NN
      NEL_OLD = NEL

   1  READ (IO,*,IOSTAT=IOS) NDIM,NN
      CALL IN_TEST(IO,IOS,*1,*999)

      READ (IO,*) (ID,(GC(J,I+NN_OLD),J=1,NDIM),I=1,NN)
      WRITE(*,'(A,I4,A)')  '>> ',NN,' nodes read'

      READ (IO,*) NEL
      DO IEL=1,NEL
        IEL_N = IEL + NEL_OLD
        READ (IO,*) ID,NEN,(NUMS(J,IEL_N),J=1,NEN),IMAT
        DO J=1,NEN
          NUMS(J,IEL) = NUMS(J,IEL_N) + NN_OLD   !-- offsets
        ENDDO
        DO J=NEN+1,INUMS
          NUMS(J,IEL_N) = 0    !- blank-out
        ENDDO
        NUMS (INUMS  ,IEL_N) = NEN      !- NEN   (eg. 8)
        NUMS (INUMS-1,IEL_N) = NDIM     !- NDIM  (eg. 2)
        NUMS (INUMS-2,IEL_N) = 1        !- ITYPE = 1
        NUMS (INUMS-3,IEL_N) = IMAT     !- mat #  (use NEN)
        NUMS (INUMS-4,IEL_N) = IEL      !- IP -x
        NUMS (INUMS-5,IEL_N) = IEL      !- IQ -y
        NUMS (INUMS-6,IEL_N) = IEL      !- IS -z
      ENDDO
      CLOSE (IO)
      WRITE(*,'(A,I4,A)')  '>> ',NEL,' elements read'
  999 CONTINUE
      NN = NN + NN_OLD
      NEL= NEL+NEL_OLD
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE R_READGM(IO,NODOF,NN,GC,IGC,NEL,NUMS,INUMS)
C
C     this routine reads the element geometry puts it into 'long' form
C     19-1-93 moved from L1A.FOR to MESH1.FOR for mesh generation
C
      REAL GC(IGC,*), COORD(20,3)
      REAL    TOP(90), BOT(90), DEPTH(90), BRDTH(90)
      INTEGER NUM(20), NUMS(INUMS,*)

c...... NODOF for the basic mesh should come from here ! .......

      NZE = 1
   1  READ (IO,*,IOSTAT=IOS) NDIME, NOD, ITYPE
        CALL IN_TEST(IO,IOS,*1,*999)
   2  READ (IO,*,IOSTAT=IOS) NXE,NXS 
        CALL IN_TEST(IO,IOS,*2,*999)
   3  READ (IO,*,IOSTAT=IOS) (TOP(I),I=1,NXS+1)
        CALL IN_TEST(IO,IOS,*3,*999)
   4  READ (IO,*,IOSTAT=IOS) (BOT(I),I=1,NXE+1)
        CALL IN_TEST(IO,IOS,*4,*999)
   5  READ (IO,*,IOSTAT=IOS) NYE,NYS 
        CALL IN_TEST(IO,IOS,*5,*999)
   6  READ (IO,*,IOSTAT=IOS) (DEPTH(I),I=1,NYE+1)
        CALL IN_TEST(IO,IOS,*6,*999)

      IF (NDIME.EQ.3) THEN
   7  READ (IO,*,IOSTAT=IOS) NZE
        CALL IN_TEST(IO,IOS,*7,*999)
   8  READ (IO,*,IOSTAT=IOS) (BRDTH(I),I=1,NZE+1)
        CALL IN_TEST(IO,IOS,*8,*999)
      ENDIF
C---------------- form the global co-ord array GCOORD ------------
      NDIM = NDIME    !--- ie 3D elements in '3D'
      NN  = 0 
      IEL = 0
      DO IP=1,NXE
        DO IQ=1+NYS*MIN(MAX(IP-NXS,0),1),NYE
          DO IS=1,NZE
            IEL = IEL + 1
            DO J=NOD+1,INUMS
              NUMS(J,IEL) = 0   !--- clear the other columns
            ENDDO
            NUMS (INUMS  ,IEL) = NOD       !-- nodes per elem
            NUMS (INUMS-1,IEL) = NDIME     !-- NODOF of element
            NUMS (INUMS-2,IEL) = 1         !-- element type
            NUMS (INUMS-3,IEL) = 1         !-- material number ??
            NUMS (INUMS-4,IEL) = IP        !- maybe these belong
            NUMS (INUMS-5,IEL) = IQ        !- in their own table 
            NUMS (INUMS-6,IEL) = IS        !- for ease of access
            IF (NODOF.EQ.2)
     +      CALL SLOGE2(IP,IQ,NXS,NYS,NYE,TOP,BOT,DEPTH,COORD,20,NUM)
            IF (NODOF.EQ.3)                !-- 8nq only
     +      CALL SLOGET(IP,IQ,IS,NXS,NYS,NYE,NZE,TOP,BOT,DEPTH,BRDTH,
     +      COORD,20,NUM,NOD)      !-- 8nb/14nb/20nb
            DO I=1,NOD
              NUMS(I,IEL) = NUM(I)
              NN = MAX (NUM(I),NN)
              DO J=1,NDIM
                GC (J,NUM(I)) = COORD(I,J)
              ENDDO
            ENDDO
          ENDDO         ! loop       -z
        ENDDO          ! the        -y
      ENDDO           ! elements   -x
      NEL = IEL
      WRITE (*,'(A,I5,A,I5)') '>> NN=',NN, ' NEL=',NEL
      RETURN
  999 PRINT*,' *** Warning: mesh reading failed'
      RETURN
      END
C ------------------------------------------------------------------
      SUBROUTINE R_GETELT (IO,NUMS,INUMS,NEL,NDIM) 
C                
C     this reads in the element material types in x,y,z from/to
C     format - storing the material property type in NUMS after
C     the number of nodes and node numbers 
C     -- data in the format /mat type/,/xfrom,xto/, etc
C     19-1-93  .. moved to MESH1.FOR from L1A.FOR
C
      INTEGER NUMS(INUMS,*), XF(5),XT(5)

      DO KK=1,9999
    1   READ(IO,*,IOSTAT=IOS) MATYPE, (XF(I),XT(I),I=1,NDIM)
        CALL IN_TEST(IO,IOS,*1,*999)
        IC = 0
        DO IEL=1,NEL
          NOD = NUMS (INUMS,IEL)  
          IFLAG = 1
          DO J=1,NDIM
            IPQR = NUMS(INUMS-3-J,IEL)
            IF (IPQR.LT.XF(J).OR.IPQR.GT.XT(J) ) IFLAG = 0   !-Not OK
          ENDDO
          IF (IFLAG.EQ.1) THEN
            IC=IC+1
            NUMS (INUMS-3,IEL) = MATYPE    
          ENDIF
        ENDDO   !---- loop elements
        WRITE(*,'(A,I4,A,I4,A,I4)')
     +  '>> ',KK,' :',IC,' Elements of type',MATYPE
      ENDDO
  999 RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE MUNG_4NQ (IO,GC,IGC,NUMS,INUMS,P,NDIM,NN,NEL)
C
C     This glues 4 consecutive 3nt's into a 4nq (eg. TEAPOT.G)
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*), P(*)
      IC = 0
      DO IEL=1,NEL,4
        IC = IC +1
        NUMS (1,IC) = NUMS(1,IEL  )
        NUMS (2,IC) = NUMS(2,IEL  )
        NUMS (3,IC) = NUMS(1,IEL+2)
        NUMS (4,IC) = NUMS(2,IEL+2)
        NUMS (INUMS,IC) = 4
      ENDDO
      NEL = IC
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE MIRROR_EVENS (IO,GC,IGC,NUMS,INUMS,P,NDIM,NN,NEL)
C
C     This reverses the order of the nodes of alternate elements (eg HEAD.G)
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*), P(*)
      DO IEL=1,NEL,2
        NEN = NUMS(INUMS,IEL)
        DO I1 = 1,NEN/2
          I2 = NEN+1 - I1
          NUM  = NUMS(I1,IEL)
          NUMS (I1,IEL) = NUMS(I2,IEL)
          NUMS (I2,IEL) = NUM
        ENDDO
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE R_NODMOV(IO,GC,IGC,NDIM,NN)
C
C     This moves nodes based on their 'old' and 'new' coords
C
      PARAMETER (MAX_NODES = 100)
      REAL GC(IGC,*), COLD(3), CNEW(3)
      INTEGER NUMLST(MAX_NODES) 

      WILD = 230964.
    1 DO K=1 , 999999
        READ(IO,*,IOSTAT=IOS) (COLD(I),I=1,NDIM) ,(CNEW(I),I=1,NDIM)
        CALL IN_TEST(IO,IOS,*1,*999)
        CALL FIND_NODES (NN,GC,IGC,COLD,NDIM,NUMLST,MAX_NODES,NNOD)
        DO I=1,NNOD
          DO J=1,NDIM
            IF (ABS(COLD(J)-WILD)/WILD.GT.1.E-4) 
     +      GC (J,NUMLST(I)) = CNEW(J)
          ENDDO
        ENDDO
        WRITE(*,'(3(A,I3))')   '>> ',K,':',NNOD,' Nodes moved'
      ENDDO
  999 RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE DEL_MATS (IO,GC,IGC,NN,NUMS,INUMS,NEL,P,NDIM)
C
C     this will delete any elements of a given material type
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*), P(*)

      DO KK = 1,999
    1   READ (IO,*,IOSTAT=IOS) ITYPE
        CALL IN_TEST (IO,IOS,*1,*999)
        NEL1 = NEL    
        IC = 0
        DO IEL=1,NEL
          NOD  = NUMS(INUMS  ,IEL)
          IMAT = NUMS(INUMS-3,IEL)
          IF (IMAT.NE.ITYPE) THEN               !-- If a 'good' element ..
            IC = IC + 1 
            DO J = 1,INUMS
              NUMS(J,IC)=NUMS(J,IEL)            !-- store it at IC
            ENDDO
          ENDIF
        ENDDO
        NEL = IC
        WRITE(*,'(A,I4,A,I4,A)')  '>> ',KK,':',
     +               NEL1-NEL,' elements deleted'
      ENDDO
  999 RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE DEL_O_NODES (IO,GC,IGC,NN,NUMS,INUMS,NEL,P,NDIM)
C
C     this will delete any 'unreferenced' nodes
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*), P(*)

      NN1  = NN           !-- original mesh sizes  (for print-out)

      DO I=1,NN
        P(I) = 0          !-- null the pointer list
      ENDDO

      DO IEL=1,NEL  !-------- light-up ALL valid nodes -----------------
        NOD  = NUMS(INUMS  ,IEL)
        IMAT = NUMS(INUMS-3,IEL)
        DO J = 1,NOD           
          IF (IMAT.NE.0) P(NUMS(J,IEL)) = 1
        ENDDO
      ENDDO
      CALL RESEQ_P (P,NN,NN_NEW)  

      CALL UPDATE_GC (GC,IGC,NDIM,NN, P)
      CALL UPDATE_NUMS (NUMS,INUMS,NEL, P)

      WRITE(*,'(A,I4,A)')  '>> ',NN-NN_NEW,' nodes deleted'
      NN = NN_NEW
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE DEL_COIN_NODES (IO,GC,IGC,NN,NUMS,INUMS,NEL,P,NDIM)
C
C     This will 'condense' all nodes that share the same position
C     using FIND_NODE to get the node number of a node
C
      REAL GC(IGC,*)
      INTEGER P(*), NUMS(INUMS,*), PERCENT

      DO I=1,NN
        P(I) = 0     !--- initialy no known nodes
      ENDDO

C---------------- index the node to their 'first occurence' ------------
      DO I=NN,1,-1    !-- only reversed for 'style' :-) 
c        PRINT*,I,'/',NN,char(27)//'[A'
        PERCENT = 100.*REAL(I)/REAL(NN)
        IF (MOD(PERCENT,5).EQ.0) THEN
          WRITE(*,'(I4,A)')  PERCENT,' %'//CHAR(27)//'[A' 
        ENDIF
        CALL FIND_NODE (GC,IGC,NDIM,GC(1,I),1,NN,INODE)
        P(I) = INODE      !-- the 'lowest' Node # at this point
      ENDDO

      CALL UPDATE_NUMS  (NUMS,INUMS,NEL, P)

      CALL  DEL_O_NODES (IO,GC,IGC,NN,NUMS,INUMS,NEL,P,NDIM)
C     WRITE(*,'(A,I4,A)')  '>> ',NN-NN_NEW,' nodes condensed'
c     NN = NN_NEW
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FIND_NODES (NN,GC,IGC,POINT,NDIM,NUMLIST,INUMLIST,IC)
C
C     given a co-ordinate, this returns its node number(s) in NUMLIST
C
      REAL GC(IGC,*),POINT(*)
      INTEGER NUMLIST(INUMLIST)

      IC = 0
      IFROM = 1
      DO I=1,INUMLIST
        CALL FIND_NODE (GC,IGC,NDIM,POINT,IFROM,NN,INODE)
        IF (INODE.GT.0) THEN         !-- node found
          IC = IC + 1
          NUMLIST(IC) = INODE
          IFROM = INODE + 1                   
        ELSE
          RETURN                     !-- No more nodes to find
        ENDIF
      ENDDO
      CALL MYERROR (1,'Too many nodes found in FIND_NODES')
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FIND_NODE (GC,IGC,NDIM,XY,IFROM,ITO,INODE)
C
C     This searches the node-list GC (IFROM to ITO) to find the node
C     number for a given coordinate (INODE)

      REAL GC(igc,*), XY(*)

      WILD = 230964.
      TOL = 1.E-4
      INODE = 0
      DO I=IFROM,ITO
        DO J=1,NDIM
          IF (ABS(XY(J)-  WILD ).GT.TOL .AND.
     +        ABS(XY(J)-GC(J,I)).GT.TOL) GOTO 1    !-- give-up !
        ENDDO
        INODE = I
        RETURN        !-- OK we have found the right node
    1   CONTINUE   !-- carry on looking
      ENDDO

      INODE = 0       !-- The node was NOT found
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE RESEQ_P (P,NN, IC)
C
C     this resequences the non-zero pointers in P()
C
      INTEGER P(*), NN, IC
      IC = 0
      DO I=1,NN
        IF (P(I).GT.0) THEN
          IC = IC + 1
          P(I) = IC
        ENDIF
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE UPDATE_GC (GC,IGC,NDIM,NN, P)
C
C     This re-orders the nodal co-ordinates in GC based on the NEW order
C     held in P
C
      REAL GC(IGC,*)
      INTEGER P(*)
      IC = 0
      DO INODE=1,NN
        IF (P(INODE).GT.0) THEN
          IC = IC + 1
          DO J=1,NDIM
            GC(J,P(INODE)) = GC(J,INODE)
          ENDDO
        ENDIF
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE UPDATE_NUMS (NUMS,INUMS,NEL, P)
C
C     This re-orders the node-steering in NUMS based on the NEW order
C     held in P
C
      INTEGER NUMS(INUMS,*), P(*)
      DO IEL=1,NEL
        NOD = NUMS(INUMS,IEL)
        DO J=1,NOD
          INODE = NUMS(J,IEL)
          NUMS(J,IEL) = P(INODE)   !-- if ZERO delete elem ?
        ENDDO
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SCALE (IO,GC,IGC,NN,NDIM)
C
C     Scale all the coordinates  in  x,y,(z)
C
      REAL GC(IGC,*) , XS(3)

    1 READ (IO,*,IOSTAT=IOS) (XS(J),J=1,NDIM)
      CALL IN_TEST (IO,IOS,*1,*999)

      DO I=1,NN
        DO J=1,NDIM
          GC(J,I) = GC(J,I) * XS(J)
        ENDDO
      ENDDO
  999 RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SHIFT (IO,GC,IGC,NN,NDIM)
C
C     Shift all the coordinates by a read-in offset  (hi-lighting?)
C
      REAL GC(IGC,*) ,XS(3)

    1 READ (IO,*,IOSTAT=IOS) (XS(J),J=1,NDIM)
      CALL IN_TEST (IO,IOS,*1,*999)

      DO I=1,NN
        DO J=1,NDIM
          GC(J,I) = GC(J,I) + XS(J)
        ENDDO
      ENDDO
  999 RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE X2Y2Z (GC,IGC,NN,NDIM)
C
C     This simply mungs the x coordinate to y, y to z and z to x
C
      REAL GC(IGC,*)
      DO I=1,NN
        X = GC(1,I)
        Y = GC(2,I)
        Z = GC(3,I)
        GC(1,I) = Z
        GC(2,I) = X
        GC(3,I) = Y
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE X2Y (GC,IGC,NN,NDIM)
C
C     This simply mungs the x coordinate to y, y to -x
C
      REAL GC(IGC,*)
      DO I=1,NN
        X = GC(1,I)
        Y = GC(2,I)
        GC(1,I) =  Y
        GC(2,I) = -X
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE WRAP_Y (GC,IGC,NN,NDIM)
C
C     This MUNGs the co-ordinates from a 'cylindrical' to 'cartesian'
C
      REAL GC(IGC,*)

      DO I=1,NN
        RAD = -GC(2,I)               !--- Radii are all -ve :-)
        ANG =  GC(3,I) * 3.14159265/180.
        GC(2,I) = RAD * SIN(ANG)     !---- apply cylindical transf.
        GC(3,I) = RAD * COS(ANG)
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE WRAP_TWIST_Y (IO,GC,IGC,NN,NDIM)
C
C     This Twists the x-z coordinates up the y-axis
C
      REAL GC(IGC,*)

    1 READ (IO,*,IOSTAT=IOS) TWIST_LEN
      CALL IN_TEST (IO,IOS,*1,*999)
      DO I=1,NN
        ANG= GC(2,I) * TWIST_LEN  * 3.14159265/180.
        XO = GC(1,I)
        ZO = GC(3,I)
        GC(1,I) =   XO * COS(ANG) + ZO * SIN(ANG)
        GC(3,I) = - XO * SIN(ANG) + ZO * COS(ANG)
      ENDDO
  999 RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE WRAP_CIRCLE_Y (GC,IGC,NN,NDIM)
C
C     This MUNGs the co-ordinates into circle in the x-z planes
C     points on the x/z axes are unaffected .. the rest get shrunk!
C
      REAL GC(IGC,*)

      DO I=1,NN
        XO = GC(1,I)
        ZO = GC(3,I)
        RADBX = MAX (ABS(XO),ABS(ZO))
        RADSP = SQRT (XO**2 + ZO**2 + 1.E-10)
        GC(1,I) = XO * RADBX / RADSP * 2.
        GC(3,I) = ZO * RADBX / RADSP * 2.
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE WRAP_SPHERE (GC,IGC,NN,NDIM)
C
C     This MUNGs the co-ordinates into a sphere
C     points on the x/y/z axes are unaffected .. the rest get shrunk!
C
      REAL GC(IGC,*)

      DO I=1,NN
        XO = GC(1,I)
        YO = GC(2,I)
        ZO = GC(3,I)
        RADBX = MAX(ABS(XO),ABS(YO),ABS(ZO))
        RADSP = SQRT (XO**2 + YO**2 + ZO**2 + 1.E-10)
        GC(1,I) = XO * RADBX / RADSP * 2.
        GC(2,I) = YO * RADBX / RADSP * 2.
        GC(3,I) = ZO * RADBX / RADSP * 2.
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE R_OPTIONS (IO,IOPS)
C
C     This reads in various program control options as an integer array
C     .. such as Output options for DISPS,ACTIONS,STRESSES etc.
C     .. and control such as Max.# of iterations, 'No-tension' supports
C         flag >limit / echo 'old'/'new' ??
      INTEGER IOPS(*)
      DO I=1,9999
    1   READ(IO,*,IOSTAT=IOS) IOP,IVAL
        CALL IN_TEST(IO,IOS,*1,*999)
        IOPS(IOP) = IVAL
      ENDDO
      CALL MYERROR (2,'-- should never get to this line') !-code = 3 ?
  999 RETURN
      END
C-----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE WR_DANPLOT(IO,FILE,GC,IGC,NUMS,INUMS,NDIM,NN,NEL)
C
C     This writes out the mesh in 'Danplot' format?
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*), IO
      CHARACTER FILE*(*)
      WRITE(IO,'(I5)') NDIM,NN
      DO I=1,NN
        WRITE(IO,'(I5,3F12.5)') I,(GC(J,I),J=1,NDIM)
      ENDDO
      WRITE(*,'(A,I4,A)')  '>> ',NN,' nodes written'
      WRITE(IO,'(I5)') NEL
      DO IEL=1,NEL
        NEN = NUMS(INUMS,IEL)
        IMAT= NUMS(INUMS-3,IEL)
        WRITE(IO,'(99I6)') IEL,NEN,(NUMS(J,IEL),J=1,NEN),IMAT
      ENDDO
      WRITE(*,'(A,I4,A)')  '>> ',NEL,' elements written'
      RETURN
      END  
C----------------------------------------------------------------------
      SUBROUTINE WR_OFF (IO,FILE,GC,IGC,NUMS,INUMS,NDIM,NN,NEL)
C
C     This writes out the mesh in 'Object_File_Format' format
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*), IO

      CHARACTER FILE*(*)
      WRITE(IO,'(I5)') NN,NEL,12345
      DO I=1,NN
        WRITE(IO,'(3F12.5)') (GC(J,I),J=1,NDIM)
      ENDDO
      WRITE(*,'(A,I4,A)')  '>> ',NN,' nodes written'
      DO IEL=1,NEL
        NEN = NUMS(INUMS,IEL)
        IMAT= NUMS(INUMS-3,IEL)
        WRITE(IO,'(99I6)') NEN,(NUMS(J,IEL),J=1,NEN)
      ENDDO
      WRITE(*,'(A,I4,A)')  '>> ',NEL,' elements written'
      RETURN
      END  
C----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------


c     options (dreal, intl)
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
c  18-4-93  : munged WIDTHS in DANBLOCKS to give true 3D from 2D
c  3- 4-92  MESH4B : 'circle_the_square'
c  3- 4-93  added 2D-to-3D  to *DANBLOCKS  
c  1- 4-93   >>>> added *DANBLOCKS ! <<<<
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
     +              ,COL*7   !-- ANSI colours
     +             ,COL1*7   !-- ANSI colours
     +             ,COL2*7   !-- ANSI colours
c    +      ,TITLES(10)*80   ! The main titles
c    +            ,LINE*255  ! 'generic' line of text

      CHARACTER  CMNAM@*100  ! command line function
c    +          ,FDATE@*40   ! Date function
c    +           ,TIME@*8    ! Time function

      INTEGER U0             ! .INP : the original data file
     +       ,U1             !  --- : the current data file being read
c    +       ,U6             ! .DP2 : The Output 'keyword' mesh file
     +       ,U7             ! .PL  : The Danplot output file
      EXTERNAL COL
c------------------ DATA stataments ------------------------------------
      DATA    U0/40/,  U7/10/
      DATA NEL,NN/2*0/            !-- initially NO elements etc.
      DATA NDIM/0/                !-- default is 0-dimensions !

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


C---------------------- Control Options --------------------------------
      IF (KEYWORD.EQ.'*EOF') THEN
        GOTO 888                          !-- exit the program nicely

      ELSEIF (KEYWORD.EQ.'*CONTROL') THEN
c .... eg. debug output options 

      ELSEIF (KEYWORD.EQ.'*TWO_DIMENSIONAL') THEN         !- N/A
        NDIM = 2

      ELSEIF (KEYWORD.EQ.'*THREE_DIMENSIONAL') THEN       !- N/A
        NDIM = 3

C-------------------------- reading-in the mesh ------------------------
c... also acts to 'add-in' elements from another file
      ELSEIF (KEYWORD.EQ.'*IMPORT_OFF') THEN
        CALL R_OFF (U1,GC,IGC,NN,NDIM,NUMS,INUMS,NEL)

      ELSEIF (KEYWORD.EQ.'*IMPORT_PL') THEN               !- synonym
        CALL R_DANPLOT (U1,GC,IGC,NDIM,NN,NUMS,INUMS,NEL)

      ELSEIF (KEYWORD.EQ.'*READ_DANPLOT') THEN            !- synonym
        CALL R_DANPLOT (U1,GC,IGC,NDIM,NN,NUMS,INUMS,NEL)

      ELSEIF (KEYWORD.EQ.'*BASIC_SLOPE_MESH') THEN
        CALL R_READGM (U1,GC,IGC,NN,NDIM,NUMS,INUMS,NEL)

      ELSEIF (KEYWORD.EQ.'*ELEMENT_TYPES') THEN    !- as IP,IQ,IS for now
        CALL R_GETELT(U1,NUMS,INUMS,NEL,NDIM)

C---------------------- removing elements + nodes ----------------------
      ELSEIF (KEYWORD.EQ.'*DELETE_MATERIALS') THEN
        CALL DEL_MATS    (U1,GC,IGC,NN,NUMS,INUMS,NEL,P,NDIM)

      ELSEIF (KEYWORD.EQ.'*DELETE_ORPHAN_NODES') THEN
        CALL DEL_O_NODES (U1,GC,IGC,NN,NUMS,INUMS,NEL,P,NDIM)

      ELSEIF (KEYWORD.EQ.'*DELETE_COINCIDENT_NODES') THEN
        CALL DEL_COIN_NODES (U1,GC,IGC,NN,NUMS,INUMS,NEL,P,NDIM)

C----------------------- modifing elements -----------------------------
      ELSEIF (KEYWORD.EQ.'*DANBLOCKS') THEN
        CALL DANBLOCKS (U1,GC,IGC,NN,NDIM,NUMS,INUMS,NEL,P)

      ELSEIF (KEYWORD.EQ.'*4NQS_FROM_4_3NTS') THEN   !-- 'teapot.g'
      CALL MUNG_4NQ  (U1,GC,IGC,NN,NDIM,NUMS,INUMS,NEL,P)

      ELSEIF (KEYWORD.EQ.'*MIRROR_EVEN_ELEMENTS') THEN   !- mung (OBS)
      CALL MIRROR_EVENS (U1,GC,IGC,NN,NDIM,NUMS,INUMS,NEL,P)

c----------------------- moving the mesh -------------------------------
      ELSEIF (KEYWORD.EQ.'*MOVE_NODES') THEN        ! from --> to
        CALL R_NODMOV(U1,GC,IGC,NDIM,NN)

      ELSEIF (KEYWORD.EQ.'*SHIFT_MESH') THEN        ! by a constant
        CALL SHIFT (U1,GC,IGC,NN,NDIM)

      ELSEIF (KEYWORD.EQ.'*SCALE_MESH') THEN        ! by 3 factors
        CALL SCALE (U1,GC,IGC,NN,NDIM)

      ELSEIF (KEYWORD.EQ.'*X_TO_Y_TO_Z') THEN       ! toggle up-dir
        CALL X2Y2Z  (GC,IGC,NN,NDIM)

      ELSEIF (KEYWORD.EQ.'*X_TO_Y') THEN            ! flip x and y
        CALL X2Y    (GC,IGC,NN,NDIM)

      ELSEIF (KEYWORD.EQ.'*WRAP_AROUND_Y') THEN     ! treat Z as é
        CALL WRAP_Y (GC,IGC,NN,NDIM)

      ELSEIF (KEYWORD.EQ.'*CIRCLE_A_SQUARE') THEN   ! move nodes onto an arc
        CALL WRAP_CIRCLE_SQUARE (U1,GC,IGC,NN,NDIM)

      ELSEIF (KEYWORD.EQ.'*WRAP_CIRCLE_Y') THEN     ! move all nodes onto arcs
        CALL WRAP_CIRCLE_Y (GC,IGC,NN,NDIM)

      ELSEIF (KEYWORD.EQ.'*WRAP_SPHERE') THEN       ! move all nodes onto spheres
        CALL WRAP_SPHERE (GC,IGC,NN,NDIM)

      ELSEIF (KEYWORD.EQ.'*WRAP_TWIST_Y') THEN      ! twist x-z planes up Y
        CALL WRAP_TWIST_Y (U1,GC,IGC,NN,NDIM)

C----------------------- Outputing the mesh ----------------------------
c... Now is the time to define the NEW Danplot styles for node and 
c...  element tables  :
c  ie. with/without node #
c      element #, material property number, NEN+NDIME+ITYPE ?
c    so really need column headers to tell me what is happening next

      ELSEIF (KEYWORD.EQ.'*WRITE_DANPLOT') THEN
        OPEN(U7,FILE=FILE(1:INDEX(FILE,'.D'))//'MSH')
        CALL WR_DANPLOT (U7,GC,IGC,NDIM,NN,NUMS,INUMS,NEL)

      ELSEIF (KEYWORD.EQ.'*WRITE_OFF') THEN
        OPEN(U7,FILE=FILE(1:INDEX(FILE,'.D'))//'GEO')
        CALL WR_OFF (U7,FILE,GC,IGC,NUMS,INUMS,NDIM,NN,NEL)

C----------------------------------------------------------------------
        ELSE
          CALL MYERROR(1,'unknown keyword:'//KEYWORD)
        ENDIF      !--- end of possible keywords

C----------------------------------------------------------------------
      WRITE(*,'(T40,A,I5,A,I5,A)')                 !- a debug option
     +     'Now:',NN,' Nodes and',NEL,' Elements'
      ENDDO
C------------------------- end of program ----------------------------
  888 STOP '> > >       Program completed OK  :-)       < < <' 
  999 STOP  '***** CRASHED OUT ******'
      END

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      FUNCTION COL(icode)
C
C     test function to output ANSI colour control strings
C     ... doesn't really work ! try statement function instead ?

      CHARACTER COL*7
      WRITE (COL,'(A,I2.2,A)')    CHAR(27)//'[1;',ICODE,'m'
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GET_KEYWORD (U1,U0,KEYWORD)
C
C     This gets the next keyword from the data file
C     '*EOF' is returned when all data is exhausted 
C      ie. U1  is reduced back to below U0 = the original data file
C
      CHARACTER KEYWORD*(*), LINE*80
      INTEGER U1,U0
    1 READ (U1,'(A)',IOSTAT=IOS) LINE
      CALL IN_TEST(U1,IOS,*1,*999)  
      IF (U1.LT.U0) GOTO 999
      IF (LINE(1:1).NE.'*') GOTO 1
      KEYWORD  = LINE(1:INDEX(LINE,' '))   !-- strip off any comments
      CALL UPCASE (KEYWORD)                 !-- Uppercase is easier to handle
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

c      ELSEIF (LINE(1:1).EQ.'!') THEN !----An 'echoed' comment line -----
c        WRITE (*,'(A)') LINE            ! ( only killed for WATFOR !)
c        RETURN 1
      ELSEIF (LINE(1:1).EQ.'%') THEN !----An 'echoed' comment line -----
        WRITE (90,'(A)') LINE        !       to the '.PLT' file
        RETURN 1
      ELSEIF (LINE(1:1).eq.'*') THEN !------- A new Keyword -----------
        BACKSPACE (IO)                  ! ie. point at the '*Keyword'
        RETURN 2                        ! 'HARD' to abort 'current read'
      ELSEIF (LINE(1:1).eq.'<') THEN !------- read_from_another_file ---
        IO = IO + 1
c... carful about the use of spaces
        OPEN (IO,FILE = LINE(2:INDEX(LINE,'  ')),IOSTAT=IOSTAT)
        IF (IOSTAT.NE.0) THEN       
          CALL MYERROR (2,'missing data file-'//LINE)
          RETURN 2  !-- crash out   ( Dead anyway )
        ENDIF
        RETURN 1   !--- this line was missing ! 14-9-92
c- - - - - - - - - - - - - - - - - - - 
      ELSE                    !--------- 'rubbish' found in file -------
        CALL MYERROR (2,'Rubbish found in the data file ='//LINE)
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
        RETURN 2               ! ie. ALL data has been exhausted
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
C-----------------------------------------------------------------------
C    > > > >   Basic Mesh Reading / Writing routines  < < < < <  
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C.. note the OLD 'SLOPE' stuff is really obsolete and should be put out
C   to graze as it is too long and old .. still allow its use however
C   ( maybe directly into FE5LIB test jobs for compactness )

C-----------------------------------------------------------------------
      SUBROUTINE R_OFF (IO,GC,IGC,NN,NDIM,NUMS,INUMS,NEL)
C
C     This imports (adds) a mesh in 'Object File Format'
C     ( These are 2d polygons (mostly 3 node triangles) in a 3d mesh
C     .. no need to test for comments as they never have any !
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*), NUM(30)

      NN_OLD  = NN
      NEL_OLD = NEL
    1 READ (IO,*,IOSTAT=IOS)  NN_NEW,NEL_NEW,NEDGES
      CALL IN_TEST(IO,IOS,*1,*999)
      READ (IO,*)  ((GC(J,I+NN_OLD),J=1,3),I=1,NN)
      WRITE(*,'(A,I4,A)')  '>> ',NN_NEW,' nodes read'
      NN = NN_OLD + NN_NEW


      NDIME = 2   !-- always 2D facets
      ITYPE = 9   !-- use type 9 for OFF format (ie. do disps, and 'flat')

      DO IEL=1,NEL
        READ (IO,*) NOD,(NUM(J),J=1,NOD)

        DO J=1,NOD
          NUM(J) = NUM(J) + NN_OLD   !-- offsets
        ENDDO
        CALL ADD_ELEMENT 
     +  (NUMS,INUMS,NEL_OLD+IEL, NUM, NOD,NDIME,ITYPE, IMAT,IEL,1,1)

      ENDDO
      CLOSE (IO)
      IO = IO - 1     !-- exit from this file  :-)
      WRITE(*,'(A,I4,A)')  '>> ',NEL,' elements read'
  999 CONTINUE
      NEL= NEL_NEW
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
        NOD = NUMS(INUMS,IEL)
        IMAT= NUMS(INUMS-3,IEL)
        WRITE(IO,'(99I6)') NOD,(NUMS(J,IEL),J=1,NOD)
      ENDDO
      WRITE(*,'(A,I4,A)')  '>> ',NEL,' elements written'
      RETURN
      END  

C-----------------------------------------------------------------------
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
c.. now 'GET_ELEMENT'
          NOD = NUMS (INUMS,IEL)  
          IFLAG = 1
          DO J=1,NDIM
            IPQR = NUMS(INUMS-3-J,IEL)
            IF (IPQR.LT.XF(J).OR.IPQR.GT.XT(J) ) IFLAG = 0   !-Not OK
          ENDDO
c..... I don't like this .. better to GET_ELEMENT then PUT_ELEMENT
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
C-----------------------------------------------------------------------
C    > > > > >  Various Geometry MUNGing routines   < < < < <
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      SUBROUTINE MUNG_4NQ (IO,GC,IGC,NN,NDIM,NUMS,INUMS,NEL,P)
C
C     This glues 4 consecutive 3nt's into a 4nq (eg. TEAPOT.G)
C
C .. nice to use GET/PUT_ELEMENT but would make this a bit too long

      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*), P(*)
      IC = 0
      DO IEL=1,NEL,4
        IC = IC +1
        NUMS (1,IC) = NUMS(1,IEL  )     !- copy
        NUMS (2,IC) = NUMS(2,IEL  )     !- copy
        NUMS (3,IC) = NUMS(1,IEL+2)     !- move
        NUMS (4,IC) = NUMS(2,IEL+2)     !- move
        NUMS (INUMS,IC) = 4
      ENDDO
      NEL = IC
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE MIRROR_EVENS (IO,GC,IGC,NN,NDIM,NUMS,INUMS,NEL,P)
C
C     This reverses the order of the nodes of alternate elements (eg HEAD.G)
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*), P(*)
      DO IEL=1,NEL,2
        NOD = NUMS(INUMS,IEL)
        DO I1 = 1,NOD/2
          I2 = NOD+1 - I1
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
C        ( 'old' routine from my 1991 Spud-can days)
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
      SUBROUTINE WRAP_CIRCLE_SQUARE (IO,GC,IGC,NN,NDIM)
C
C     This MUNGs the co-ordinates lying on a given square into a circle
C
      REAL GC(IGC,*), XC

      DO ILOOP = 1,99999

    1 READ (IO,*,IOSTAT=IOS) XC,YC, R_SQ,R_C
      CALL IN_TEST (IO,IOS,*1,*999)


      IC = 0
      DO I=1,NN
        XP = GC(1,I)
        YP = GC(2,I)

        XR = XP - XC      !- the 2 'radii'
        YR = YP - YC 

        R_HI = MAX ( ABS(XR),ABS(YR) )     !- must = R_SQ

        IF ( ABS(R_HI-R_SQ).LT.1.E-4 ) THEN

          IC = IC + 1
          RAD = SQRT (XR**2 + YR**2)
          GC(1,I) = XC + R_C * XR/RAD      ! (note Dir. cos's)
          GC(2,I) = YC + R_C * YR/RAD 

        ENDIF
      ENDDO    !- the nodes

      ENDDO
  999 WRITE(*,'(I4,A,I5,A)')
     +  ILOOP,' :',IC,' nodes moved'
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
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C    > > > > >  Element moving/deleting routines  < < < < <
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      SUBROUTINE DANBLOCKS (IO,GC,IGC,NN,NDIM,NUMS,INUMS,NEL,P)
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     This Major routine will subdivide ANY element ( no tri's yet) into
C     smaller elements of ANY type
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     3-4-93  now can fill 8nq's with 20nb !!  :-)
C     ... however first we need to set the z-dimensions to zero
C     ... then upadte NDIM to 3
C
C.... It would be nice to give some user-feedback as this can take a 
C     long time

      PARAMETER (MAX_NOD = 27,MDIM=4)
      PARAMETER (ICOORD= MAX_NOD,IDER=MDIM)

      REAL      GC (IGC,*)        !- Nodal coord database
     +     ,WIDTHS (0:20,5)       !- input sub-element widths
     +   ,WIDTHS_L (0:20,5)       !- sub-element widths in locals
     +  ,COORD_OLD (ICOORD,MDIM)  !- coords of the parent element
     +    ,COORD_L (ICOORD,3)     !-  local coords of the daughter elements
     +    ,COORD_G (MDIM)         !- global coords of a daughter node
     +       ,SAMP (1,MDIM)       !- sampling point of the 'new' nodes
     +        ,FUN (MAX_NOD)      !- shape functions of the parent
     +        ,DER (IDER,MAX_NOD) !- sf derivs of the parent (dummy)

      INTEGER  NUMS (INUMS,*)     !- Element database
     +         ,NUM (MAX_NOD)     !- node numbers of the new elements
     +           ,P (*)           !- node pointers (unused)
     +           ,L (MDIM)        !- toggles for the sub-directions
     +         ,NXE (MDIM)        !- number of sub-elements in each direction

      NEL_NEW  = NEL      !- new element numbers (incremented)
      NN_NEW   = NN       !- new node numbers    (updated)

      DO ILOOP = 1,99999

    1 READ (IO,*,IOSTAT=IOS) JMAT , NDIME_NEW,NOD_NEW,ITYPE  
      CALL IN_TEST (IO,IOS,*1,*999)

c------- go into '3d'  if we were only in '2d' before -------
c.. maybe as a routine of its own ? (as x/y to x/y/z mappings ?)
        IF (NDIME_NEW.GT.NDIM) THEN
          DO I=1,NN
            DO J=NDIM+1,NDIME_NEW
              GC(J,I) = 0.
            ENDDO
          ENDDO
          NDIM = NDIME_NEW
        ENDIF
c------------------------------------------------------------

c---------  get the local coords of the new elems. ----------
      CALL WTHATN (NOD_NEW,NDIME_NEW,ITYPE,COORD_L,ICOORD)

      DO I=1,NDIME_NEW       !-- read the edge spacings
    2   READ (IO,*,IOSTAT=IOS)   NXE(I),(WIDTHS(J,I),J=1,NXE(I))
        CALL IN_TEST (IO,IOS,*2,*999)   

        WIDTHS(0,I) = 0.    !- null the first coord in the list

C.... careful .. I dont want to 'non-dimensionalize' if 2d-> 3d

        WIDTH_TOT = 0.
        DO J=0,NXE(I)    !-- get cumultive widths
          WIDTH_TOT   = WIDTH_TOT + WIDTHS(J,I)
          WIDTHS_L(J,I) = WIDTH_TOT 
        ENDDO
        DO J=0,NXE(I)   
          WIDTHS_L(J,I) = WIDTHS_L(J,I) / WIDTH_TOT !-- scale to 0. -> +1.
          WIDTHS_L(J,I) = WIDTHS_L(J,I) * 2. -1.    !-- into locals -1 -> +1
c... triangles ???
        ENDDO
      ENDDO

C-------------------- loop the 'original' elements ---------------------
      DO IEL = 1,NEL
c..  I could pre-scan to get a list of all relevant elements ?
        CALL GET_ELEMENT (NUMS,INUMS,IEL, NUM, 
     +          NOD_OLD,NDIME_OLD,ITYPE,IMAT,IUSER1,IUSER2,IUSER3)

        IF (JMAT.EQ.IMAT) THEN               !-- found a material

          CALL GET_COORD (NUM,NOD_OLD,NDIM,GC,IGC,COORD_OLD,ICOORD)  

c-----------------------------------------------------------------------

          CALL ADD_ELEMENT (NUMS,INUMS,IEL, NUM, NOD_OLD,NDIME_OLD
     +      ,ITYPE  ,-1,IUSER1,IUSER2,IUSER3)    !- kill the parent

c------------ set up the indecs to loop ALL daughter elements ----------
c... careful with Triangles (cf. DANPLOT's 'sub-facets' )
          N_SUB_EL = 1
          DO IPOS_L=1,NDIME_NEW
            L(IPOS_L) = 0                      !-- null the pointers
            N_SUB_EL = N_SUB_EL * NXE(IPOS_L)  !-- 
          ENDDO

C--------------- loop the in-plane daughter elements -------------------
c.. then inside loop the 'new' dimension

          DO I_SUB_EL = 1,N_SUB_EL
      
c-------- now get the local (hence global) coords of the new nodes -----
          DO INODE = 1,NOD_NEW
            DO J=1,NDIME_NEW
              SAMP(1,J)   = 
     +        WIDTHS_L(L(J),J) + (WIDTHS_L(L(J)+1,J) - WIDTHS_L(L(J),J))
     +     * ( COORD_L(INODE,J)+1.)/2.
            ENDDO

            CALL GSF (NDIME_OLD,NOD_OLD,ITYPE ,DER,IDER,FUN ,SAMP,1,1)

            DO J=1,NDIME_NEW
              IF (J.LE.NDIME_OLD) THEN
c------ first the 'in-plane' coordinates (eg X, Y) -------
                COORD_G(J) = 0.         !- ie. MTV_MULT  :-)
                DO K=1,NOD_OLD
                  COORD_G(J) = COORD_G(J) + FUN(K) * COORD_OLD(K,J)
                ENDDO   !- nodes
              ELSE
c----- then the 'out-of-plane' coordinates (eg Z) --------
                COORD_G(J) = WIDTHS (L(J),J)   !- direct copy :-)
              ENDIF

C----------------------------------------------------------
            ENDDO   !- the dimensions of the new element
               
c----------------- build the nodes of the new element ------------------ 
            CALL FIND_OR_ADD_NODE (GC,IGC,NDIM,COORD_G,1,NN_NEW,JNODE)
            NUM(INODE) = JNODE
          ENDDO     !-- (loop nodes of the new element)

c--------------------- add-in the new element ---------------------------
          IF (N_SUB_EL.EQ.1) THEN
            IEL_NEW = IEL                !- can replace the old one
          ELSE
            NEL_NEW = NEL_NEW + 1        !- or build the new elements
            IEL_NEW = NEL_NEW            !- at the end of the list
          ENDIF

          CALL ADD_ELEMENT 
     +    (NUMS,INUMS,IEL_NEW, NUM, NOD_NEW,NDIME_NEW,ITYPE, IMAT,
     +     IUSER1,IUSER2,IUSER3)

C------------- end of this new element .. increment pointers -----------
C.. not quite like this for triangles
C.. anyway maybe better as its own routine

          IPOS_L = NDIME_NEW
  123     L(IPOS_L) = L(IPOS_L) + 1               !-- increment
          IF (L(IPOS_L) .GE.NXE(IPOS_L)) THEN     !-- the pointers
            L(IPOS_L) = 0
            IPOS_L = IPOS_L -1
            IF (IPOS_L.GT.0) GOTO 123             !-- next toggle
          ENDIF
          ENDDO           !-- the new elements
        ENDIF          !-- only this material
      ENDDO         !-- original element loop

      WRITE(*,'(I4,A,I5,A,I4,A)')
     +  ILOOP,' :',NEL_NEW-NEL,' subdivided into', N_SUB_EL,' each'
      ENDDO         !-- end of main loop

  999 CONTINUE

      NEL = NEL_NEW
      NN  = NN_NEW
      END            

C-----------------------------------------------------------------------
      SUBROUTINE DEL_MATS (IO,GC,IGC,NN,NUMS,INUMS,NEL,P,NDIM)
C
C     This will delete any elements of a given material type
C     ... usually those of TYPE = 0  'air' elements
C  --> better to sub-call the actual delete routine so I can call from 
C      elsewhere  (eg. DANBLOCKS)
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*), P(*)

      DO KK = 1,999
    1   READ (IO,*,IOSTAT=IOS) ITYPE
        CALL IN_TEST (IO,IOS,*1,*999)
        NEL1 = NEL    
        IC = 0
        DO IEL=1,NEL
c.. please use GET_ELEMENT to get MAT the PUT_ELEMENT to re-save
          IMAT = NUMS (INUMS-3,IEL)
          IF (IMAT.NE.ITYPE) THEN               !-- If a 'good' element ..
            IC = IC + 1 
            DO J = 1,INUMS     !- only needed if IC is NOT = IEL :-)
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

      NN_OLD  = NN           !-- original mesh sizes  (for print-out)

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
      CALL RESEQ_P (P,NN,NN_NEW)           ! == a 1D NF array

      CALL UPDATE_GC (GC,IGC,NDIM,NN, P)
      CALL UPDATE_NUMS (NUMS,INUMS,NEL, P)

      WRITE(*,'(A,I4,A)')  '>> ',NN_OLD-NN_NEW,' nodes deleted'
      NN = NN_NEW
      END

C-----------------------------------------------------------------------
      SUBROUTINE DEL_COIN_NODES (IO,GC,IGC,NN,NUMS,INUMS,NEL,P,NDIM)
C
C     This will 'condense' all nodes that share the same position
C     using FIND_NODE to get the node number at a given co-ord
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
C-----------------------------------------------------------------------
C  > > >   A set of routines for finding / re-ordering node #s   < < < 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

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
      END

C-----------------------------------------------------------------------
      SUBROUTINE FIND_NODES (NN,GC,IGC,POINT,NDIM,NUMLIST,INUMLIST,IC)
C
C      *OBSOLETE as routines can multiply call FIND_NODE
C     Given a co-ordinate, this returns its node number(s) in NUMLIST
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
      SUBROUTINE FIND_OR_ADD_NODE (GC,IGC,NDIM,XY,IFROM,ITO,INODE)
C
C     This calls FIND_NODE to get the node number at a given coordinate
C     if now node was found then it is added to the end of the list

      REAL GC(IGC,*), XY(*)
        CALL FIND_NODE (GC,IGC,NDIM,XY,IFROM,ITO,INODE)
        IF (INODE.LE.0) THEN
          ITO = ITO + 1
          DO J=1,NDIM
            GC(J,ITO) = XY(J)
          ENDDO
          INODE = ITO
        ENDIF
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
C     held in P  ... This will ONLY work if the P's are in ASCENDING 
C     order  .. This must NOT be used if the P are sorted (eg. by 
C     distance form a point or by element order
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
C     This re-orders the node-steering in NUMS based on the NEW node order
C     held in P
C
C.. should really use GET_ELEMENT and PUT_ELEMENT  !

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
C----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
c     include 'c:\users\dan\lib\shape.for'
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C   > > > > > > 'OLD' style basic mesh forming routines < < < < < < <
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
c.. should really be archived out of the way !
C-----------------------------------------------------------------------
      SUBROUTINE R_READGM (IO,GC,IGC,NN,NDIM,NUMS,INUMS,NEL)
C
C     this routine reads the element geometry puts it into 'long' form
C     19-1-93 moved from L1A.FOR to MESH1.FOR for mesh generation
C
      REAL GC(IGC,*), COORD(20,3)
      REAL    TOP(90), BOT(90), DEPTH(90), BRDTH(90)
      INTEGER NUM(20), NUMS(INUMS,*)

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

c....... hmm I could even have 1d elements !
          IF (NDIM.EQ.2) THEN
            CALL SLOGE2 (IP,IQ,NXS,NYS,NYE,TOP,BOT,DEPTH,COORD,20,NUM)
          ELSEIF (NDIM.EQ.3) THEN           !-- 8nq only
            CALL SLOGET (IP,IQ,IS,NXS,NYS,NYE,NZE,TOP,BOT,DEPTH,BRDTH,
     +      COORD,20,NUM,NOD)      !-- 8nb/14nb/20nb
          ELSE 
            CALL MYERROR(2,'NDIM not valid (not 2 or 3)')  
          ENDIF
        NDIME = NDIM  !------------- update the element table ----------
        IMAT  = 1
        ITYPE = 1
        CALL ADD_ELEMENT 
     +  (NUMS,INUMS,IEL, NUM, NOD,NDIME,ITYPE, IMAT,IP,IQ,IS)

         DO I=1,NOD        !----- update the coordinate table --------- 
           INODE = NUM(I)
           NN = MAX (INODE,NN)
           DO J=1,NDIM
             GC (J,INODE) = COORD(I,J)
           ENDDO
         ENDDO

          ENDDO         ! loop       -z
        ENDDO          ! the        -y
      ENDDO           ! elements   -x
      NEL = IEL
      WRITE (*,'(A,I5,A,I5)') '>> NN=',NN, ' NEL=',NEL
      RETURN
  999 CALL MYERROR(2, 'mesh reading failed')
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SLOGET (IP,IQ,IS,NXS,NYS,NYE,NZE,TOP,BOT,DEPTH,BRDTH
     +       ,COORD,ICOORD,NUM,NEN)
C
C    header to SLOGE3 (20nb), SLOG14 (14nb) and SLOG08 (8nb/11nb)
C    .. needs fixing so that it can call SLOGE2 also 
C
      REAL TOP(*),BOT(*),DEPTH(*),BRDTH(*),COORD(ICOORD,*)
      INTEGER NUM(*)
      IF (NEN.EQ.20) THEN
      CALL SLOGE3 (IP,IQ,IS,NXS,NYS,NYE,NZE,TOP,BOT,DEPTH,BRDTH
     +       ,COORD,ICOORD,NUM)
      ELSEIF (NEN.EQ.14)THEN
      CALL SLOG14 (IP,IQ,IS,NXS,NYS,NYE,NZE,TOP,BOT,DEPTH,BRDTH
     +       ,COORD,ICOORD,NUM)
      ELSEIF (NEN.EQ. 8 .OR. NEN.EQ.11) THEN
      CALL SLOG08 (IP,IQ,IS,NXS,NYS,NYE,NZE,TOP,BOT,DEPTH,BRDTH
     +       ,COORD,ICOORD,NUM)
      IF (NEN.EQ.11) THEN
        DO I=9,11
          NUM(I) = I    !-- for 11nb's hack the last 3 nodes to 'dummies'
        ENDDO
        ENDIF
      ELSE
        PRINT*,'  OOPS not a 8,11,14 or 20 node brick (SLOGET)'
        STOP
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------

      SUBROUTINE SLOGE2 (IP,IQ,NXS,NYS,NYE,TOP,BOT,DEPTH,COORD,ICOORD
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
C-----------------------------------------------------------------------
      SUBROUTINE SLOGE3 (IP,IQ,IS,NXS,NYS,NYE,NZE,TOP,BOT,DEPTH,BRDTH
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
      IF (IQ.LE.NYS.AND.IP.GT.NXS) THEN
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
     +       - 3*NZE*NYEH -2*(NZE+NYEH) -1
     +       + (IP-NXS-1)*(4*NYEH*NZE+3*(NYEH+NZE)+2)
      ENDIF
      N1=3*NZE*NYEH+2*(NZE+NYEH)+1
      N2=N1+(NZE+1)*(NYEH+1)
      NUM(15) = BASE+ (3*NZE+2)*(IQ1-1)+ 2*(IS-1)+ 1
      NUM(10) = NUM(15)+ (2*NZE+1)- (IS-1)
      NUM( 3) = NUM(15)+ (3*NZE+2)
      NUM(16) = BASE+N1+ (NZE+1)*(IQ1-1)+  (IS-1)+ 1
      NUM( 4) = NUM(16)+ (NZE+1)
      NUM(17) = NUM(15)+ N2
      NUM(11) = NUM(10)+ N2
      NUM( 5) = NUM( 3)+ N2
      DO 5,I=1,12
    5 NUM(ITO(I)) = NUM(IFROM(I)) + 1
C     WRITE(*,'(20I4)') (NUM(I),I=1,20)
      IF (IQ.LE.NYS) THEN
        FAC1 = (BOT(IP  )-TOP(IP  ))/ (DEPTH(NYS+1)-DEPTH(1))
        FAC2 = (BOT(IP+1)-TOP(IP+1))/ (DEPTH(NYS+1)-DEPTH(1))
        COORD( 3,1) = TOP(IP  )+(DEPTH(IQ+1)-DEPTH(1))*FAC1
        COORD(15,1) = TOP(IP  )+(DEPTH(IQ  )-DEPTH(1))*FAC1
        COORD(17,1) = TOP(IP+1)+(DEPTH(IQ  )-DEPTH(1))*FAC2
        COORD( 5,1) = TOP(IP+1)+(DEPTH(IQ+1)-DEPTH(1))*FAC2
      ELSE
        COORD( 3,1) = BOT(IP)
        COORD(15,1) = BOT(IP)
        COORD(17,1) = BOT(IP+1)
        COORD( 5,1) = BOT(IP+1)
      ENDIF
      COORD( 3,2) = DEPTH(IQ+1)
      COORD(15,2) = DEPTH(IQ)
      COORD(15,3) = BRDTH(IS  )
      COORD(13,3) = BRDTH(IS+1)
      COORD(13,1) = COORD(15,1)
      COORD(19,1) = COORD(17,1)
      COORD( 1,1) = COORD( 3,1)
      COORD( 7,1) = COORD( 5,1)
      COORD(17,2) = COORD(15,2)
      COORD(13,2) = COORD(15,2)
      COORD(19,2) = COORD(15,2)
      COORD( 5,2) = COORD( 3,2)
      COORD( 1,2) = COORD( 3,2)
      COORD( 7,2) = COORD( 3,2)
      COORD(17,3) = COORD(15,3)
      COORD( 3,3) = COORD(15,3)
      COORD( 5,3) = COORD(15,3)
      COORD(19,3) = COORD(13,3)
      COORD( 1,3) = COORD(13,3)
      COORD( 7,3) = COORD(13,3)
      DO 8,I=1,12
      DO 8,J=1,3
    8 COORD(IE(2,I),J) = .5*(COORD(IE(1,I),J)+COORD(IE(3,I),J))
      RETURN
      END
C-----------------------------------------------------------------------
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
      IF (IQ.LE.NYS.AND.IP.GT.NXS) THEN
        PRINT*,'OUT OF MESH RANGE..... ie this element aint ''ere'
        STOP
      ENDIF
      IF (IP.LE.NXS)NYEH = NYE
      IF (IP.GT.NXS)NYEH = NYE-NYS
      F1 = (NYEH+1) * (NZE+1)
      F2 = (NYE +1) * (NZE+1)
      IF (IP.LE.NXS) THEN
        IQ1=IQ
      BASE = (IP-1) * F2
      ELSE
        IQ1 = IQ-NYS
        BASE = NXS* F2  + (IP-1-NXS) * F1
      ENDIF
C   ------- simply step around the element,to number it-------
C     PRINT*,BASE,F1,F2
      NUM( 6) = BASE+ (NZE+1) * (IQ-1) + IS
      NUM( 5) = NUM( 6)+ 1
      NUM( 2) = NUM( 5)+ NZE
      NUM( 1) = NUM( 2)+ 1
      NUM( 7) = NUM( 6)+ F1
      NUM( 8) = NUM( 7)+ 1
      NUM( 3) = NUM( 8)+ NZE
      NUM( 4) = NUM( 3)+ 1
C     WRITE(*,'(8I4)')(NUM(I),I=1,8)
C -------- set the 4 x-values of the 'trapezoidal' faces -----
      IF( IQ.LE.NYS) THEN
        FAC1 = (BOT(IP  )-TOP(IP  ))/ (DEPTH(NYS+1)-DEPTH(1))
        FAC2 = (BOT(IP+1)-TOP(IP+1))/ (DEPTH(NYS+1)-DEPTH(1))
        COORD( 1,1) = TOP(IP  )+(DEPTH(IQ+1)-DEPTH(1))*FAC1
        COORD( 5,1) = TOP(IP  )+(DEPTH(IQ  )-DEPTH(1))*FAC1
        COORD( 8,1) = TOP(IP+1)+(DEPTH(IQ  )-DEPTH(1))*FAC2
        COORD( 4,1) = TOP(IP+1)+(DEPTH(IQ+1)-DEPTH(1))*FAC2
      ELSE
        COORD( 1,1) = BOT(IP)
        COORD( 5,1) = BOT(IP)
        COORD( 8,1) = BOT(IP+1)
        COORD( 4,1) = BOT(IP+1)
      ENDIF
      COORD( 1,2) = DEPTH(IQ+1)
      COORD( 5,2) = DEPTH(IQ  )
      COORD( 1,3) = BRDTH(IS+1)
      COORD( 2,3) = BRDTH(IS  )
C --------  copy  to nodes of the same x,y,z ------
      DO 4,I = 1,5,4
      COORD( 8-I,1) = COORD(9-I,1)
    4 COORD( 1+I,1) = COORD(I  ,1)
      DO 5,I=2,4
      COORD(I  ,2) = COORD(1,2)
    5 COORD(I+4,2) = COORD(5,2)
      COORD(4,3) = COORD(1,3)
      COORD(5,3) = COORD(1,3)
      COORD(8,3) = COORD(1,3)
      COORD(3,3) = COORD(2,3)
      COORD(6,3) = COORD(2,3)
      COORD(7,3) = COORD(2,3)
      RETURN
      END
C-----------------------------------------------------------------------
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
      IF (IQ.LE.NYS.AND.IP.GT.NXS) THEN
        PRINT*,'OUT OF MESH RANGE..... ie this element aint ''ere'
        STOP
      ENDIF
      IF(IP.LE.NXS)NYEH = NYE
      IF(IP.GT.NXS)NYEH = NYE-NYS
C     PRINT*,NYE,NZE
      F1 = 2*NYEH*NZE + NYEH+ NZE
      F2 = 2*NYE *NZE + NYE + NZE
      IF(IP.LE.NXS)THEN
        IQ1 = IQ
        BASE = (IP-1)*(2*F2+1)
      ELSE
        IQ1 = IQ-NYS
        BASE = NXS* (2*F2+1) +NYS* (2*NZE+1) + (IP-1-NXS) * (2*F1+1)
      ENDIF
C---------- simply step around the element to number it ----------------
C     PRINT*,BASE,F1,F2
      NUM( 6) = BASE+ (2*NZE+1)*(IQ1-1)+ IS
      NUM( 2) = NUM( 6)+ 1
      NUM(13) = NUM( 2)+ NZE
      NUM( 5) = NUM(13)+ NZE
      NUM( 1) = NUM( 5)+ 1
      NUM(12) = NUM( 2)+ F1
      NUM(10) = NUM(12)+ NZE
      NUM( 9) = NUM(10)+ 1
      NUM(11) = NUM( 9)+ NZE
      NUM( 8) = NUM(12)+ F1
      NUM( 4) = NUM( 8)+ 1
      NUM(14) = NUM( 4)+ NZE
      NUM( 7) = NUM(14)+ NZE
      NUM( 3) = NUM( 7)+ 1
C -------- set the 4 x-values of the 'trapezoidal' faces -----
      IF (IQ.LE.NYS) THEN
        FAC1 = (BOT(IP  )-TOP(IP  ))/ (DEPTH(NYS+1)-DEPTH(1))
        FAC2 = (BOT(IP+1)-TOP(IP+1))/ (DEPTH(NYS+1)-DEPTH(1))
        COORD( 2,1) = TOP(IP  )+(DEPTH(IQ+1)-DEPTH(1))*FAC1
        COORD( 1,1) = TOP(IP  )+(DEPTH(IQ  )-DEPTH(1))*FAC1
        COORD( 4,1) = TOP(IP+1)+(DEPTH(IQ  )-DEPTH(1))*FAC2
        COORD( 3,1) = TOP(IP+1)+(DEPTH(IQ+1)-DEPTH(1))*FAC2
      ELSE
        COORD( 2,1) = BOT(IP)
        COORD( 1,1) = BOT(IP)
        COORD( 4,1) = BOT(IP+1)
        COORD( 3,1) = BOT(IP+1)
      ENDIF
      COORD( 1,2) = DEPTH(IQ+1)
      COORD( 2,2) = DEPTH(IQ  )
      COORD( 1,3) = BRDTH(IS+1)
      COORD( 5,3) = BRDTH(IS  )
C --------  copy  to nodes of the same x,y,z ------
      DO 4,I=5,8
    4 COORD( I,1) = COORD(I-4,1)
      COORD( 3,2) = COORD( 1,2)
      COORD( 5,2) = COORD( 1,2)
      COORD( 7,2) = COORD( 1,2)
      COORD( 4,2) = COORD( 2,2)
      COORD( 6,2) = COORD( 2,2)
      COORD( 8,2) = COORD( 2,2)
      COORD( 2,3) = COORD( 1,3)
      COORD( 3,3) = COORD( 1,3)
      COORD( 4,3) = COORD( 1,3)
      COORD( 6,3) = COORD( 5,3)
      COORD( 7,3) = COORD( 5,3)
      COORD( 8,3) = COORD( 5,3)
C--------------- put face nodes at the 'mid-face' ---------------
      DO 8,I=1,6
        DO 8,J=1,3
        SUM=0
          DO 9,K=1,4
    9     SUM = SUM+COORD(IF(I,K),J)
    8   COORD(IFN(I),J) = SUM/4.
C     DO 10,I=1,14
C  10 WRITE(*,'(3F10.2)')(COORD(I,J),J=1,3)
      RETURN
      END

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C    > > > > > > > >     Element Database routines < < < < < < < < <
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      SUBROUTINE R_DANPLOT (IO,GC,IGC,NDIM,NN,NUMS,INUMS,NEL)
C
C     This imports a mesh in 'Danplot format'
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*), NUM(30)

      NN_OLD  = NN
      NEL_OLD = NEL

   1  READ (IO,*,IOSTAT=IOS) NDIM,NN_NEW
      CALL IN_TEST(IO,IOS,*1,*999)

      READ (IO,*) (ID,(GC(J,I+NN_OLD),J=1,NDIM),I=1,NN_NEW)
      WRITE(*,'(A,I4,A)')  '>> ',NN_NEW,' nodes read'
      NN = NN_OLD + NN_NEW

      NDIME = NDIM             !- always the same
      ITYPE = 1                !- all we can guess !
      READ (IO,*) NEL_NEW
      DO IEL=1,NEL_NEW
        READ (IO,*) IDUMMY,NEN,(NUM(J),J=1,NEN),IMAT

        DO J=1,NEN
          NUM(J) = NUM(J) + NN_OLD   !-- offsets
        ENDDO

        CALL ADD_ELEMENT 
     +  (NUMS,INUMS,NEL_OLD+IEL, NUM, NEN,NDIME,ITYPE, IMAT,IEL,1,1)

      ENDDO
c     CLOSE (IO)
      WRITE(*,'(A,I4,A)')  '>> ',NEL_NEW,' elements read'
      NEL= NEL_OLD + NEL_NEW
  999 CONTINUE
      END
C----------------------------------------------------------------------
      SUBROUTINE WR_DANPLOT(IO,GC,IGC,NDIM,NN,NUMS,INUMS,NEL)
C
C     This writes out the mesh in 'Danplot' format?
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*), IO
      WRITE(IO,'(I5)') NDIM,NN
      DO I=1,NN
        WRITE(IO,'(I5,3F12.5)') I,(GC(J,I),J=1,NDIM)
      ENDDO
      WRITE(*,'(A,I4,A)')  '>> ',NN,' nodes written'
      WRITE(IO,'(I5)') NEL
      DO IEL=1,NEL
        NOD  = NUMS (INUMS,IEL)
        IMAT = NUMS (INUMS-3,IEL)
        WRITE(IO,'(99I6)') IEL,NOD,(NUMS(J,IEL),J=1,NOD),IMAT
      ENDDO
      WRITE(*,'(A,I4,A)')  '>> ',NEL,' elements written'
      RETURN
      END  
C-----------------------------------------------------------------------
      SUBROUTINE ADD_ELEMENT (NUMS,INUMS,IEL, NUM, NOD,NDIME,ITYPE
     +      ,IMAT,IUSER1,IUSER2,IUSER3)
C
C     This adds an element reference to the NUMS database
C     : an element of NOD nodes in NDIME dimensions (eg 2D) of 
C       element 'type' ITYPE,
C       IUSER(1-3) are user-defined properties (eg. element groups.)
C
C... must test to see wether there is enough space for this element !

      INTEGER NUMS(INUMS,*) , NUM(*)

      DO J=1,NOD
        NUMS(J,IEL) = NUM(J)
      ENDDO

      NUMS (INUMS  ,IEL) = NOD      
      NUMS (INUMS-1,IEL) = NDIME    
      NUMS (INUMS-2,IEL) = ITYPE    
      NUMS (INUMS-3,IEL) = IMAT     
      NUMS (INUMS-4,IEL) = IUSER1   
      NUMS (INUMS-5,IEL) = IUSER2   
      NUMS (INUMS-6,IEL) = IUSER3  

      DO J=NOD+1,INUMS-7
        NUMS(J,IEL) = 0    !- null any other columns
      ENDDO
      END

C-----------------------------------------------------------------------
      SUBROUTINE GET_ELEMENT (NUMS,INUMS,IEL, NUM, NOD,NDIME,ITYPE
     +      ,IMAT,IUSER1,IUSER2,IUSER3)
C
C     This extratcs an element reference to the NUMS database
C     : an element of NOD nodes in NDIME dimensions (eg 2D) of 
C       element 'type' ITYPE,
C       IUSER(1-3) are user-defined properties (eg. element groups.)
C
      INTEGER NUMS(INUMS,*) , NUM(*)

      NOD     = NUMS (INUMS  ,IEL) 
      NDIME   = NUMS (INUMS-1,IEL) 
      ITYPE   = NUMS (INUMS-2,IEL) 
      IMAT    = NUMS (INUMS-3,IEL) 
      IUSER1  = NUMS (INUMS-4,IEL) 
      IUSER2  = NUMS (INUMS-5,IEL) 
      IUSER3  = NUMS (INUMS-6,IEL) 

      DO J=1,NOD
        NUM(J) = NUMS(J,IEL)
      ENDDO
      END

C-----------------------------------------------------------------------
      SUBROUTINE GET_COORD (NUM,NOD,NDIM,GC,IGC,COORD,ICOORD)    !- coords
C
C     This simply extracts the nodal coordinates from the GC table
C     .. nice but its much easier not to bother & do directly 

      REAL GC(IGC,*),COORD(ICOORD,*)
      INTEGER NUM(*)
      DO I=1,NOD
        DO J=1,NDIM
          COORD(I,J) = GC(J,NUM(I))
        ENDDO
      ENDDO
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------


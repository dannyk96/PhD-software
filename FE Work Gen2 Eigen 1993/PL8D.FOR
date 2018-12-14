c     OPTIONS (fullcheck,undef)

C--------------------------------------------------------------------
      PROGRAM     P L O T T E R _ 7 c
C--------------------------------------------------------------------


C     this is Dan Kidger's general purpose plotting package


C---------- revision history ----------------------------------------

C     -- Version 2              23/9/91                Dr. Dan Kidger
C     -- was PL2.FOR  ... now auto -cycling added 27/9/91
C     -- light-source shading added  and  'c=5' (?)
C     -- palete changing ','
C     -- ver 4.0 read all load steps at once
C     -- screen border toggle '*'
C     -- polygon decompostion
C  22-10-91 DO-loops changed to DO..ENDDO structures 
C  23-10-91 4 node quad in '3D' added, + color backface facets 
C     'cursor keys' redraw bug fixed/ mouse action changed
C      color map '~' option/ '.' and ',' joined
C  24-10-91 screen based menu started.
C  30-10-91 on-screen time per frame '"'
C   8-11-91 2d 'flat', more menuing...
C  20-11-91 PCX fixing.. esp. palette !
C  20-11-91 dialog from top of the screen
C  25-11-91 check mouse event BEFORE key event, add abort mid-picture
C  25-11-91 sub-facet division (via shape funcs) started...
C  27-11-91 ..continued
C   5-12-91    ..continued
C   9-12-91       .. completed!
C  12-12-91 ** Gouraud shader added for disp. contouring
C  19-12-91  'G' to 'glass' a material property 
C  19-12-91  'M' to Mirror the model.
C  19-12-91  Strain-contouring begun
C   7- 1-92  Command line arguements for filename and macros added
C   8- 1-92  <<PL7>> major menu restructuring
C  25- 3-92  PL7A .. strain contouring/ 'palette' menu/ PL7A.MNU etc.
C  21- 4-92  PL7C .. GFACE restructuring.
C  28- 4-92  add coments to ALL initial data ... / .CFG config files :-)
C  30- 4-92  (push eye posn, light posn.), picture # to CV() ?
C   4- 5-92  multiplle windows !
C  22- 5-92  'PL7E' FSTRIP changes for speed!
C  26- 5-92  transform and sort only 'visible' facets
C   1- 6-92  multiple plots.. link load#/view# to pic#
C   5- 6-92  shading fixed/ 'zebra' contours/ 'material look-up table'
C   9- 6-92  -ve 'I' options moved to (+ve) J :-)
C   9- 6-92  output options on 'O' <including laser-printer:-) >
C  22- 6-92  draw'n'save animation (via PCX/memory)
C  23- 6-92  ']' to get 'cut-out models patch.. -ve mats. = <ld+1
C            '[' to get built up models patch.. +ve mats  = >ld
C   2- 7-92  'NFF' file-input support added !
C  14 -7-92  FTN77 spurious-line-diagonal BUG fixed !

C   4-12-92 WOW.. what a long time ... Postscript Output :-)

C   30-3-93 PL8C_M  'extra lines' if '*.PL2'
C   31-3-93 properly mended GOURAUD for 'clipped' contours
C    1-4-93 reversed the default colours to black on white 

C-------- INF= Max nodes, MEL= Max Elements, MN,MDF = max elem = 20nb
C-------  IFCETS= Max. facets
      INTEGER*4 IGDISPS,IBASE, IB1,IB2,IB3
      PARAMETER ( INF=10000, MEL=5000, IGDISPS=INF*10L )
c     PARAMETER ( INF=20000, MEL=13000, IGDISPS=INF*30L )

      PARAMETER ( MN=20, MDF=3,IFCETS= 5940*2+1 ) ! for N1452.ho etc.
      PARAMETER ( MOPS=400, M_CV=50,MIMAGES=200 )

C -- array sizes, NUMS (element steering), GC,etc.(Nodal coords+disps) 

      PARAMETER (INUMS=MEL,IGCRD=3, ISC=MN)
      REAL GCOORD(3,INF), GDISPS(MDF,IGDISPS)
     +   , SCOORD(3,INF)
     +   , X(10),Y(10)   !--- misc. real arrays

      REAL  XTRA_L(4,30)      !--- for 'added' lines on the drawing

C --------- the following need their subscripts reversing ! -----
     +   , SC_F(ISC,3), SC_SF(ISC,3), DI_F(ISC,3), LC_SF(ISC,3)
     +   , ACOL(ISC), VEC(3), C_F(ISC,3), C_R(6), TIME(20)

      COMMON /DEST/IOUT_D     !-- eg. 5=postscript

C-------------- the 'contour colour' options ---------------------------
c     COMMON /CCOLS/CCOL    !-- pass to DRAW_POLY ? NO! pass to GOURAUD
      INTEGER CCOL(6)

C------ nicer if NUMS was stoered the 'other-way' too ! --------
      INTEGER NUMS(INUMS,MN+3), U3, IOP_ANIM(10), CV(M_CV)  

     +       ,FACETS(3,IFCETS), FS(ISC), TABLE(11,MOPS),MENUP(40)
     +  ,IC_MATS(20)
      LOGICAL LMOUSE, KEY_WAITING@, EXIT,BEEP,FLASH

C ------------- viewing parameters, etc. -------------------------------
      COMMON /PALETTE/PAL
      INTEGER PAL(3,0:255),PL2(IFCETS)  ! colour palette, FACET pointers
     + , ID_FACE(6)   ! 'to-draw' facets
      REAL EYE(3), COA(3), UP  (3), LIGHT(3), PM(4,4), LD,LDD
     +    ,GPT(3), SPT(3), SPT2(3)
      CHARACTER LINE*132, FILE_DAT*80, FILE2*15, FILE_PCX*12
     +  ,CKEY*1,NTXT*4 ,FILE_ROOT*5
     +  ,ESC*1, TXT(MOPS)*21, C_T*6, CMNAM@*60,PATH*60,CURDIR@*50
C-----------------------------------------------------------------------
C ----- explicit 'short-integers' variable types for ftn77 calls
      INTEGER*4 BUFFER, IMAGE(MIMAGES), VSCREEN, PL1(IFCETS), NBYTES_R
     + ,ADDRESS
      INTEGER*2 IRESX,IRESY,IRESC,IRESCP,IXM,IYM,IBUTTON, RES(6) ! hmm?
     +         ,IFAIL,KEY,  IBITPL
      REAL*4 ZBUF(IFCETS)
C-----------------------------------------------------------------------
C---------------- some initial data values -----------------------------
      DATA U3/13/, FLASH/.TRUE./,BEEP/.FALSE./,EXIT/.FALSE./
      DATA N_X_L /0/ 
C------------------- viewing parameters --------------------------------
      DATA AR      /-1./            ! aspect ratio for the screen
     +   , UP      /0.,-1.,0./      ! 'up' vector = y-axis   
     +   , FEYE    /20./            ! viewing angle          +
     +   , DEYE    /-3./            ! viewing distance       +
     +   , FACT    / 1./            ! disp scale factor      +
     +   , XML,YML /-20.,160./      ! angle of light vector  +
     +   , XME,YME /-30., 15./      ! angle of the eye       + Integers?
     +   , SSC     /4./             ! image scale factor
     +   , XDS,YDS,ZDS /1.,1.,1./   ! anisotropic disp. scale factors
     +   , IREDRAW /1/              ! flag to redraw the title 'menu'?
     +   , NEL,NN,ITYPE /0,0,0/     ! initially no elements or nodes
     +   , C_R /0.,0.,0.,1.,1.,-1./ ! contour scale factors
     +   , ITEM2C  /1/              ! - token for a CV() entry :-)
     +   , IOP_ANIM/10*0/           ! animation token list
     +   , CCOL/14,1,14,1,-1,-1/    ! colours:use,from,to,offset/infill,wrap
     +   , ID_FACE/6*0/             ! facet# to draw (0='selective')
     +   , IMAGE/MIMAGES*-1/        ! saved images
     +   , IMAGEN/0/                ! # of auto-save images
c    +   , IOUT_D/1/                ! output destination (1=screen)
     +   , LD_CE/0/                 ! auto-glass:1=built-up,2=cut-out
C --------------- colour info vector -----------------------------------
      DATA (CV(I),I=1,10) /
     +   1    !  1.   face colour option 1=mats, 2=facet dirn, etc.
     +,  0    !  2.   actual face colour if CV(1) option = 0 ?Mat table!
     +, -1    !  3.   ()    ^ or 'secondary face color' (top layer)
     +,  1    !  4.   'standard' edge colour
     +, -1    !  5.    sub-edge colour
     +, -1    !  6.      ( material edge colour )
     +, -1    !  7.      ( geometry change colour )
     +, -1    !  8. 
     +, -1    !  9.    node colour
     +,  2 /  ! 10.    node size (pixels or 1/640ths. of image x ?)
      DATA (CV(I),I=11,20) /
     +  -1    ! 11.    node # colour
     +, -1    ! 12.    disp vector colour
     +, -1    ! 13.    element # colour
     +, -1    ! 14.      contour line type -1,0,1,(2=coloured)
     +,  2    ! 15.      contour fill type -1,0,1,(2=col,7/8=zebra)
     +, -1    ! 16.   undeformed facet colour
     +,  2    ! 17.   depth sort (1=nearest,2=furthest,else =central)
     +, 14    ! 18.  ** number of colours to use for contouring **
     +, -1    ! 19.  * # sub-facets per element (per side)
     +,  0 /  ! 20.  * screen background colour
      DATA (CV(I),I=21,30) /
     +  -1    ! 21.    graticule colour
     +,  1    ! 22.    # x-graticule lines
     +,  1    ! 23.    # y-graticule lines (?) 
     +, -1    ! 24.    axes colour
     +, -1    ! 25.       ( axes lengths )
     +, -1    ! 26.       ( flag fo -ve axes too )
     +, -1    ! 27.    axes label colour (x,y,z)
     +, -1    ! 28.    window border colour
     +, -1    ! 29.       ( redraw . stamp colour) ?? skipped?
     +, -1 /  ! 30.       ( load-step colour ) .. or have a 'status' box
      DATA (CV(I),I=31,40) /
     +   0    ! 31.    ! shrink factor (+ve=%-age,-ve = in pixels)
     +,  1    ! 32.    ! backface-poly cull 0=none,1=back,2=front
     +, -1    ! 33.    ! 
     +, -1    ! 34.    ! default menu = #1  (full screen)
     +, -1    ! 35.    ! (last menu ?)
     +, -1    ! 36.    ! 
     +, -1    ! 37.    ! multi-view type (-1=none,1=load#,2=view)
     +, -1    ! 38.    ! view number (-1 = via XME,YME)
     +, -1    ! 49.    ! 
     +, -1 /  ! 50.    ! 
      DATA (CV(I),I=41,50) /       ! windowing etc.
     +   1    ! 41.    ! picture #
     +,  1    ! 42.    ! # pics. in x-direction
     +,  1    ! 43.    ! # pics. in y-direction           
     +,  0    ! 44.    ! image size - x              = res(1)    ?
     +,  0    ! 45.    !            - y              = res(2)    ?
     +,  0    ! 46.    !            - offset in x    = res(4)    ?
     +,  0    ! 47.    !            - offset in y    = res(5)    ?
     +,  0    ! 48.    !%image centre offset - x  (ie. 'panning' )
     +,  0    ! 49.    !                     - y 
     +, -1 /  ! 50.    ! 

C--------------------- output options ----------------------------------
C.... should ALL be as ONE output dest. option ! 'cept bufferin'
      DATA IOP_BUF/0/, IROT/0/

C---------------------- initial constants ------------------------------
      IOUT_D = 1
      DTR = 3.14159265/180.
      ESC = CHAR(27)
      C_T =  esc//'[2;2H'     ! just needs to be a position on the screen 
      DO i=1,20
        IC_MATS(I) = I+1   ! set material-cols to be 'themselves'
      ENDDO
      PATH=curdir@()         ! default directory is the current one ?
c      i = leng(path)
c      path(I+1:i+5) ='\*.*'   ! append with a wildcard
C-------------------------- read in the menus --------------------------
      CALL READ_MENU (TABLE,TXT,MOPS,MENUP, 640,480)
      OPEN (19,file='debug')

C----------------------- try being in VGA allways ? --------------------
      CALL VGA@()

      IRESX = 640          ! these should also be set at the top !
      IRESY = 480
      IRESC = 16
      IRESCP = NINT (LOG(REAL(IRESC))/LOG(2.))

        RES(4) = 0
        RES(5) = 0
        RES(1) = IRESX
        RES(2) = IRESY
        IRESCP = NINT (LOG(REAL(IRESC))/LOG(2.))
        CALL SET_PAL ( 2, 0,0,0,CV(18)) ! set the default colour palette

C----------------------- reset mouse -----------------------------------
C... possibly go into a VGA mode first ??.......
      CALL INITIALISE_MOUSE@()        !
      CALL MOUSE_SOFT_RESET@(LMOUSE)  ! =.false. if mouse is not present
      CALL SET_MOUSE_POSITION@ (iresx-2,12)
      CALL SET_MOUSE_BOUNDS@ (1,1,IRESX-2,IRESY-2)
      CALL DISPLAY_MOUSE_CURSOR@()    !

C------------------- interpret the command line ------------------------
      CALL   CMNAMR()
      LINE = CMNAM@()

      DO WHILE (LINE.ne.' ')
        IF (LINE(1:1).NE.'-') THEN    ! read a filename (or wildcard)
          IF (line(2:2).eq.'\'.or.line(3:3).eq.':')then
            file_dat=line(2:leng(line)) ! already have a pathname
          else
            FILE_DAT = PATH(1:LENG(PATH)) //'\'// LINE(1:LENG(LINE))
          endif

        ELSEIF(LINE(1:2).eq.'-c') THEN  !  initial config file
          OPEN (11,FILE= LINE(3:leng(LINE)),IOSTAT=iostat)
          IF (IOSTAT.ne.0) THEN
            PRINT*,'***WARNING:  Config file not found'
          ELSE
            DO WHILE (IOSTAT.EQ.0)
              READ(11,*,IOSTAT=IOSTAT) I,CV(I)
            ENDDO
            CLOSE(11)
          ENDIF
        ELSEIF (LINE(1:2).eq.'-k'       ! command macro 
     +    .or.  LINE(1:2).eq.'-K' )THEN  

          IF(LINE(1:2).eq.'-K') THEN          ! get from a file
            CALL OPENRW@ (LINE(3:leng(LINE)),   IO, IERROR)
            CALL READF@  (LINE(3:3),IO,60L,NBYTES_R,IERROR)
            CALL CLOSEF@(IO,IERROR)
            II = NBYTES_R +2
          ELSE
            II = LENG(LINE) + 1
            LINE(II:II) = '~'
          ENDIF

          DO I=3,II
            J = ICHAR(LINE(I:I))
            IF (J.EQ.ICHAR('~')) J=13
            IF (J.EQ.ICHAR('_')) J=32
            CALL FEED_KEYBOARD@( J, IERROR)
            CALL DOSERR@ (IERROR)
          ENDDO
        ENDIF
        LINE = CMNAM@()  ! get the next command-line-arguement
      ENDDO


C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C        S T A R T    O F    T H E    M A I N    L O O P
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 1111 CONTINUE
      CALL SLEEP@(.1)

      IBUT1 = 0
      IBUT2 = 0
      KEY2 = -99
c     IF (LMOUSE) CALL GET_MOUSE_POSITION@(IXM,IYM,IBUTTON)
      IF (LMOUSE) CALL GET_MOUSE_BUTTON_PRESS_COUNT@(0,IBUT1)
      IF (LMOUSE) CALL GET_MOUSE_BUTTON_PRESS_COUNT@(1,IBUT2)

C---------- find which menu button has been pressed --------------------
      IF (IBUT1+IBUT2.gt.0) THEN
        CALL GET_MOUSE_POS@(IXM,IYM)
        CALL Q_MENU (IXM,IYM, TABLE,MENUP,CV(34), ITEM,KEY,KEY2,XM,YM)

C-------------- if none try the main menu (#2) -------------------------
C .. hmm.. I don't like this .. better to post the 'true' main menu
C .. if the mouse is neither in a menu or in a 'picture window'
        IF(KEY.eq.0)
     +  CALL Q_MENU (IXM,IYM, TABLE,MENUP,   2, ITEM,KEY,KEY2,XM,YM)

C----------------- Ok ... try to find the picture # --------------------
        IF (KEY.eq.0) THEN
          XM = REAL(IXM) / RES(1)
          YM = REAL(IYM) / RES(2)
          IF (XM.LT.1.OR.YM.LT.1) THEN   ! 'within' the drawing area
            KEY = ICHAR('W')
          ENDIF
        ENDIF

C------------------ beep if we are posting a menu ----------------------
C .. note if I 'miss' a button then 'beep' a warning ? !

        if (key.ne.0) then  !------- a 'hit'----------
          CALL HIDE_MOUSE_CURSOR@()
          CALL ISWAP ( TABLE(6,ITEM),TABLE(7,ITEM) )
          CALL POST_WIDGET (TABLE,TXT(ITEM),ITEM)
          CALL DISPLAY_MOUSE_CURSOR@()
          IF (BEEP) call sound@ (1024,1)

        else ! ----------------- a 'miss' ------------------------------
          call sound@ ( 256,2)
        endif

C--------------- wait for button release -------------------------------
 1003   CALL GET_MOUSE_POSITION@(IXM,IYM,IBUTTON) 
        IF (IBUTTON.ne.0) goto 1003
C------------------ raise button ---------------------------------------
        if (key.ne.0) then
          CALL HIDE_MOUSE_CURSOR@()
           CALL ISWAP ( TABLE(6,ITEM),TABLE(7,ITEM) )
           CALL POST_WIDGET (TABLE,TXT(ITEM),ITEM)
          CALL DISPLAY_MOUSE_CURSOR@()
        ENDIF
      ELSE  ! (--end of mouse-driven menus--)

C------------ no mouse button pressed so scan keyboard -----------------
        CALL GET_KEY1@(KEY)
         IF (NEL.eq.0) THEN
           KEY = ichar ('#') ! force data read if no elements 
           KEY2 = 1          ! ( as a 're-read' of the data ! )
         ENDIF

         IF (KEY.EQ.0.AND.IOP_ANIM(1).NE.0) KEY = 32  ! animation redraw
         IF (KEY.ne.0) goto 1001  ! goto event handler

C.. Ok post the title menu :-)  since we have nothing else to do
        IF (CV(34).gt.1.and.iredraw.eq.1) THEN ! I don't really like this
        DO I = MENUP(2), MENUP(2+1) -1
          CALL POST_WIDGET (TABLE,TXT(I),I)
        ENDDO
        iredraw =0
        ENDIF
        GOTO 1111
      ENDIF
C ---------------------------------------------------------------
 1001 continue
C... if key = 'space' turn into '_'
      IF (KEY.EQ.32) KEY = ICHAR('_')
      iredraw=1
C----------------------------------------------------------------------
C-----------------------  key event handling --------------------------
C----------------------------------------------------------------------
      CKEY = char(0)
      IF (KEY.GE.0.and.KEY.LE.255) CKEY = CHAR(KEY)

      IF (KEY.eq.0) THEN
C-----------------> draw image ----------------------------------------- 
      ELSE IF (CKEY.eq.'_') THEN
        GOTO 901   !      Jump straight to the image draw bit !
C -------- choose input file name -------------------------------------
      ELSEIF (CKEY.EQ.'#')THEN
C... key2 = 0 new data, =1 re-read = 2, append data etc.
       IF (KEY2.eq.0.or.KEY2.eq.-99) FILE_DAT ='-------'  ! assume no file 

        IOSTAT = -1
        DO WHILE (IOSTAT.NE.0)
          PRINT*,C_T//'FILE="',FILE_DAT(1:LENG(FILE_DAT)),'"'
          OPEN (U3,FILE=FILE_DAT,STATUS='OLD',iostat=iostat)
          IF (IOSTAT.NE.0) THEN
c            CALL HIDE_MOUSE_CURSOR@()
C            CALL TEXT_MODE@()

            i = leng(path)
            path(I+1:i+5) ='\*.*'   ! append with a wildcard

            CALL SET_VIDEO_DAC@(  1,5S,5S,20S) 
            CALL SELECT_FILE@ (PATH,file_dat,*999)
            CALL SET_ALL_DACS()
          ENDIF

c-------- now set PATH='the root name of the file's directory' ---------
            i = leng(file_dat)
            do while(FILE_DAT(I:I).ne.'\')
              i = i-1
            enddo
            idir_name = i-1
            PATH = FILE_DAT(1:idir_name)    ! so PATH = directory name
c-----------------
        ENDDO

C------------------ process ZIP files :-) (not yet !) ------------------
        IF (INDEX(FILE_DAT,'.ZIP').gt.0) THEN
          CALL CISSUE ('pkunzip -v '//FILE_DAT,ifail)
          print*,'ifail=', ifail
          print*,'which file ?'
          READ(*,'(a)') FILE2 
          CALL CISSUE ('pkunzip -o '//FILE_DAT//' \ '//FILE2,ifail)
          FILE_DAT = 'c:\'//FILE2
          OPEN (U3,FILE=FILE2,STATUS='OLD',iostat=ifail)
        ENDIF

C======================================================================
C ======== read basic data ============================================
C======================================================================
      call clock@(time(10))

C-----------------------------------------------------------------------
      PRINT*,'reading geometry.....'
C-----------------------------------------------------------------------

      CALL UPCASE (FILE_DAT)
                                       IDT = 0    ! default = 'unknown'
      IF (INDEX(FiLE_DAT,'.DIS').ne.0) IDT = 1   !- my .PL format(IAN's) 
      IF (INDEX(FiLE_DAT,'.MES').ne.0) IDT = 1   !- my .PL format(IAN's) 
      IF (INDEX(FILE_DAT,'.PL' ).ne.0) IDT = 1   !- my .PL format 
      IF (INDEX(FILE_DAT,'.PL2').ne.0) IDT = 1   !- my .PL + xtra lines
      IF (INDEX(FILE_DAT,'.G'  ).ne.0) IDT = 2   !- 'Objext File Format'
      IF (INDEX(FILE_DAT,'.GEO').ne.0) IDT = 2   !- 'Objext File Format'
      IF (INDEX(FILE_DAT,'.HO' ).ne.0) IDT = 3   !- David Ho's format 
      IF (INDEX(FILE_DAT,'.NFF').ne.0) IDT = 4   !- NFF format
      IF (IDT.eq.0) THEN
        PRINT*,'file=',FILE_DAT,'is of unknown type !'
      ENDIF
       
      IF (INDEX(FILE_DAT,'.PL2').ne.0) THEN
        PRINT*,'reading Extra lines.....'
        READ (U3,*) N_X_L , ((XTRA_L(J,I),J=1,4),I=1,N_X_L)
      ENDIF

      IF (IDT.eq.4) THEN
           CALL READ_NFF (U3,IDT,NN,NDIM,NDIME,ITYPE,
     +     GCOORD,3,INF,NEL,NUMS,INUMS,MEL,NLDS,CV,PAL)   ! 'NFF format'
      ELSE
c           read(u3,*)junk --?? !
      IF (IDT.eq.1.or.IDT.eq.2.or.IDT.eq.3)
     +                 CALL READ_NODES (U3,IDT,NN,NDIM,GCOORD,3,INF,NEL)
      NDIME = NDIM  ! default elems dimension = nodes' dimension
      IF (IDT.eq.1.or.IDT.eq.2.or.IDT.eq.3)
     +           CALL READ_ELEMS (U3,IDT,NEL,NDIME,ITYPE,NUMS,INUMS,MEL)
      NODOF = NDIM  ! default disps dimension = nodes' dimension
      IF (IDT.eq.1.or.IDT.eq.2.or.IDT.eq.3)
     +           CALL READ_LOADS (U3,IDT,NLDS,GDISPS,3,IGDISPS,NODOF,NN)
      ENDIF
      CLOSE(U3)

C-------- viewing parameters .. object position + view direction -------
      CALL RANGE (GCOORD,IGCRD ,NN,NDIM, COA,DIAG)
      IF(NDIM.EQ.2)THEN                 ! if 2d view straight 'on'
        XME = 0.     ! x-eye position
        YME = 0.     ! y-eye position
        CV(32) = 0   ! no backface_poly_cull
      ENDIF
C------- set a 'nice' default displacement scale factor ---------------
C         .. based on the 'final' load step..
      LD   = NLDS
      FACT = 0.
      IF (NLDS.GT.0) THEN  
        DISM  = 1.E-10
        IBASE = INTL(NLDS-1) * NN
        DO I=1,NN
          DISM1 = 0
          DO J=1,NODOF
            DISM1 = DISM1 + GDISPS(J,IBASE+I)**2
          ENDDO
          DISM = MAX(DISM,DISM1)
        ENDDO
        DISM = SQRT(DISM)
        FACT = DIAG / DISM / 25.   ! should be 'device-independant'!
      ENDIF

C--------- turn the NUMS into a list of external facets 'FACETS' -------
      call clock@(time(11))
      PRINT*,'Please wait as the FACETS are assembled'
      CALL FSTRIP2 (NEL,NUMS,INUMS,FACETS,PL2,IFCETS,NFCETS,NDIME)
c-----------------------------------------------------------------------
      call clock@(time(12))

C------------- set up default painter's list ---------------------------
C .... this maybe will dissappear for 2d (ie. clip but skip RSORT@() )
      IF (NDIM.EQ.2) THEN
        DO I=1,IFCETS
          PL1(I) = I
        ENDDO
      ENDIF
C======================== ENTRY TO GRAPHICS ============================
C... should really only need to return to graph mode
C... IRESX, etc. should be set to their default values at the top!
      CV(34) = 1  ! default menu = #1 (whole screen)
c     CALL VGA@
      CKEY='_'     ! draw the image 

C--------------------- screen mode -------------------------------------
      ELSEIF (CKEY.EQ.'/') THEN
        IRESX = -1
        IF (KEY2.gt.0) THEN
          LINE=' 320 200 256'//
     +         ' 640 480  16'//' 640 480 256'//
     +         ' 800 600  16'//' 800 600 256'//
     +         '1024 768  16'//'1024 768 256'
        READ(LINE,'(99I4)') (IRESx,IRESy,IRESC,i=1,KEY2)
        print*,iresx,iresy,iresc
        print*,'--- press <RETURN> ----'
        read*
        ENDIF

        IERR = -1
        DO WHILE (IERR.NE.0) 
          IF (IRESX.lt.0) then
            PRINT*, 'IRESX,IRESY,IRESC =? (640,480,16)'
            READ*,IRESX,IRESY,IRESC
          endif
          CALL GRAPHICS_MODE_SET@ (IRESX,IRESY,INTL(IRESC),IERR)
          DO I=1,15                    !
            CALL SET_PALETTE@ (I,I)    ! set the palette :-)
          ENDDO                        !
          CALL SET_ALL_DACS ()      !
        ENDDO
        IF (.not.(iresx.eq.640.and.iresy.eq.480.and.iresc.eq.16)) then
          cv(34)=-1   ! switch-off menues in SVGA !
        endif
        RES(1) = IRESX
        RES(2) = IRESY
        IRESCP = NINT (LOG(REAL(IRESC))/LOG(2.))
        PRINT*,IRESX,IRESY,IRESC, IRESCP
        CALL SET_MOUSE_BOUNDS@(1,1,IRESX-1,IRESY-1)
c       CALL READ_MENU (TABLE,TXT,MOPS,MENUP, iresx,iresy)
c       ckey='_'
C-------------------> clear darwing area -------------------------------
      ELSE IF (CKEY.eq.'!') THEN
        call clear_screen_area@ (0,0,res(1),res(2),0)
C-----------------------> beep on/off ----------------------------------
      ELSE IF (CKEY.eq.'B') THEN
        beep = .not.BEEP
C-----------------------> flash on/off ---------------------------------
      ELSE IF (CKEY.eq.'f') THEN
        flash = .not.flash
C-----------------------> no operation :-)  ----------------------------
      ELSE IF (CKEY.eq.'Z') THEN
C------------------------> post a menu ---------------------------------
      ELSE IF (CKEY.eq.'m') THEN
        IF (BEEP) call sound@ (1024,1)    ! 'beep' if desired !
        IF (BEEP) call sound@ (1194,1)
        IF (KEY2.LT.0) THEN
          PRINT*,C_t//'Which menu number ?'
          READ*,KEY2
        ENDIF
        CV(34) = KEY2 
        CALL HIDE_MOUSE_CURSOR@()
          IF (.not.(iresx.eq.640.and.iresc.eq.16))then
            CALL VGA@()
            DO I=1,15                    !
              CALL SET_PALETTE@ (I,I)    ! set the palette :-)
            ENDDO                        !
            CALL SET_ALL_DACS ()      !
            IRESX = 640     !--- force back into VGA mode if menu!
            IRESY = 480
            IRESC  = 16
            IRESCP =  4

            RES(1) = IRESX
            RES(2) = IRESY
            PRINT*,IRESX,IRESY,IRESC, IRESCP
            CALL SET_MOUSE_BOUNDS@(1,1,IRESX-1,IRESY-1)
          ENDIF

        DO I = MENUP(CV(34)), MENUP(CV(34)+1) -1
          CALL POST_WIDGET (TABLE,TXT(I),I)
        ENDDO
        CALL DISPLAY_MOUSE_CURSOR@()
c.. show mouse if applicable
        IF (CV(34).EQ.1) THEN
          CALL HIDE_MOUSE_CURSOR@()
        ELSE
          CALL DISPLAY_MOUSE_CURSOR@()
        ENDIF
c.......... force 'window' mode if currnetly too big to fit ............
        IF (RES(1).gt.479) THEN
            RES(1) =  479
            RES(2) =  479
        ENDIF
        IF (KEY2.eq.1) THEN
          RES(1) = IRESX
          RES(2) = IRESY
        ENDIF

c---------------------------> exit -------------------------------------
      ELSE IF(CKEY.EQ.'q') THEN
        CALL TEXT_MODE@
        STOP '--STOP--'
      ELSE IF(KEY.EQ.27) THEN
        CALL TEXT_MODE@
        STOP '--escape--'
C------------------------> dos shell -----------------------------------
      ELSE IF(CKEY.EQ.'d') THEN
        PAUSE '--PAUSE--'

C------------------------> animation -----------------------------------
      ELSEIF (CKEY.EQ.'A') THEN

        IF (KEY2.ge.0) THEN      ! mouse driven preset options ?
c         IOP(ANIM(KEY2)) = ??
        ELSE 
          PRINT*,'C_T'//'Enter the <animation parameters>'
          PRINT*,'eg .. 1,0,  10,5,   0,-5,  1,  1,32'
          PRINT*,'="on","no disp", eye x+10,y+5,light x+0,y+5'//
     +           '--> to 320PCX as 1:32 pics'
          READ*, (IOP_ANIM(i),i=1,9)
          IF (IOP_ANIM(7).lt.10)THEN
            print*,'enter the animation root name (eg. FRED)'
            read(*,'(a)') FILE_ROOT
          endif
          ipic_n = iop_anim(8)   !-- the first frame number ??
          IF (iop_anim(2).ne.0) THEN                
c           LD  = 0  ! start at zero ?
            LDD = REAL(NLDS) / REAL(iop_anim(9)-1)
          ELSE
            LDD = 0
          ENDIF
        ENDIF
C------------- 'facet# to always/never draw' ---------------------------
      ELSEIF (CKEY.EQ.'|')THEN
        write(*,'(a,a,6i2)')C_T,'facet codes(1..6)?=',(ID_FACE(i),i=1,6)
        read*,(ID_FACE(i),i=1,6)

C------------- backface polygon cull -----------------------------------
      ELSEIF (CKEY.EQ.'\')THEN
        IF (KEY2.ge.0) THEN
          IF (KEY2.LE.3) CV(32) = KEY2            ! b_p_c
          IF (KEY2.LE.6) CV(17) = KEY2-3          ! depth sort
        ELSE 
          CV(32) = MOD(CV(32)+1,3) ! just toggle if keyb.
        ENDIF 
C------------- facet sub-division --------------------------------------
      ELSEIF (CKEY.EQ.'@')THEN
        PRINT*,C_T,'How many sub-facets per facet (',CV(19),' ) ?'
        READ*,CV(19)
        CKEY = '_'
C------------- input macro ---------------------------------------------
      ELSEIF (CKEY.eq.'k') THEN    ! (rarely used !)
        PRINT*,C_T//'Enter macro string (_=space,~=CR)'
        READ(*,'(a)') LINE
        DO I=1,leng (LINE)
          IF (LINE(I:I).eq.'~') THEN
            CALL FEED_KEYBOARD@(13, ierror)
          ELSEIF (LINE(I:I).eq.'_') THEN
            CALL FEED_KEYBOARD@(32, ierror)
          ELSE
            CALL FEED_KEYBOARD@(ICHAR(LINE(I:I)), ierror)
          ENDIF
        ENDDO
C----------------------- colors, etc (+node draw) ----------------------
      ELSEIF (CKEY.EQ.'I') THEN
       IF (KEY2.eq.-99) then
         print*,C_T,'Item Number ?'
         read*,key2
       ENDIF           ! pass ITEM2C as a token to 'C' option :-)
       ITEM2C  = KEY2
C------- use 'J' to pointer to "op-codes" not '-ve cols'! --------------
      ELSEIF (CKEY.EQ.'J') THEN
       IF (KEY2.eq.-99) then
         print*,C_T,'Item Number ?'
         read*,key2
       ENDIF           ! pass ITEM2C as a token to 'C' option :-)
       ITEM2C  = -KEY2
C------------------> set an item's color -------------------------------
C           (or take other appropriate action)
      ELSEIF (CKEY.EQ.'C') THEN
       IF (KEY2.eq.-99) then     
         print*,C_T,'colour value ?'
         read*,key2
       ENDIF
       IF (ITEM2C.GT.0) THEN
         CV(ITEM2C) = KEY2     ! just set the CV() vector

C ....> special cases of ITEM2C .. are handled individually
C-------------------> glassing materials -------------------------------
        ELSEIF (ITEM2C.EQ.-2) then
          DO I=1,NEL
            NOD  =          NUMS(I,1)
            IMAT =         (NUMS(I,NOD+2))
            IF(ABS(IMAT).eq.KEY2) NUMS(I,NOD+2) = -IMAT  !-- flip
          ENDDO
C----------> shade palette based on an existing colour -----------------
        ELSEIF (ITEM2C.EQ.-3) then
          IR = PAL(1,KEY2)
          IG = PAL(2,KEY2)
          IB = PAL(3,KEY2)
          CALL SET_PAL ( 1,ir,ig,ib,CV(18))
        ELSEIF (ITEM2C.EQ.-4) then              !-----redefine a colour
          print*,C_T,'enter the new rgb values'
     +              ,PAL(1,KEY2),PAL(2,KEY2),PAL(3,KEY2)
          read*,ir,ig,ib
          CALL SET_PAL (-KEY2,ir,ig,ib,CV(18))
C----------- overscan colour ----------------
        ELSEIF (ITEM2C.EQ.-5) then
          CALL SET_PALETTE@ (17,KEY2)
C----------- load step # patch --------------
        ELSEIF (ITEM2C.EQ.-6) then
          LD = KEY2
C----------- save image to memory --------------
        ELSEIF (ITEM2C.EQ.-7) then
         CALL GET_SCREEN_BLOCK@ (0S,0S,RES(1),RES(2),IMAGE(KEY2))
C----------- load iamge to memory --------------
        ELSEIF (ITEM2C.EQ.-8) then
        IF(IMAGE(KEY2).ne.-1L) THEN
          CALL DOSERR@(IFAIL)
        ELSE
          CALL SOUND@( 256,2)
          print*,'****ERROR: IMAGE',KEY2,' has not been saved' 
        ENDIF
        CALL DOSERR@(IFAIL)
C----------- auto save image to memory --------------
        ELSEIF (ITEM2C.EQ.-9) then
        IMAGEN = IMAGEN+1
         CALL GET_SCREEN_BLOCK@ (0S,0S,RES(1),RES(2),IMAGE(IMAGEN))
        ENDIF
C------------------------ shade palette --------------------------------
C --> shade range specification
      ELSEIF (CKEY.EQ.';') THEN
        IF (KEY2.eq.1) THEN                       ! blue --> yellow
           line = '2 15 1.     0 0 .7   1 1 0 '
           call set_shade (line,pal)
        ELSEIF (KEY2.eq.2) THEN                   ! green--> red
           line = '2 15 1.     0 .5 0   1 0 0 '
           call set_shade (line,pal)
        ELSEIF (KEY2.eq.3) THEN                   ! blue->black->red
           line = '2  8 1.     0 1 1    0 0 0 '
           call set_shade (line,pal)
           line = '8 15 1.     0 0 0    1 0 0 '
           call set_shade (line,pal)
        ELSE
           line = '2 15 1.     0 0 .7    1 1 0 '
           call read_edited_line@ (line,0,1,15,ifail)
           call set_shade (line,pal)
        ENDIF
C --> individual color specification  <- used for pre-def palettes etc.
      ELSEIF (CKEY.EQ.'.') THEN
        IF (KEY2.ne.-99)THEN
          ICOL = KEY2 
        ELSE
          PRINT*,C_T//'enter code (0-15,-1 for all) & rgb values(0-255)'
          PRINT*,    '--> eg. 0,10,0,0 for a red background'
          read*,icol,ir,ig,ib
        ENDIF
        call set_pal (icol,ir,ig,ib,CV(18))

C------------------- palette cycling -----------------------------------
      ELSEIF(CKEY.EQ.'P') THEN
        DO I=1,5*3
          CALL SET_PAL (10,IR,IG,IB,CV(18))
        ENDDO

C----------------> 'hand' face color selection ------------------------
      ELSEIF (CKEY.EQ.'c') THEN
        IF (KEY2.ne.-99) THEN      ! menu pointer via 'co-ord' menu
          CV(1) = KEY2
        ELSE
          PRINT*,C_T//'enter facet color option (',CV(1),')'
          read*,CV(1) 
        ENDIF
        CKEY='_'

C------------------- contour range setting -----------------------------
      ELSEIF (CKEY.EQ.'^') THEN
        IF (KEY2.eq.1) THEN    !-- auto set range
          C_R(3)= c_r(1)
          C_R(4)= c_r(2)
          C_R(5)= 1.
          C_R(6)=-1.
        ELSE      !--- manualy set the contour range
c... we need to auto-set both the 'max' AND the 'min' contour values
c...  -and be able to type these limits in
C       print*, C_T//'(dis2=',dis2,' dis1=',dis1,') disp=', dis_p,' =?'
        print*,C_T
        print*,'# contours =', CV(18)
        print*,      'actual  range =',C_R(1),'<---->',C_R(2)
        print*,      'current range =',C_R(3),'<---->',C_R(4),
     +               ' power=',C_R(5)
        print*,'enter new # contours, min, max, & power, (-ve for abs)'
C.. use read_edited line ?
        read*, CV(18), C_R(3), C_R(4), C_R(5)
        C_R(6) = SIGN (1.,C_R(5))
        C_R(5) =  ABS    (C_R(5))
        ENDIF
c       CKEY='_'

C---------- contour method (=Gouraud/lines or averaged) ----------------
      ELSEIF (CKEY.EQ.'g') THEN
        IF (KEY2.ge.0.and.KEY2.le.9) THEN    !-- cont. method ?
          CV(18) = key2
        ELSEIF (KEY2.eq.10) THEN     !-- # colours
          PRINT*,'# contours?=',CV(18) 
          READ*, CV(18) 
c         READ*, CCOL(1)
        ELSEIF (KEY2.eq.11) THEN     !-- color range + offset
          PRINT*,'palette from/to &offset?=',CCOL(2),CCOL(3),CCOL(4)
          READ*,CCOL(2),ccol(3),ccol(4)

c        ELSEIF (KEY2.eq.12) THEN     !-- 'infill' type
c          CCOL(5) = -1                  ! no zebra infills
c        ELSEIF (KEY2.eq.13) THEN     
c          CCOL(5) =  0                  ! alternate black
c        ELSEIF (KEY2.eq.14) THEN     
c          CCOL(5) =  1                  ! alternate see-thru
c
c        ELSEIF (KEY2.eq.15) THEN     !-- 'outside' range colour  ___
c          CCOL(6) = -1                  ! = max/min          ___/
c        ELSEIF (KEY2.eq.16) THEN        !                        . .
c          CCOL(6) =  0                  ! = black            . ./
c        ELSEIF (KEY2.eq.17) THEN        !                        . .
c          CCOL(6) =  1                  ! = see-thru         . ./ 
c        ELSEIF (KEY2.eq.17) THEN     
c          CCOL(6) =  2                  ! = cycle round =MOD() / / / 
c        ELSEIF (KEY2.eq.18) THEN     
c          CCOL(6) =  3                  ! ='bounce'round       \ / \
         ENDIF

C--------- screen image store/restore ----------------------------------
       ELSEIF (CKEY.eq.'}') THEN    !------ forwards
         DO j=1,999
         DO I=iop_anim(8),iop_anim(9)
          CALL RESTORE_SCREEN_BLOCK@    (0S,0S,IMAGE(I),0,IFAIL)
          IF (KEY_WAITING@()) GOTO 1111
         ENDDO
         enddo
       ELSEIF (CKEY.eq.'{') THEN   !---- forwards and backwards
         do j=1,999
         DO I=iop_anim(8),iop_anim(9)-1
          CALL RESTORE_SCREEN_BLOCK@    (0S,0S,IMAGE(I),0,IFAIL)
          IF (KEY_WAITING@()) GOTO 1111
         ENDDO                                        
         DO I=iop_anim(9),iop_anim(8)-1,-1
          CALL RESTORE_SCREEN_BLOCK@    (0S,0S,IMAGE(I),0,IFAIL)
          IF (KEY_WAITING@()) GOTO 1111
         ENDDO
        enddo

C------------------> glassing materials --------------------------------
        ELSEIF (CKEY.EQ.'[') then    !------ built-up models------
          LD_CE = 1
        ELSEIF (CKEY.EQ.']') then   !------ excavation models-------
          LD_CE = 2

C----------------- output destination choice ---------------------------
C ........ should really be a typed numeric option
C ---> PCX screen dump .. 'just the current drawing area ?'

      ELSEIF (CKEY.EQ.'O')THEN
       IF (KEY2.eq.-99) then
         print*,C_T,'Output destination ?'
         print*,'(1=screen,3=PCX,4=laser,5=PS,6=pro-p,8=HP7550)'
         read*,key2
       ENDIF       

      IF (key2.eq.1) THEN               !---screen----
        IOUT_D = KEY2
        print*,'image size (640,480) ?'
        read*,RES(1),RES(2)
        AR = -1.   !.. but NOT=-1.if 320x200 mode etc.

c--- this should not realy be here !!!
      ELSEIF (key2.eq.2) THEN           !---(PCX screen dump) --
        CALL GET_SCREEN_BLOCK@(0,0,639,479,BUFFER)
        PRINT*,C_T//'PCX screen dump File name: ?'    
        READ(*,'(A)')FILE_PCX
        CALL SCREEN_BLOCK_TO_PCX@(FILE_PCX,BUFFER,IFAIL)
        CALL RETURN_STORAGE@(BUFFER)
        CALL PCX_FIX (FILE_PCX)

      ELSEIF (KEY2.EQ.3)THEN             ! redraw to a PCX
        IOUT_D = KEY2
        PRINT*,'X,Y,(D) RESOLUTION ?  (1200,1200,4)'
        READ*,RES(1),RES(2),IBITPL
        PRINT*,C_T//'PCX File name: ?'
        READ(*,'(A)')FILE_PCX

        CALL CREATE_SCREEN_BLOCK@(RES(1),RES(2),IBITPL,VSCREEN)
        CALL OPEN_VSCREEN@(VSCREEN,IFAIL)
        AR = -1.

      ELSEIF (KEY2.EQ.4) THEN    ! ---> laser-printer output
        IOUT_D = KEY2
        FILE_PCX='LASER'
        PRINT*,'printer resolution ? (75,100,150,300)'
        READ*,IPRES
        call select_pcl_printer@ (0,'A4',IPRES,RES(1),RES(2))
        CALL OPEN_GPRINT_FILE@     (FILE_PCX,ierror)
        AR = -1.

      ELSEIF (KEY2.EQ.5) THEN    ! ---> Postscript output
        IOUT_D = KEY2
        PRINT*,C_T//'PostScript File name: ?'
        READ(*,'(A)')FILE_PCX
        CALL SELECT_PS_PRINTER (FILE_PCX,RES(1),RES(2))
        AR =  1.

      ELSEIF (KEY2.EQ.6) THEN    ! ---> pro-printer output
        IOUT_D = KEY2
        CALL OPEN_GPRINT_DEVICE@     (1,IFAIL)
        RES(1) = 960
        RES(2) = 576
        AR = -1.3

        CALL DYNT@('SELECT_DOT_MATRIX@',ADDRESS)  ! does it exist ?
        IF (ADDRESS.gt.0) THEN
           CALL SELECT_DOT_MATRIX@ (0S,960S,576S)
        ENDIF

      ELSEIF (KEY2.eq.8) THEN         ! ---> HP7550 plot output
        IOUT_D = KEY2
        PRINT*,C_T//'HP7550 File name:'
        READ(*,'(A)') FILE_PCX
        CALL OPEN_PLOT_FILE@ (FILE_PCX,IIFAIL)
        RES(1) = 10870                        !-store the sizes in a 
        RES(2) =  7600                        ! table? so can modify
        AR = 1.

        ENDIF  ! end of ouput devices
        CKEY='_'

C-------------> set the current window number --------------------------
      ELSEIF (CKEY.EQ.'W') THEN  ! (from the mouse position)
        cv(41) = INT(XM*CV(42)) + INT(YM*CV(43)) * cv(42) + 1

C-------------> screen window size and position ------------------------
      ELSEIF (CKEY.EQ.'w') THEN
        PRINT*,C_T//'enter window size, x,y (460,460)'
        read*,res(1),res(2)
        res(4) = (479 -res(1)) / 2 
        res(5) = (479 -res(2)) / 2 
C-------------------> buffered screen output ---------------------------
      ELSEIF (CKEY.EQ.'b') THEN
        IOP_BUF = MOD(IOP_BUF+1,2)

C ------------------------ s c a l i n g -------------------------------
C --> menu-driven  eye-angle/eye-pan/light-source chenaging

      ELSEIF (CKEY.EQ.'V')THEN
        IF (KEY2.eq. 1) YME = YME - 5.     ! rotate the eye
        IF (KEY2.eq. 2) XME = XME + 5.     ! by 5 degress
        IF (KEY2.eq. 4) XME = XME - 5.
        IF (KEY2.eq. 5) YME = YME + 5.
        IF (KEY2.eq. 3) THEN
          PRINT*,C_T//'Enter the EYE longitude and latitude'
          PRINT*,'XME,YME=',xme,yme
          READ*,XME,YME
        ELSEIF (KEY2.ge.6.and.KEY2.le.9) THEN
           cv(38) = KEY2 - 5         ! a 'preset' code for CAMERA
        ENDIF

        IF (KEY2.eq.11) cv(47) = CV(47) + 10.  ! shift the image
        IF (KEY2.eq.12) CV(46) = CV(46) - 10.  ! by 10% screen width
        IF (KEY2.eq.14) CV(46) = CV(46) + 10.
        IF (KEY2.eq.15) CV(47) = CV(47) - 10.
        IF (KEY2.eq.13) THEN
          PRINT*,C_T//'Enter the PAN offset (0,0?)'
          PRINT*,'X_offset,Y_offset=',cv(46),cv(47)
          READ*,CV(46),CV(47)
        ENDIF       

        IF (KEY2.eq.21) YML = YML + 10.   ! rotate the light
        IF (KEY2.eq.22) XML = XML - 10.   ! by 10 degrees
        IF (KEY2.eq.24) XML = XML + 10.  
        IF (KEY2.eq.25) YML = YML - 10.
        IF (KEY2.eq.23) THEN
          PRINT*,C_T//'Enter the LIGHT longitude and latitude'
          PRINT*,'XML,YML=',xmL,ymL
          READ*,XML,YML
        ENDIF

        IF (KEY2.eq.31) THEN     ! change the perspective !
C         FEYE = ATAN (TAN(FEYE)
          DEYE = DEYE / SQRT(2.)   ! wide-angle
          SSC  = SSC / SQRT(2.)
        ELSEIF (KEY2.EQ.32) THEN  
          DEYE = DEYE * SQRT(2.)   ! narrow-angle
          SSC  = SSC * SQRT(2.)
        ENDIF

C ---> typed eye-position
      ELSEIF (CKEY.EQ.'t')THEN 
        PRINT*,C_T//'Enter the EYE longitude and latitude'
        PRINT*,'XME,YME=',xme,yme
        READ*,XME,YME
C ---> cursor key eye-position move
      ELSEIF (KEY.eq.331) THEN
        XME = XME + 5.
      ELSEIF (KEY.eq.333) THEN
        XME = XME - 5.
      ELSEIF (KEY.eq.328) THEN
        YME = YME - 5.
      ELSEIF (KEY.eq.336) THEN
        YME = YME + 5.
C ---> image scaling
      ELSEIF (CKEY.EQ.'+') THEN 
         SSC = SSC * SQRT(2.)
      ELSEIF (CKEY.EQ.'-') THEN  
         SSC = SSC / SQRT(2.)
      ELSEIF (CKEY.EQ.'*') THEN  
        WRITE(*,'(A,E8.2,A)') 'image scale ?(',SSC,' )'
        READ*,SSC
C ---> displacemnt scaling
      ELSEIF (CKEY.EQ.'<') THEN
        FACT = FACT / SQRT(2.)
      ELSEIF (CKEY.EQ.'>') THEN 
        FACT = FACT * SQRT(2.)
      ELSEIF  (CKEY.EQ.'s') THEN
        WRITE(*,'(A,E8.2,A)') 'displacement scale ?(',FACT,' )'
        READ*,FACT
C------------------> typed anisotropic displacemnt scaling
      ELSEIF (CKEY.EQ.'d') THEN
        WRITE(*,'(A,3F8.2,A)') '(',XDS,YDS,ZDS,
     +                             ' ) new x,y,z disp. scale ?'
        READ*,XDS,YDS,ZDS
C ---> typed centre -of-attention
      ELSEIF (CKEY.EQ.'œ') THEN  
        WRITE(*,'(A,3F10.4)') 'COA ?',COA(1), COA(2), COA(3)
        READ*,COA(1), COA(2), COA(3)
        CKEY = '_'
C ---> (screen x,y flipping)
      ELSEIF (CKEY.EQ.'%') THEN 
        IROT=1-IROT

C-------------------- 'Mirroring the model' ----------------------------
C---> push the work to a subroutine to save space ? :-)
C... also allow the option to 'remove' a block to 'un-mirror' 

      ELSEIF (CKEY.EQ.'R') THEN
        IF (KEY2.le.0) THEN
          PRINT*, 'Enter the direction (1=x,2=y,3=z) and the coord (0.)'
          READ*,IDIR,VAL
        ELSE
          IDIR = KEY2
           II=-1               ! default: Mirror wrt. the min coord
          IF (idir.ge.4.and.IDIR.le.6)then
           idir = idir-3       
           II=1                ! Mirror wrt. the MAX coord
          endif
        ENDIF
c ... if direction = 0 then *remove* elements !
c.. at the moment I am merely reducing the # facets.. should also
c  'trim' the element (and node &disps) lists !

        IF (IDIR.EQ.7) THEN     !  un-mirror.. = 'halve the model'
          NFCETS = NFCETS *1./2.
        ELSEIF (IDIR.EQ.8) THEN !  'quarter the model' eg. (spud-can)
          NFCETS = NFCETS *3./4.
        ELSEIF (IDIR.EQ.9) THEN !  'eighth the model'  ??
          NFCETS = NFCETS *7./8.

        ELSEIF (IDIR.GE.1.AND.IDIR.LE.3) THEN
        IF (II.lt.0) VAL = GCOORD(IDIR,NN+1)  ! min x,y,z is stored at NN+1
        IF (II.gt.0) VAL = GCOORD(IDIR,NN+2)  ! max x,y,z is stored at NN+2
        COA(IDIR) = VAL    ! move COA to the 'mirror' plane

        DO I=1,NN           ! copy the nodes
          DO J=1,NDIM
            IF (J.NE.IDIR) THEN
              GCOORD (J,NN+I) = GCOORD(J,I)        ! copy
            ELSE
              GCOORD (J,NN+I) = 2*VAL - GCOORD(J,I)  ! reflect 
            ENDIF
          ENDDO     
        ENDDO
        CALL RANGE (GCOORD,IGCRD, 2*NN,NDIM, COA,DIAG)  ! get DIAG,etc.

C--------------------- copy the displacement info ----------------------
C ----> assume the NN is constant, so copy from the top down 
 
        DO K=NLDS-1,0,-1     ! loop load-steps
          IB1 = INTL(K) * NN * 2
          IB2 = INTL(K) * NN
          DO I=NN,1,-1           ! loop nodes
            DO J=1,3                 ! loop freedoms
                GDISPS (J,IB1 +I   )  =  GDISPS(J,IB2 +I)  ! move + 
              IF (J.NE.IDIR) THEN
                GDISPS (J,IB1 +I+NN)  =  GDISPS(J,IB2 +I)  ! copy
              ELSE                         
                GDISPS (J,IB1 +I+NN ) = -GDISPS(J,IB2 +I)  ! reflect
              ENDIF
            ENDDO     
          ENDDO
        ENDDO

        DO I=1,NEL                 !----- copy the elements
          NOD  = NUMS(I,1)
          NUMS(NEL+I,1) = NOD                 ! # nodes
          DO J=2,NOD+1
            NUMS(NEL+I,J) = NUMS(I,J) + NN    ! node number
          ENDDO
          NUMS(NEL+I,NOD+2) = NUMS(I,NOD+2)   ! material type
        ENDDO

        CALL RANGE (GCOORD,IGCRD, 2*NN,NDIM, COA,DIAG)  ! get DIAG,etc.

        DO I=1,NFCETS                            !-----copy facets-----
          FACETS (1,NFCETS+I) = FACETS(1,I)+NEL
          FACETS (2,NFCETS+I) =-FACETS(2,I)      ! ie. reflect them too!
          FACETS (3,NFCETS+I) = FACETS(3,I)+NFCETS   !'touching pointer'
          IF (FACETS(3,I).eq.0) FACETS(3,I+NFCETS) = 0  ! fix 'boundary'
        ENDDO

        NN     = NN  * 2
        NEL    = NEL * 2
        NFCETS = NFCETS * 2
        ENDIF  !--- block for IDIR =1,=2, or =3

C----------------- 'linkage of load step # to picture #-----------------
      ELSEIF (CKEY.EQ.'F') THEN 
        IF (KEY2.eq.1) THEN         !  link load step #
          IF (CV(37).eq.1) THEN
            CV(37) = -1         ! switch-off
          ELSE
            CV(37) = 1          ! switch-on 
          ENDIF
        ELSEIF (KEY2.eq.2) THEN     !  link view direction
          IF (CV(37).eq.2) THEN
            CV(37) = -1         ! switch-off
            CV(38) = -1         ! view = non-auto (ie. can rotate)
          ELSE
            CV(37) = 2          ! switch-on 
            CV(38) = 1          ! first view = down x ? 
          ENDIF
        ENDIF
C-------------------- sum/difference load steps-------------------------
      ELSEIF (CKEY.EQ.'L') THEN

        IF (KEY2.eq.1) then   !------- summation ----------
          DO K=2,NLDS              ! loop load-steps (miss first :-)
            IB1 = INTL(K-2) * NN    
            IB2 = INTL(K-1) * NN
            DO I=1,NN               ! loop nodes
              DO J=1,3                 ! loop freedoms
                GDISPS (J,IB2+I) = GDISPS(J,IB2+I)+GDISPS(J,IB1+I)  
              ENDDO     
            ENDDO
          ENDDO

        ELSEIF (KEY2.eq.2) then  !------ diference ---------
          DO K=NLDS,2,-1              ! loop load-steps (miss first)
            IB1 = INTL(K-2) * NN 
            IB2 = INTL(K-1) * NN
            DO I=1,NN               ! loop nodes
              DO J=1,3                 ! loop freedoms
                GDISPS (J,IB2+I) = GDISPS(J,IB2+I)-GDISPS(J,IB1+I)  
              ENDDO     
            ENDDO
          ENDDO
        ENDIF

C-------------------- load step number ---------------------------------
      ELSEIF (CKEY.EQ.'l') THEN
        IF (KEY2.eq.1)then
          LD = MAX (0.,LD-1.)
        ELSEIF (KEY2.eq.2)then
          LD = MIN (REAL(NLDS),LD+1.)
        ELSE
          PRINT*,C_T//'enter load step 1.-->',NLDS,'(',LD,')'
          READ*,LD
        ENDIF
C------------------- element shrinking ---------------------------------
      ELSEIF (CKEY.EQ.'S') THEN
        print*,C_T//'ishrink=', cv(31),' ?'
        read*,cv(31)
C------------------- x,y,z up-direction toppling -----------------------
C      ( for Ian Williams et al. who think lop-sided :-)
      ELSEIF (CKEY.EQ.'a') THEN
C..   better as a subroutine to save space ??
        T      = COA(1)
        COA(1) = COA(2)
        COA(2) = COA(3)
        COA(3) = T
        DO I=1,NN  + 2   ! ie. all the nodes AND the x,y,z min/max
                   T  = GCOORD(1,I)  
          GCOORD(1,I) = GCOORD(2,I)
          GCOORD(2,I) = GCOORD(3,I)
          GCOORD(3,I) = T
        ENDDO
        DO IBASE = 0,INTL(NLDS)*NN-1,NN 
          DO I = 1,NN
                           T  = GDISPS(1,IBASE+I)  
            GDISPS(1,IBASE+I) = GDISPS(2,IBASE+I)
            GDISPS(2,IBASE+I) = GDISPS(3,IBASE+I)
            GDISPS(3,IBASE+I) = T
          ENDDO
        ENDDO
C-----------------> else invalid key-press detected <-------------------
      ELSEIF (KEY.NE.0) THEN
        PRINT*,C_T//'*** INVALID key stroke, code= -',KEY,'-',ckey,'-'
        CKEY=' '
      ENDIF
C-----------------> if CKEY='_' then we generate a refresh  
      IF ( IBUT2.ne.0) CKEY='_'                    ! redraw
      IF ( IBUT1.ne.0.or.KEY_WAITING@()) CKEY=' '  ! don't redraw
      IF ( .not.(CKEY.eq.'_'.or.KEY.eq.13))  GOTO 1111   ! go back

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ---- OK either we have moved or a new load case .. so draw -----------
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

  901 CONTINUE
      call clock@(time(1))

C----------------- create transformation matrix ------------------------
      CALL XM2EYE (XME,YME,EYE)      ! convert to a unit vector
      CALL XM2EYE (XML,YML,LIGHT)    !    "    "  "  "     "
       
      T = SIGN (1.,ABS(YME)-90.001)     ! 'fly -over-the-top'
      DO I=1,3
        UP(I) = SIGN( UP(I),T )
      ENDDO
      CALL CAMERA (EYE,COA,UP, DEYE*DIAG,FEYE, PM,cv(38))  ! create PM()

c------- set the aspect ratio for different ouput devices --------------

C............. set up the windowing parameters .........................
      iwinx =  mod(cv(41)-1  , cv(42) )
      iwiny =     (cv(41)-1) / cv(42)
      iwin_xw = res(1)/ cv(42)
      iwin_yw = res(2)/ cv(43)
      iwin_xo = iwin_xw * iwinx  + res(4)  !.. the picture offsets
      iwin_yo = iwin_yw * iwiny  + res(5)  !..

      SC_X = SSC * iwin_xw         ! x,y,z scaling factors (screen)
      SC_Y = SSC * iwin_xw * ar
      SC_Z = SSC * iwin_xw
      XCEN = iwin_xo + iwin_xw/2. + iwin_xw/2.*CV(46)/100. ! image-centre
      YCEN = iwin_yo + iwin_yw/2. + iwin_xw/2.*CV(47)/100.
      IF (IOP_BUF.eq.1) then   
        XCEN = XCEN - iwin_xo  ! if buffered then remove the offset
        YCEN = YCEN - iwin_yo  ! (then add when copying to the screen)
      ENDIF

C-------------- set up the factors of each loadstep to add -------------

      ILD2 = MAX(1,(MIN(INT(LD+1),NLDS)))
      ILD1 = ILD2 - 1
      FAC2 = LD - REAL(ILD1)
      FAC1 =  1. - FAC2
      IF (ILD1.EQ.0) THEN
        ILD1 = 1
        FAC1 = 0.
      ENDIF
      IB1  = INTL(NN) * (ILD1-1)
      IB2  = INTL(NN) * (ILD2-1)
      XDS1 = XDS * FACT * FAC1
      YDS1 = YDS * FACT * FAC1
      ZDS1 = ZDS * FACT * FAC1
      XDS2 = XDS * FACT * FAC2
      YDS2 = YDS * FACT * FAC2
      ZDS2 = ZDS * FACT * FAC2

C------------ auto-glassing of built-up/excavated models ---------------

      IF (LD_CE.eq.1) THEN  !--built-up models  (ie. mat. 'appears')
        DO I=1,NEL
          NOD  =          NUMS(I,1)
          IMAT =     ABS (NUMS(I,NOD+2))
          IF (IMAT.NE.1.and.IMAT.GT.LD+1.5) IMAT = -IMAT  !<-- mend this!
          NUMS(I,NOD+2) = IMAT  
        ENDDO
      ELSEIF (LD_CE.eq.2) THEN  !--excavated models  (mat. 'disappears')
        DO I=1,NEL
          NOD  =          NUMS(I,1)
          IMAT =     ABS (NUMS(I,NOD+2))
          IF (IMAT.NE.1.and.IMAT.LT.LD+0.5+1.01) IMAT = -IMAT
          NUMS(I,NOD+2) = IMAT
        ENDDO
      ENDIF

C------------ loop the facets to get the 'draw-list' -------------------
c .. I could mark the 'touch' column as -ve if the adjoining material 
c    is of a different material type?--> so can 'always' draw?

      NFCETS_D = 0
      DO IFC=1,NFCETS
        IEL  = FACETS(1,IFC)
        IFACE= ABS(FACETS(2,IFC))
        ITOUCH=FACETS(3,IFC)
        NOD = NUMS(IEL,1)
        IMAT =NUMS(IEL,2+NOD)

        IF (IMAT.le.0) GOTO 701    ! 'skip' (as 'invis')

        IF (ID_FACE(IFACE).eq.1) GOTO 701   ! never 'draw' this side ?
        IF (ID_FACE(IFACE).eq.2) GOTO 777   ! always 'draw' this side ?

        IF (ITOUCH.eq.0) GOTO 777   !  'draw' (as it is a 'boundary')

        IEL2 = FACETS(1,ITOUCH)  ! get info on its toucher'
        IMAT2 = NUMS(IEL2,2+NUMS(IEL2,1))
        IF (IMAT2.gt.0) GOTO 701  ! 'skip' (cos 'touch' is present !)

c---------- now get & transform this facets' nodes --------------------
  777 CONTINUE
      CALL GFACE (NOD,NDIME,ITYPE,IFACE,FS,NN_F,NN_FT, 2 )

C .. ok turn FS into the 'global node numbers .. not just the 'local' !
      I_SCREEN = 0
C      Z_MIN,Z_MAX,Z_MEAN ???
      DO II=1,NN_FT
        I = NUMS(IEL,1+FS(II))
        FS(II) =I
        GPT(1) = GCOORD(1,I)
        GPT(2) = GCOORD(2,I)
        GPT(3) = GCOORD(3,I)
        IF (NLDS.GT.0) THEN
          GPT(1) = GPT(1) +XDS1 *GDISPS(1,IB1+I) +XDS2 *GDISPS(1,IB2+I)
          GPT(2) = GPT(2) +YDS1 *GDISPS(2,IB1+I) +YDS2 *GDISPS(2,IB2+I)
          GPT(3) = GPT(3) +ZDS1 *GDISPS(3,IB1+I) +ZDS2 *GDISPS(3,IB2+I)
        ENDIF
        CALL TRANSFORM (GPT,SPT, PM, IROT, XCEN,YCEN,SC_X,SC_Y,SC_Z)

        DO J=1,3
          SCOORD(J,I) = SPT(J)  ! copy screen position into SCOORD
          SC_F (II,J) = SPT(J)  !  and into 'this' facet too
        ENDDO  

C---------------------- now test for 'off' screen' ---------------------
c.... this should be 'optional' .. or at least fixed for 'buffering'

      IF (IOP_BUF.eq.0) THEN  !--- clip to full image
        IF (SPT(1).gt.iwin_xo.and.SPT(1).lt.iwin_xo+iwin_xw) THEN
         IF (SPT(2).gt.iwin_yo.and.SPT(2).lt.iwin_yo+iwin_yw) I_SCREEN=1
        ENDIF
      ELSE  !-- clip to the 'buffered' image (0->iwin_xw)  ..
        IF (SPT(1).gt.0..and.SPT(1).lt.iwin_xw) THEN
          IF (SPT(2).gt.0..and.SPT(2).lt.iwin_yw) I_SCREEN=1
        ENDIF
      ENDIF

      ENDDO            !-- end-of-loop-a-facets-nodes
      IF (I_SCREEN.eq.0) GOTO 701   !-- skip as no_nodes are on_screen

C------------------ back-face polygon culling --------------------------
      C3 = 0.
      DO I=1,NN_F
        IP1 = MOD (I,NN_F) + 1
        C3 = C3 + (SC_F(I,2)+SC_F(IP1,2)) * (SC_F(IP1,1)-SC_F(I,1))/2.
      ENDDO
      C3 = -C3 * AR  ! weight wrt. the screen topology
      IF (CV(32).EQ.1.AND.C3.GT.0) GOTO 701   ! cull 'fronts'
      IF (CV(32).EQ.2.AND.C3.LT.0) GOTO 701   ! cull 'backs'

C--------- OK this facet IS to be drawn.. so get its 'z'-value ---------
      NFCETS_D = NFCETS_D +1    ! = '# of Facets to draw'
      PL2(NFCETS_D) = IFC       !-- 'index' of this facet in PL1 


       ZVAL = SC_F(1,3)
       IF (cv(17).eq.2) then                    ! furthest point
         DO I=1,NN_FT
           ZVAL = MIN(ZVAL,SC_F (I,3))
         ENDDO
       ELSEIF (cv(17).eq.1) then                ! nearest point
         DO I=1,NN_FT
           ZVAL = MAX(ZVAL,SC_F (I,3))
         ENDDO
       ELSE                                     ! centroid point 
         DO I=1,NN_FT
           ZVAL = ZVAL + SC_F (I,3)
         ENDDO
         ZVAL = ZVAL / REAL(NN_FT)
       ENDIF

       ZBUF(NFCETS_D) = ZVAL

  701 CONTINUE !.......... end-of 'skip-to' (='ignore this facet')
      ENDDO  !............ end of the facet loop (IFC)

      call clock@(time(2))  !--  time to end of transformation

C---------------- depth sort the 'drawable' facets ---------------------
      CALL RSORT@(PL1,ZBUF,INTL(NFCETS_D)) 

      call clock@(time(3))  !-- time to (after) depth-sorting the facets 

C-----------------------------------------------------------------------





C=======================================================================
C========================  draw the new image ==========================
C=======================================================================

C----------- open virtual screen if buffered ---------------------------
      IF(IOP_BUF.NE.0) THEN
        CALL CREATE_SCREEN_BLOCK@ (iwin_xw,iwin_yw,IRESCP,VSCREEN)
        CALL OPEN_VSCREEN@(VSCREEN,IFAIL) 
      ENDIF
C------------------------ blank window ---------------------------------
C .. 'skip' for pen-plotters !
      IF (CV(20).NE.-1.and.IOUT_D.ne.8) call fill_rectangle@ 
     + (iwin_xo, iwin_yo, iwin_xo+iwin_xw, iwin_yo+iwin_yw, CV(20))

C------------------------ draw the picture frame -----------------------
      IF (CV(21).NE.-1)   ! colour may also be the 'picture-# !'
     +  CALL GRATIC (RES(1),RES(2),cv(22),cv(22), CV(21))

C--------------------- now plot the axes -------------------------------
C ... nice to draw the 'back-planes too'                             :-)
C ... ie. loop the '6' sides of the surrounding cuboid
C ... then b_p_c each, then loop each's 'squares' & draw wire-frame/fill
C ... to the nearest 10**size.. ie. square size =1mm/5m/20km etc.
      IF(CV(24).NE.-1)THEN
      DO I=1,3                ! loop the 3 axes
        DO J=1,3
          GPT(J) = COA(J)
          IF (I.EQ.J) GPT(J) = GPT(J) + DIAG/2.  ! get the axis
        ENDDO
        CALL TRANSFORM (COA,SPT,  PM, IROT, XCEN,YCEN,SC_X,SC_Y,SC_Z)
        CALL TRANSFORM (GPT,SPT2, PM, IROT, XCEN,YCEN,SC_X,SC_Y,SC_Z)
        CALL DR_ARROW  (SPT,SPT2,10,15,CV(24))

          IF (CV(27).NE.-1) THEN   ! label the axis
            NTXT = CHAR (ICHAR('w') + I )  ! ie 'x','y','z'
            CALL DRAW_TEXT@ (NTXT,INTS(SPT2(1)),INTS(SPT2(2)),CV(27))
          ENDIF
        ENDDO
      ENDIF

C---------------- set min/max contour values etc. ----------------------
      C_R(1) =  1.e12  ! set the max/min contour values
      C_R(2) = -1.e12  ! to very high/low number

C------------------- loop facets and plot them -------------------------

      DO I_F = 1, NFCETS_D  ! =just the NFCETS that are to be drawn!

        IFC = PL2 (PL1(I_F))
        IEL  = FACETS (1,IFC)
        IFACE= FACETS (2,IFC)
        NOD  = NUMS(IEL,1)
        IMAT = NUMS(IEL,NOD+2)
        ICMAT= IC_MATS(IMAT) 
        IF (IMAT.EQ.0) GOTO 81
        CALL GFACE (NOD,NDIME,ITYPE,IFACE,FS,NN_F,NN_FT, 2)
        DO I=1,NN_FT
          FS(I) = NUMS(IEL,1+FS(I))
        ENDDO

C------------- get screen positions of the nodes of a facet ------------
      DO I=1,NN_FT    !--  hmm in a sub call could just pass 
        DO J=1,3      !--               SCOORD(J,FS(i) ?
          SC_F(I,J) = SCOORD (J,FS(I))
        ENDDO
      ENDDO


C-------------------- 'shrink' algorithm -------------------------------
      IF (CV(31).ne.0) THEN
      DO J=1,3    ! no need for screen 'z' ?
        GPT(J) = 0.
        DO I=1,NOD
          GPT(J) = GPT(J) + GCOORD (J,NUMS(IEL,1+I))
        ENDDO
        GPT(J) = GPT(J) / REAL (NOD)     ! the element 'centroid'
      ENDDO
      CALL TRANSFORM (GPT,SPT, PM, IROT, XCEN,YCEN,SC_X,SC_Y,SC_Z)

      DO J=1,2                 ! no need for screen 'z' ?
        IF (cv(31).GT.0) then     ! as %-age shrinking
          DO I=1,NN_FT
            SC_F(I,J) = SC_F(I,J) - (SC_F(I,J)-SPT(J)) * cv(31)/100.
          ENDDO
        ELSEIF (cv(31).lt.0) THEN  ! as # of pixels to shribk by
          DO I=1,NN_FT
            SC_F(I,J) = SC_F(I,J) - SIGN (REAL(CV(31)),SC_F(I,J)-SPT(J))
          ENDDO
        ENDIF
      ENDDO
      ENDIF

C-------------------get disps of this facet-----------------------------
      DO I=1,NN_FT  ! (I don't always need to do this ?)
        II  = FS(I)
        IB1 = NN * (ILD1-1)
        IB2 = NN * (ILD2-1)
        DO J=1,3
           C_F(I,J) = GCOORD(J,II)
           IF (NLDS.GT.0) THEN  ! also should 'null' DI_F at the top :-)
             DI_F(I,J) = FAC1 *GDISPS(J,IB1+II) + FAC2 *GDISPS(J,IB2+II)
           ENDIF
        ENDDO
      ENDDO
C------------------------ sub-faceting ---------------------------------
C-------------- loop sub-facets ----------------------------------------
C-------- get the geometry of the sub-facet ----
C..>  not really CV(19)**2 as in triangles = sum CV(19)
C... therefore should call a subroutine to set 'up' this sub-facet
C... and another within the loop to return the local co-ords of one
C... of the points (given a 'patch' of what it looks like)

      N_SF = CV(19)
      DO I_SF = 1, N_SF * N_SF

c...................................................................
C ... skip to main menu if a key pressed  ** mouse press ?? **
      EXIT = (KEY.ne.13.and.KEY_WAITING@()) 
      IF (EXIT) GOTO 1119
c...................................................................

C---------------------- no sub-faceting --------------------------------
c----- get the local co-ords of this facets' nodes using WTHATN
c.........! note I am forcing 2D !
        IF (N_SF.eq.-1) THEN   
        IF (.not.( IDT.eq.2.or.idt.eq.4).and.CV(1).ge.14) 
     +        CALL WTHATN (NN_F,2,1,LC_SF,ISC)
c..... note 'skip' if IDT=2/4 (='OFF/NFF'format) as pentagon! etc. ?
c..........  or we are NOT contouring :-)
          NN_SF = NN_F
          DO I=1, NN_SF
            DO J=1,3   
              SC_SF(I,J) = SC_F(I,J)   
            ENDDO
          ENDDO

C--------------------------- sub-faceting ------------------------------
        ELSE
          CALL GET_FAC (LC_SF,ISC, I_SF,N_SF, NN_SF, NN_FT)  ! get LC_SF
          DO I=1,NN_SF              ! loop the nodes of a sub-facet
            CALL SAMPLE (SC_F,ISC,DI_F,ISC,  NN_FT,2,3,
     +                            VAL,VEC, LC_SF,ISC,I,10)
C.. maybe get the face normals at this point too ?
            SC_SF(I,1) = VEC(1)    ! copy the interpolated nodal co-ords
            SC_SF(I,2) = VEC(2)    ! into SC_SF
            SC_SF(I,3) = VEC(3)
          ENDDO
        ENDIF

C------------------ back-face sub-facet culling ------------------------
C .. trash this for now !
      DO J=1,3
        SC_SF (NN_SF+1,J) = SC_SF (1,J)
      ENDDO

C-----------------------------------------------------------------------
C-------------------- select face color option  ------------------------
C-----------------------------------------------------------------------
c ... maybe loop around this twice .. for 2 'layers' eg. conts+mats

      IF (CV(1).LE.12) THEN    !--- non-interpolated colouring -----

      IF (CV(1).EQ.-1) THEN                    ! See-thru
        ICOL = -1 
      ELSEIF (CV(1).EQ.0) THEN                 ! Black (other colours ?) 
        ICOL = cv(2)
      ELSEIF (CV(1).EQ.1) THEN                 ! Material # (indexed)
        ICOL = ICMAT
      ELSEIF (CV(1).EQ.2) THEN                 ! Facet direction
        ICOL = ABS(IFACE) + 1
      ELSEIF (CV(1).EQ.3) THEN                 ! Inside/Outside
        IF (C3.LT.0) ICOL = 2
        IF (C3.GT.0) ICOL = 3
      ELSEIF (CV(1).EQ.4) THEN                 ! Glass front :-)
        IF (C3.LT.0) ICOL = -1                 ! only if no b_p_c (!)
        IF (C3.GT.0) ICOL =  2
      ELSEIF (CV(1).EQ.5) THEN                    ! Element #
        ICOL = INT(IEL / REAL(NEL) *CV(18)) + 2
      ELSEIF (CV(1).EQ.6) THEN                    ! Mean Node #
        COL = 0
        DO I=1,NN_FT
          COL = COL + FS(I)
        ENDDO
        ICOL = COL /REAL(NN_FT) /REAL(NN) * CV(18) + 2
      ELSEIF (CV(1).EQ.7) THEN                    ! sub-facet #
        ICOL = 1+ I_SF
      ELSEIF (CV(1).EQ.8) THEN                    ! Chessb.  sub-facets 
        ICOL = MOD(I_SF,2) + 2
      ELSEIF (CV(1).EQ.9) THEN                    ! Chessboard elements
        ICOL = MOD (IEL,2) + 2

C---------------- 'flat' light-source shading of facet -----------------
      ELSEIF (CV(1).GE.10.AND.CV(1).LE.11)THEN

C... this is FLAT shading so use the screen coords in SC
C... and the 'projected-area' method  :-)
C... note. 'Other' 'Flat' data eg. mean element pore pressures etc.
        C1 = 0. 
        C2 = 0.
        C3 = 0.
        DO I=1,NN_SF
          J  = MOD (I,NN_F) + 1
          C1 = C1+(SC_SF(I,3)+SC_SF(J,3))*(SC_SF(J,2)-SC_SF(I,2)) /2.
          C2 = C2+(SC_SF(I,1)+SC_SF(J,1))*(SC_SF(J,3)-SC_SF(I,3)) /2.
          C3 = C3+(SC_SF(I,2)+SC_SF(J,2))*(SC_SF(J,1)-SC_SF(I,1)) /2.
        ENDDO
        C4 = ABS (C1*C1 + C2*C2 + C3*C3 )
        IF (C4.GT.1.E-12) THEN
        C4 = SQRT(c4)
          A = (LIGHT(1)*C1+LIGHT(2)*C2+LIGHT(3)*C3) /C4  ! light
          B = (                           (-1.)*C3) /C4  ! eye
          COL = SIGN(A,A*B)    ! ie. on the 'same' side as the light
c... hmm should all the colour range to be user-defined
          ICOL = NINT(MIN(MAX(2.,2.+CV(18)*COL),CV(18)+1.))    ! clip
        ELSE
          ICOL = -1   ! zero-area face so skip colouring it
        ENDIF

C --> should insert here the patch for different material palettes
        IF (CV(1).EQ.11) ICOL = (15-ICOL)*16 +IMAT+1 ! eg. palette #2

      ENDIF   ! end of facet flat-colour options

C------------------'flat'-colour in the face ---------------------------
      IF (CV(1).NE.-1) CALL DRAW_POLY (SC_SF,ISC,NN_SF,ICOL,2) 
      ENDIF

C-----------------------------------------------------------------------
C----------- contouring of co-ord/disp/strain etc. ---------------------
c... hmm light-source-shading should use the 'displaced' co-ords
c...     coord/disp/strains will use the 'original' geometry !

      IF (CV(1).GE.14) THEN
C--------------------- light-source shading ----------------------------
        IF (CV(1).lt.20) then    !--- light-source shading

        DO I=1,NN_SF   !----- loop the nodes of a sub-facet (NN_SFT ?)
      
      CALL SAMPLE(SC_F,ISC,DI_F,ISC,NN_FT,2,3,VAL,VEC,LC_SF,ISC,I,CV(1))

        T = sqrt(vec(1)**2+vec(2)**2+vec(3)**2) ! noramalise face normal
        IF (t.gt.0.) THEN
          A = (LIGHT(1)*VEC(1)+LIGHT(2)*VEC(2)+LIGHT(3)*VEC(3)) /T
          B = (                                   (-1.)*VEC(3)) /T  !eye
          COL = SIGN(A,A*B)        ! ie. on the 'same' side as the light

C         write(19,'(7f12.4)') col,a,b,t,vec(1),vec(2),vec(3)
c... hmm should all the colour range to be user-defined
          ACOL(I) = MIN (MAX(1.,(1.+CV(18)*COL)),cv(18)+1.) ! clip to 2.--> 15 ?
        ELSE
          ACOL(I) = 9   ! if zero-area draw in grey ?(col=9)
        ENDIF
        IF (CV(1).EQ.15)  ACOL(I) = ACOL (I) + (IMAT-1) * 16  ! NO!
      ENDDO
C------------------ contouring of strains, disps, etc. -----------------
      ELSE
        DO I=1,NN_SF   !----- loop the nodes of a sub-facet (NN_SFT ?)
      
      CALL SAMPLE(C_F,ISC,DI_F,ISC,NN_FT,2,3,VAL,VEC,LC_SF,ISC,I,CV(1))
        DIS = val
        C_R (1) = MIN (C_R(1),DIS)         !  min contour value found
        C_R (2) = MAX (C_R(2),DIS)         !  max contour value found

        IF (C_R(6).gt.0) DIS = ABS(DIS)    ! take absolute if desired
        DIS = (dis-c_r(3))/(c_r(4)-c_r(3)) ! scale to 0.--> 1.

c       ....logarithic scaling  .... (skip if 'log factor' = 1 )
        IF (int(c_r(5)*10).ne.10)  DIS = DIS ** (1./c_r(5)) 
c... hmm should all the colour range to be user-defined
        N_CONTS = CV(18)

c       ACOL(I) = MIN (MAX(0.,(CONTS*DIS)),CONTS ) + 1.
c.. so valid values to contour are 0. --> N_CONTS
C.. but value (not to contour) will lie outside these limits
        ACOL(I) = DIS * N_CONTS 
      ENDDO

      ENDIF   ! end of contouring / lighting

      CALL GOURAUD (SC_SF,ISC,ACOL, N_CONTS,NN_SF,CV(15),1,ICMAT)  !-- filled
      CALL GOURAUD (SC_SF,ISC,ACOL, N_CONTS,NN_SF,CV(14),2,ICMAT)  !-- +lines

      ENDIF   ! endif of an 'interpolated' colour scheme

C-------------------- edge the sub-facet -------------------------------
      IF (CV(5).NE.-1) CALL DRAW_POLY (SC_SF,ISC,NN_SF,CV(5),1) ! edges
C.. OK have a CV(6) (say) for edging the 'exterior' only ?

C------------------ end of the sub-faceting ----------------------------
   82 CONTINUE    ! skip-to index for backface sub-facet cull (icol=-1)
      ENDDO
C---------------------- edge facet -------------------------------------
C... need to sub-sample for 'smooth' edges 
C...  -also identify the geometry/material 'edges'
      IF (CV(4).NE.-1) CALL DRAW_POLY (SC_F,ISC,NN_F,CV(4),1)

C-------------------- 'undeformed' facet -------------------------------
      IF (CV(16).NE.-1) THEN
        DO I=1,NN_F
          II = FS(I)
          DO J=1,3
            GPT(J) = GCOORD(J,II)   ! start-point
          ENDDO
          CALL TRANSFORM (GPT,SPT, PM, IROT, XCEN,YCEN,SC_X,SC_Y,SC_Z)
          DO J=1,3
            SC_SF(I,J) = SPT(J)
          ENDDO
        ENDDO
        CALL DRAW_POLY (SC_SF,ISC,NN_F,CV(16),1) !  :-)
      ENDIF                        

C----------------------- draw nodes ------------------------------------
      IF (CV(9).NE.-1) THEN
        II = NINT (cv(10)*res(1)/640.)  !-- for 'large pcx/plt's etc.
        DO I=1,NN_FT
c... note a 'point' circle doesn't work! .. patch with 'draw point'?
c... or let CV(9) be the diameter (not radius) so move centre ?
          CALL FILL_ELLIPSE@
     +     (INTS(SC_F(I,1)),INTS(SC_F(I,2)),II,II,CV(9))
        ENDDO
      ENDIF
C-------------------- node numbering -----------------------------------
C .... need an option to draw 'centred' on the node
      IF (CV(11).NE.-1) THEN
        DO I=1,NN_FT       ! ie. 'all' the nodes on this facet
          WRITE (NTXT,'(I4)') FS(I)
          CALL TRIM (NTXT)           ! ie. right-justify
          CALL DRAW_TEXT@ (NTXT,INTS(SC_F(I,1)),INTS(SC_F(I,2)),CV(11))
        ENDDO
      ENDIF

C------------------element numbering -----------------------------------
      IF (CV(13).NE.-1) THEN
        WRITE (NTXT,'(I4)') IEL
        CALL TRIM (NTXT)  !       right-justify
        IXP= 0
        IYP= 0            ! find the centroid of the face
        DO J=1,NN_F          ! (or NN_FT ?)
          IXP = IXP + SC_F(J,1)
          IYP = IYP + SC_F(J,2)
        ENDDO
        CALL DRAW_TEXT@ (NTXT,IXP/NN_F,IYP/NN_F,CV(13))
      ENDIF
C-------------------- displacement vectors -----------------------------
      IF (CV(12).NE.-1) THEN
        DO I=1,NN_FT
          II = FS(I)
          DO J=1,3
            GPT(J) = GCOORD(J,II)   ! start-point
          ENDDO
          CALL TRANSFORM (GPT,SPT, PM, IROT, XCEN,YCEN,SC_X,SC_Y,SC_Z)
          DO J=1,3
            IF (NLDS.GT.0)          ! end-point
     +      GPT(J) = GPT(J) + GDISPS(J,II+IB1) * FAC1 * FACT
     +                      + GDISPS(J,II+IB2) * FAC2 * FACT
          ENDDO
          CALL TRANSFORM (GPT,SPT2, PM, IROT, XCEN,YCEN,SC_X,SC_Y,SC_Z)
          CALL DR_ARROW  (SPT,SPT2,10,20,CV(12))
        ENDDO
      ENDIF
C---------------------- end of this facet ------------------------------
   81 CONTINUE
      ENDDO    

C ....  hmm .. if 'inerupted' then we won't reach here anyway !

 1119 IF (IOP_BUF.EQ.1)THEN
        IF (.NOT.EXIT)         ! inhibit if only partly drawn
     +    CALL VSCREEN_TO_SCREEN@ (INTS(iwin_xo),INTS(iwin_yo),0,IFAIL) 
        CALL CLOSE_VSCREEN@()  ! do I need to do this ? ,_ yes ! -why? 
        CALL RETURN_STORAGE@(VSCREEN)  ! remove it completely
        IF (IFAIL.NE.0) PRINT*,'Ifail=',IFAIL,' VS_2_S' 
      ENDIF
      IF (EXIT) GOTO 1111   ! jumo back to the top of the program

C--------- Draw the 'extra-lines (eg. model boundaries) ----------------
      DO I=1,N_X_L
        GPT(1) = XTRA_L(1,I)
        GPT(2) = XTRA_L(2,I)
        GPT(3) = 0.
        CALL TRANSFORM (GPT,SPT, PM, IROT, XCEN,YCEN,SC_X,SC_Y,SC_Z)
        X(1) = SPT(1)
        Y(1) = SPT(2)

        GPT(1) = XTRA_L(3,I)
        GPT(2) = XTRA_L(4,I)
        GPT(3) = 0.
        CALL TRANSFORM (GPT,SPT, PM, IROT, XCEN,YCEN,SC_X,SC_Y,SC_Z)
        X(2) = SPT(1)
        Y(2) = SPT(2)

        CALL DR_PRIM (X,Y,2,1,1)       !--- draw as a closed polygon
      ENDDO
C-------------- finished model so draw borber etc. ---------------------
C... ok loop and draw around the current window
C .. need to invert for the plotter (its origin is the bottom left)
      IF (CV(28).GE.0) THEN
        X(1) = IWIN_XO
        X(2) = IWIN_XO + IWIN_XW
        X(3) = IWIN_XO + IWIN_XW
        X(4) = IWIN_XO

        Y(1) = IWIN_YO
        Y(2) = IWIN_YO 
        Y(3) = IWIN_YO + IWIN_YW
        Y(4) = IWIN_YO + IWIN_YW
        CALL DR_PRIM (X,Y,4,CV(28),1)     !--- draw as a closed polygon
      ENDIF

C--------------------- auto-'draw-next image' options ------------------
      IF (CV(37).eq.1)THEN 
        LD = LD+1           ! load step
        CV(41) = CV(41) +1  ! pic#
        IF (CV(41).le.CV(42)*CV(43) ) GOTO 901 ! redraw
        CV(41) = 1
        LD = 1   !  reset to beginning ready for a re-draw
      ELSEIF (CV(37).eq.2)THEN 
        CV(38) = cv(38)+1   ! view direction
        CV(41) = CV(41) +1  ! pic#
        IF (CV(41).le.CV(42)*CV(43) ) GOTO 901 ! redraw
        CV(41) = 1
        CV(38) = 1  ! reset to the first picture ready for a re-draw
      ENDIF


C--------------------- end of image drawing ----------------------------
 1120 continue   ! 'exit'- label (if a key/mouse-button pressed) ?? 
      call clock@(time(4))

C-------------  dump image to printer/plotter file etc.-----------------
      IF (IOUT_D.ne.1) THEN

      IF (IOUT_D.EQ.3) THEN             ! PCX output to file
        CALL VSCREEN_TO_PCX@ (FILE_PCX,IFAIL)
        CALL PCX_FIX (FILE_PCX) 

        CALL CLOSE_VSCREEN@ ()
        CALL RETURN_STORAGE@ (VSCREEN)

      ELSEIF (IOUT_D.EQ.8) THEN         ! HP7550 to file
        CALL CLOSE_PLOTTER@()

      ELSEIF (IOUT_D.EQ.4) THEN        ! HP-LJ to file
        CALL CLOSE_GRAPHICS_PRINTER@() 
        IF (FILE_PCX.eq.'LASER') THEN      ! auto print to HP-UX :-)
          CALL CISSUE@('pr_hpux',ifail) ! sends & deletes 'LASER'
        ENDIF
        IF (FILE_PCX.eq.'PRINT') THEN      ! auto print to HP-UX :-)
          CALL CISSUE@('print print',ifail) ! sends & deletes 'LASER'
        ENDIF
      ELSEIF (IOUT_D.EQ.5) THEN        !--------- PostScipt ------------
        CALL CLOSE_PS_PRINTER () 
        IF (FILE_PCX.eq.'LASER') THEN      ! auto print to HP-UX :-)
          CALL CISSUE@('pr_hpux',ifail) ! sends & deletes 'LASER'
        ENDIF
        IF (FILE_PCX.eq.'PRINT') THEN      ! auto print to HP-UX :-)
          CALL CISSUE@('print print',ifail) ! sends & deletes 'LASER'
        ENDIF
      ELSEIF (IOUT_D.EQ.6) THEN        ! pro-printer to LPT1
        CALL CLOSE_GRAPHICS_PRINTER@() 
      ENDIF

      CALL SOUND@ (1024,1)
      RES(1) = IRESX     !-- (this is back to 'full-screen') ?
      RES(2) = IRESY     !---   (ie. no menus ?)
      IOUT_D = 1         !---- revert output to the screen 
      ENDIF
C-------------------------- s t a t u s --------------------------------
C...  show other anotation too !
C..... eg. load step # , viewing angle, light angle, disp scale, etc.

      IF (CV(34).eq.9) THEN            
        call clock@(time(5))
        CALL POST_WIDGET (TABLE,TXT(menup(9)),menup(9))

        write(line,'(a)')'   < '//
     +      FILE_DAT(idir_name+1:leng(file_dat))//' >'
          call draw_text@(line(1:leng(line)),490S, 70S,14S)
        write(line,'(a,i5)')  '   NN =',NN
          call draw_text@(line(1:leng(line)),490S,115S,14S)
        write(line,'(a,i5)')  '   NEL=',NEL
          call draw_text@(line(1:leng(line)),490S,130S,14S)
        write(line,'(a,i5)')  '#loads=',NLDS
          call draw_text@(line(1:leng(line)),490S,145S,14S)
        write(line,'(2(a,i4))')  '#facets=',NFCETS_D,'/',NFCETS
          call draw_text@(line(1:leng(line)),490S,160S,1S)

        write(line,'(a,f10.3)')  '--- timings ---'
          call draw_text@(line(1:leng(line)),490S,300S,1S)
        write(line,'(a,f10.3)')  'total =',time(5)-time(1)
          call draw_text@(line(1:leng(line)),490S,315S,15S)
        write(line,'(a,f10.3)')  'transf=',time(2)-time(1)
          call draw_text@(line(1:leng(line)),490S,330S,1S)
        write(line,'(a,f10.3)')  'd.sort=',time(3)-time(2)
          call draw_text@(line(1:leng(line)),490S,345S,1S)
        write(line,'(a,f10.3)')  'draw  =',time(4)-time(3)
          call draw_text@(line(1:leng(line)),490S,360S,1S)
        write(line,'(a,f10.3)')  'buffer=',time(5)-time(4)
          call draw_text@(line(1:leng(line)),490S,375S,1S)

        write(line,'(a,f10.3)')  'reading=',time(11)-time(10)
          call draw_text@(line(1:leng(line)),490S,390S,15S)
        write(line,'(a,f10.3)')  'FSTRIP =',time(12)-time(11)
          call draw_text@(line(1:leng(line)),490S,405S,15S)
      ENDIF
c-------------------------- animation ----------------------------------
      IF (IOP_ANIM(1).GT.0)THEN

        i = iop_anim(2)     ! = # of frames to reach final load increment
        if (i.ne.0) then    
          ld = ld + ldd     ! increment the load
          if (i.gt.0.and.ld.gt.nlds) ld = 1.    !----- wrap
          if (i.lt.0.and.ld.ge.nlds) ldd = -abs(ldd) !---- bounce down! 
          if (i.lt.0.and.ld.le.1   ) ldd =  abs(ldd) !---- bounce up! 
        endif
c.... also sinusoidal animation of a single load step ??
        XME = XME + iop_anim(3)  ! why were these once halved ??
        YME = YME + iop_anim(4)
        XML = XML + iop_anim(5)
        YML = YML + iop_anim(6)

c---- save the image if desired (to a screen_block/pcx)
        if (iop_anim(7).eq.0) then
          continue  ! a no-op :-)
        else
        iop = mod (iop_anim(7),10)    ! the 'bit-left over'

          if (ipic_n.gt.iop_anim(9)) iop_anim(1) = 0 ! kill at the end :-)
          write(*,'(f4.2)') LD
          IF(IOP.eq.1) CALL GET_SCREEN_BLOCK@
     +       (iwin_xo+0,iwin_yo+0,iwin_xo+320-1,iwin_yo+200-1,BUFFER)
          IF(IOP.eq.2) CALL GET_SCREEN_BLOCK@
     +                                       (0,0,640-1,480-1,BUFFER)
          IF(IOP.eq.3) CALL GET_SCREEN_BLOCK@
     +                                       (0,0,480-1,480-1,BUFFER)
          IF(IOP.eq.4) CALL GET_SCREEN_BLOCK@
     +   (iwin_xo,iwin_yo,iwin_xo+iwin_xw-1,iwin_yo+iwin_yw-1,BUFFER)

        print*,C_T//'picture #=',ipic_n  !-Only write this AFTER capture

        if (iop_anim(7).lt.10) then  !--- dump to a PCX file
          print*,C_T//'picture #=',ipic_n
          write(FILE_PCX,'(a,a,i2.2,a)')
     +           'AN\',FILE_ROOT(1:leng(FILE_ROOT)),ipic_n,'.pcx'
          CALL SCREEN_BLOCK_TO_PCX@ (FILE_PCX,BUFFER,IFAIL)
          CALL RETURN_STORAGE@(BUFFER)
          CALL PCX_FIX (FILE_PCX)
        else
          IF (image(ipic_n).ne.-1) CALL RETURN_STORAGE@(image(ipic_n))
          image(ipic_n) = buffer
        endif
        ipic_n = ipic_n +1   ! increment picture number for next time
        if (ipic_n.eq.iop_anim(9)) iop_anim(1) = 0 ! kill at the end :-)
        endif  ! end of 'save-image-to-pcx-or-memory'

        GOTO 901  !-- straight to re-draw ? (bypass event handling ?!)
      ENDIF


C-------> goto the start of the loop --> --> --> --> --> --> --> --> -->
      GOTO 1111
C--------------------> E X I T   P R O G R A M <------------------------
  999 continue   ! .. if crashed out of 'select_file@()'
      CALL TEXT_MODE@()
      END
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE ISWAP (I,J)
C   --- swaps the 2 integer arguements
        K = I
        I = J
        J = K
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
        if (mod(i_sf,2).eq.mod(ix,2))then    ! an 'up' pointing triangle
          LC_SF(2,1) = LC_SF(1,1) - delta
          LC_SF(2,2) = LC_SF(1,2) + delta
          LC_SF(3,1) = LC_SF(1,1) - delta
          LC_SF(3,2) = LC_SF(1,2)
        else
          LC_SF(2,1) = LC_SF(1,1)
          LC_SF(2,2) = LC_SF(1,2) + delta
          LC_SF(3,1) = LC_SF(1,1) - delta
          LC_SF(3,2) = LC_SF(1,2) + delta
        endif

      ELSE
         PRINT*,'*** WARNING: NEN=',NN_FT,' not known in GET_FAC!'
      ENDIF
      RETURN
      END

C-----------------------------------------------------------------------
       SUBROUTINE WTHATN2 (NEN,NDIME,ITYPE, LC_SF,ILC_SF)
c
c     this selectively calls WTHATN (ie. if 'different')
c     to find the local co-ordinates of an element (facet)
c

      REAL LC_SF(ILC_SF,*)
c-------------------------- the 'last time' values ---------------------
      SAVE NEN_OLD, NDIME_OLD, ITYPE_OLD
      DATA NEN_OLD, NDIME_OLD, ITYPE_OLD /0,0,0/  ! zap initially

      IF (NEN.eq.NEN_OLD.and.NDIME.eq.NDIME_OLD    ! no work to do
     +          .and.ITYPE.eq.ITYPE_OLD) RETURN

      NEN_OLD =  NEN
      NDIME_OLD= NDIME
      ITYPE_OLD= ITYPE
      CALL WTHATN (NEN,NDIME,ITYPE,LC_SF,ILC_SF)
      RETURN
      END
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
      PARAMETER (IDER=2,ISMP=2,IDERIV=3,IJAC=3,ICOORD=20,IBEE=20)
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
      CALL SHAPE (NDIME,NEN,1, DER,IDER,FUN,SMP,ISMP,1)

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
            IF (IOP.le.29) VEC(J) = VEC(J) + FUN(I) * SC  (I,J) ! coord
            IF (IOP.ge.30) VEC(J) = VEC(J) + FUN(I) * DISP(I,J) ! disp
          ENDDO
        ENDDO
        IF(IOP_2.EQ.0) val = sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
        IF(IOP_2.EQ.1) val = vec(1)
        IF(IOP_2.EQ.2) val = vec(2)
        IF(IOP_2.EQ.3) val = vec(3)

c-----------------------------------------------------------------------
      ELSEIF (IOP.ge.40.and.iop.le.59) THEN  

      CALL MATMUL (DER,ider,COORD,icoord,JAC,ijac,2,NEN,2)
c     CALL INVERT (JAC,IJAC,DET,2L)
      DET = jac(1,1)*jac(2,2) - jac(2,1) * jac(1,2)
      IF (abs(det).lt. 1.e-10) then
        val =0.
        return
      endif
      CALL TWOBY2 (JAC,IJAC,JAC1,IJAC,DET)

      CALL MATMUL (JAC1,IJAC,DER,IDER,DERIV,IDERIV,2,2,NEN)

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
        CALL FMDEPS (DEE,3,1.E6,0.3)    !-- is this plane-stress ?
        CALL MVMULT (DEE,3, EPS,3,3, STRESS)
        DO J=1,3
          EPS(J) = STRESS(J)
        ENDDO
      ENDIF
C----------------------- get the principal values ----------------------
      
      IF (IOP_2.GE.6 .and. IOP_2.LE.8) THEN   ! if INVAR: put shear in 3
        IF (IOP.LT.50) EPS(3) = EPS(3) /2.          !-- divide by 2
        CALL INVAR (EPS, EPS(6),EPS(7),EPS(8))
c       CALL INVAR (EPS, sigm,  DSBAR ,THETA )
      ELSEIF (IOP_2.EQ.9) THEN
        EPS(9) = DET
      ENDIF
      VAL = EPS (IOP_2)

c-----------------------------------------------------------------------
      ENDIF   ! end of IOP option_choice
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE READ_MENU (TABLE,TXT,MOPS,MENUP, iresx,iresy)
c
c     this reads the menu from the data file 'menu.txt'
c     and store menu label in 'txt' and info in table
c     MENUP = the start positions of each (sub-) menu
c
c     ... version #2, (pl7e)  relative positions etc.

      integer table(11,mops), menup(*)
      character txt(mops)*(*), text*3

c................ search for the 'menu' file ...........................

      open (11,file='pl7e.mnu',STATUS='old',iostat=iostat)
      if(iostat.ne.0) open (11,file='c:\users\dan\gen\pl7e.mnu'
     +                       ,STATUS='old',iostat=iostat)
      if(iostat.ne.0) open (11,file='c:\util\pl7e.mnu'
     +                       ,STATUS='old',iostat=iostat)
      if(iostat.ne.0)then
        print*,'*** MENU file PL7e.MNU could not be found !'
        return
      endif
c-----------------------------------------------------------------------
      read (11,'(6(/))')
      do i=1,mops
        read(11,*,end=91) (table(j,i),j=1,9),text,txt(i)
        do j=1,2
          table(9+j,i) = ichar(text(j:j))
        enddo

c-------------------- patch to reverse colours 1 and 2 -----------------
        do j=6,9
          if (table(j,i).eq.1) then
            table(j,i) = 0
          elseif (table(j,i).eq.0) then
            table(j,i) = 1
          endif
        enddo
c---------------- store start point of menus, etc. ---------------------
        if (table(10,i).eq.ichar('M')) then
          table(10,i) = ichar('m')
          read (text,'(1x,i2)') ii
          menup(ii) = i

          menu_xo = table(1,i)   ! set this menu's position
          menu_yo = table(2,i)
          menu_xl = table(3,i)
          menu_xl = table(4,i)
        else                     ! buttons are all 'relative'
C------------------- automatic spacing of buttons ----------------------
          if (table(1,i).eq. 0) then
             table(1,i) = table(1,i-1)         ! no x-move
          elseif (table(1,i).eq.-1) then
             table(1,i) = table(3,i-1) + 3     ! auto x-move
          else
            table(1,i) = table(1,i) + menu_xo  ! rel. to menu-pane
          endif
  
          if (table(2,i).eq. 0) then
            table(2,i) = table(2,i-1)          ! no y-move
          elseif (table(2,i).eq.-1) then
            table(2,i) = table(2,i-1) + 20     ! auto y-move
          else
            table(2,i) = table(2,i) + menu_yo  ! rel. to menu-pane
          endif

          if (table(3,i).eq.-1) table(3,i) = leng(txt(i)) * 8 + 8 !text

          do j=leng(txt(i)),1,-1   !---------- strip padding characters
            if (txt(i)(j:j).eq.'~') txt(i)(j:j) = ' '
          enddo
        endif      !-- end-of button/menu
c------------- store pointers to menus  (etc.?) ------------------------
        read (text,'(1x,i2)') ii 
        table(11,i) = ii
        if (table(10,i).eq.32) table(10,i) = -1 !-- if no-op store as -1


C------ convert tile lengths into their end-points ---------------------
        table(3,i) = table(1,i) + table(3,i) -1 
        table(4,i) = table(2,i) + table(4,i) -1 

C--- rescale the menu to the current video resolution (?)
c        table(1,i) = table (1,i) * real(iresx) / 640.
c        table(2,i) = table (2,i) * real(iresy) / 480.
c        table(3,i) = table (3,i) * real(iresx) / 640.
c        table(4,i) = table (4,i) * real(iresy) / 480.

      enddo
   91 nops = i-1 
      close (11)
      end
c-----------------------------------------------------------------------
      subroutine post_widget (table,txt,i)
c
c     this posts a menu item 'txt' at position/colors held in 
c     line 'i' of 'table'
c                                              Dan Kidger  29-11-91
c--> need to extend to include 'control' characters eg. colours, font,
c    'dials', 'sliders', etc. !

      integer table(11,*)
      character txt*(*)

      ixf  = table(1,i)
      iyf  = table(2,i)
      ixt  = table(3,i)
      iyt  = table(4,i)
      ibdr = (iyt - iyf - 14 ) /2
      icolt = table(6,i)      !- text
      icolb = table(7,i)      !- background
      icol1 = table(9,i)      !- top+left border
      icol2 = table(8,i)      !- bot+right border
c      icol1 = table(8,i)
c      icol2 = table(9,i)
      do j=0,table(5,i)-1
        call draw_line@ (ixf+j,iyt-j,ixf+j,iyf+j,icol1)
        call draw_line@ (ixf+j,iyf+j,ixt-j,iyf+j,icol1)
        call draw_line@ (ixt-j,iyf+j,ixt-j,iyt-j,icol2)
        call draw_line@ (ixt-j,iyt-j,ixf+j,iyt-j,icol2)
      enddo
      if (leng(txt).eq.0) icolt = -1
      if (icolb.ge.0)
     + call clear_screen_area@ (ixf+j,iyf+j,ixt-j,iyt-j,icolb)
c... if text color is not black then draw a black shadow around it :-)
      if (icolt.ne.1)
     + call draw_text@ (txt(1:leng(txt)),ixf+j+1, iyf+ibdr+j-2, 1)
      if (icolt.ne.-1)
     + call draw_text@ (txt(1:leng(txt)),ixf+j+2, iyf+ibdr+j-1, icolt)
      end
c-----------------------------------------------------------------------
      SUBROUTINE Q_MENU (IH,IV, TABLE,MENUP,MENUN, ITEM,KEY,KEY2,XM,YM)
c
c     this returns the menu item (and sub-postion) of a selection
c     ih,iv = mouse position (screen)
c     key1,2 = selection code (and sub-code) (= -1) if 'none')
c     xm,ym = fractional position within a menu item
c
      integer table(11,*),menup(*)

      do i = menup(menun+1)-1,menup(menun),-1
        if (iv.ge.table(2,i)) then
          if (iv.le.table(4,i)) then
            if (ih.ge.table(1,i)) then
              if (ih.le.table(3,i)) then
                key  =   table(10,i)
                key2 =   table(11,i)
                if (key.ne.-1) then
                  item = i
                  xm = real(ih-table(1,i)) / real(table(3,i)-table(1,i))
                  ym = real(iv-table(2,i)) / real(table(4,i)-table(2,i))
                  return
                endif  
              endif  
            endif  
          endif  
        endif  
      enddo
c .. nothing found so crash out 
      item=  1
      key =  0
      key2 =  -99
      xm  = .5
      ym  = .5
      return
      end
C----------------------------------------------------------------------
      SUBROUTINE SET_PAL(IOP,ir,ig,ib,ncol)
C
C     set palette option  1 = simple, 2= greyscale
C
C .. SET COLORS 1-5 = black,white,d.green,l.green, brown,red
C ... -1 = user defined , -2 = 'standard' ,-3 = 'firestorm'
C      0..16 = normal colors
C --> need to make this do all the work ?? .. or keep input in the 
C     main prog ??

      COMMON /PALETTE/PAL
      INTEGER PAL (3,0:255)
      INTEGER   PAL1(3,0:15), PAL2(3,15) 
      REAL PI
      DATA PAL1/255,255,255, 0,0,0, 0,107,0, 100,30,30, 127,127,0, 
     + 200,20,20, 90,0,160, 160,0,160, 127,127,127, 63,63,63,
     + 0,255,0, 40,40,255, 255,40,255, 255,0,0, 255,127,0, 255,255,0/

      DATA PAL2/ 9,147,227, 1,177,204, 1,204,177, 9,227,147, 23,243,116,
     +45,253,85, 71,255,57, 100,249,33, 131,236,15, 162,216,4, 191,191,0
     +,216,162,4, 236,131,15, 249,100,33, 230,0,0/      
      DATA PI/3.14159265/
C----------------------------------------------------------------------
      DO I=1,15
        CALL SET_PALETTE@ (I,I)
      ENDDO
        CALL SET_PALETTE@ (17,8)
      IF (IOP.LE.0) THEN
C .. single colour redefinition
        PAL(1,-IOP) = IR
        PAL(2,-IOP) = IG
        PAL(3,-IOP) = IB
      ELSE IF(IOP.EQ.1) THEN
C .. shade range redefinition
        IOFF= 6
        DO I=2,15
          PAL(1,I) = IR * (I+IOFF) /(15+IOFF)
          PAL(2,I) = IG * (I+IOFF) /(15+IOFF)
          PAL(3,I) = IB * (I+IOFF) /(15+IOFF)
        ENDDO
      ELSE IF(IOP.EQ.2) THEN
C ... standard 'material' colours
        DO I=0,15
          DO J=1,3
            PAL(J,I) = PAL1(J,I)
          ENDDO
        ENDDO
      ELSE IF(IOP.EQ.3) THEN
C ... contouring 'rainbow' colours
        DO I=2,15
          DO J=1,3
            PAL(J,I) = PAL2(J,I)
          ENDDO
        ENDDO

      ELSE IF(IOP.EQ.4) THEN   ! 'hot-iron' = black->red->yellow->white
        DO I=2,2+ncol-1
          hue = real(I-2) / real(ncol-1)   ! hue = 0.-->1.
          PAL(1,I) = min(max(0,nint(255.* 3.*(hue+.03    ))),255) ! note 
          PAL(2,I) = min(max(0,nint(255.* 3.*(hue-.333333))),255) ! clipping
          PAL(3,I) = min(max(0,nint(255.* 3.*(hue-.666667))),255)
        ENDDO
      ELSE IF(IOP.EQ.5) THEN    ! should really be '7' !
C .. grey scale
        DO I=2,2 + ncol-1
          PAL(1,I) = 255. * (I+1) /ncol
          PAL(2,I) = 255. * (I+1) /ncol
          PAL(3,I) = 255. * (I+1) /ncol
        ENDDO
      ELSE IF(IOP.EQ.6) THEN
C ... reverse palette ......................................
        DO I=2,(2+ncol-1)/2 +1
          DO J=1,3
            CALL ISWAP (PAL(J,I),PAL(J,ncol+3-I) )
          ENDDO
        ENDDO
      ELSE IF(IOP.EQ.7) THEN
C ... reverse black and white ......................................
        DO J=1,3
          CALL ISWAP (PAL(J,0),PAL(J,1) )
        ENDDO
      ELSE IF(IOP.EQ.8) THEN                                       
C ... reverse alternate colours ....................................
        DO I=2,2 + ncol-1
          DO J=1,3
            CALL ISWAP (PAL(J,2+I),PAL(J,2+I+1) )
          ENDDO
        ENDDO
      ELSE IF(IOP.EQ.9) THEN
C ... cycle up ....................................
        DO J=1,3
          DO I=15,3,-1 
            PAL(J,I) = PAL(J,I-1)
          ENDDO
            PAL(J,2) = PAL(J,15)
        ENDDO
      ELSE IF(IOP.EQ.10) THEN                                       
C ... cycle down ..................................
        DO J=1,3
          DO I=2,14 
            PAL(J,I) = PAL(J,I+1)
          ENDDO
            PAL(J,15) = PAL(J,2)
        ENDDO
      ELSE IF(IOP.EQ.11) THEN                                       
C ... spectrum wheel ..................................
C--- adjust saturation ??
       DO I=0,ncol-1
         RED = real(I) /real(ncol) * 2 * PI
         PAL (1,2+I) = 128. + 127. * COS (RED)
         PAL (2,2+I) = 128. + 127. * COS (RED - 2.* PI/3.)
         PAL (3,2+I) = 128. + 127. * COS (RED + 2.* PI/3 )
       ENDDO

      ELSEIF (IOP.EQ.12) THEN    ! invert all the other colours :-)
        DO I=2, 2 + ncol-1
          DO J=1,3
            PAL(J,2+I) = 255 - PAL(J,2+I)
          ENDDO
        ENDDO
      ENDIF
      IF (NCOL.EQ.14) THEN   !-- shade down the other colors :-)
        DO I=0,15
          DO K=0,15
            I2 = K*16 + I
            DO J=1,3   
              PAL(J,I2) = PAL(J,I) * (16+5-K) / (16.+5.)
            ENDDO
          ENDDO
        ENDDO
      ENDIF  
      CALL SET_ALL_DACS ()
      END

C----------------------------------------------------------------------
      SUBROUTINE SET_ALL_DACS()
C... sets all the video dacs ...      16 or 256 colours ???
      COMMON /PALETTE/PAL
      INTEGER PAL(3,0:255)
      INTEGER*1 PAL1 (3,0:255)
      INTEGER*4 ADDRESS

      CALL DYNT@('SET_VIDEO_DAC_BLOCK@',ADDRESS)  ! does it exist ?

      IF (ADDRESS.eq.0) THEN    ! s--l--o--w method (old FTN77)
      DO I=0,15
        CALL SET_VIDEO_DAC@( I, PAL(1,I)/4,PAL(2,I)/4,PAL(3,I)/4)
      ENDDO

      ELSE    ! ->-> fast -> -> method  (new FTN77)
        DO I=0,255
          DO J=1,3   
            PAL1(J,I) = PAL(J,I)/4 
          ENDDO
        ENDDO
        CALL SET_VIDEO_DAC_BLOCK@( 0, 256, PAL1)
      ENDIF

      END
C----------------------------------------------------------------------
      SUBROUTINE SET_SHADE (LINE)
C
C     this sets a shade range palette
C
      CHARACTER LINE*(*)
      COMMON /PALETTE/PAL
      INTEGER PAL(3,0:255)

c     call read_edited_line@ (line,0,1,15,ifail)
      read (line,*) icf,ict,gama, rf,gf,bf,  rt,gt,bt
c ..>  gama correction not yet implimented
      do i = icf,ict
        fact = real (ict - i ) / real (ict-icf)
        pal(1,i) = (rf * fact + rt * (1.-fact)) * 255
        pal(2,i) = (gf * fact + gt * (1.-fact)) * 255
        pal(3,i) = (bf * fact + bt * (1.-fact)) * 255
      enddo
      CALL SET_ALL_DACS ()
      END
C----------------------------------------------------------------------
      SUBROUTINE PCX_FIX (FILE_PCX)
C
C     this 'corrects' the palette in a PCX file to be the same as
C     on the screen
C
      COMMON /PALETTE/PAL
      INTEGER PAL(3,0:255)
      INTEGER*4 NBYTES_R
      CHARACTER FILE_PCX*(*), dummy*16
      CALL OPENRW@(FILE_PCX,IO,IERROR)
        CALL READF@(dummy,IO,16L,NBYTES_R,IERROR)
      DO I=0,15
        CALL WRITEF@ (CHAR(PAL(1,I)),IO,1L,IERROR)
        CALL WRITEF@ (CHAR(PAL(2,I)),IO,1L,IERROR)
        CALL WRITEF@ (CHAR(PAL(3,I)),IO,1L,IERROR)
      ENDDO
      CALL CLOSEF@(IO,IERROR)
      END
C-----------------------------------------------------------------------

      SUBROUTINE TRANSFORM (GPT,SPT, PM, IROT, XCEN,YCEN,SC_X,SC_Y,SC_Z)
C       
C     this transforms a point in GPT ('object' co-ordinates)
C                           into SPT ('screen' co-ordionates)
C     PM() is the homogeneous transformation matrix (from CAMERA q.v.)
C     if IROT = 1 then the X and Y screen axes are swapped (ie. X=up)
C     XCEN,YCEN = screen postion of the centre of the object 
C     SC_X,SC_Y,SC_Z are screen scaling factors 
C      (all the same unless the screen is anisotropic) 
C 
      REAL GPT(*), SPT(*), PM(4,4)

c....................... transform the 3D position .....................
      SPT4 =  GPT(1)*PM(4,1)+GPT(2)*PM(4,2)+GPT(3)*PM(4,3)+PM(4,4)
      SPT1 = (GPT(1)*PM(1,1)+GPT(2)*PM(1,2)+GPT(3)*PM(1,3)+PM(1,4))/SPT4
      SPT2 = (GPT(1)*PM(2,1)+GPT(2)*PM(2,2)+GPT(3)*PM(2,3)+PM(2,4))/SPT4
      SPT3 = (GPT(1)*PM(3,1)+GPT(2)*PM(3,2)+GPT(3)*PM(3,3)+PM(3,4))/SPT4
c....................... flip X and Y if desired........................
      IF (IROT.EQ.1) THEN
        T      = SPT1
        SPT1   =-SPT2
        SPT2   = T
      ENDIF
C-------------------- compute the screen coords ------------------------
      SPT(1) = XCEN + SPT1 * SC_X
      SPT(2) = YCEN + SPT2 * SC_Y
      SPT(3) =        SPT3 * SC_Z
      RETURN
      END
C-----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE CAMERA(EYE,COA,UP, DEYE,FEYE,PM,IOP)
C 
C     sets up the camera transformation matrix
C
C     based on line_of_sight EYE (unit vector)  eg. 1 1 1
C              camera_y_axis UP  (unit vector)  eg  0 1 0
C              centre_of_attention COA          eg 123. 87. -20.
C     DEYE = eye distance, FEYE = eye field of view (eg. 40.)
C     M = returned homogenous 4x4 transformation matrix
C
C     if IOP>0 then use this to re-define EYE to preset values (1-6-92)
C
      REAL EYE(*), COA(*), UP(*)  
      REAL PM(4,4),   RM(4,4), TM(4,4), MP(4,4), MM(4,4), RX(3)

      IF (IOP.gt.0) THEN
        IF (IOP.eq.1) THEN    ! down y-axis
          XME =   0.
          YME =  89.999   ! (to avoid the singularity)
        ELSEIF (IOP.eq.2) THEN    ! 'perspective'
          XME = -30.
          YME =  15.
        ELSEIF (IOP.eq.3) THEN    ! down x-axis
          XME =   0.
          YME =   0.
        ELSEIF (IOP.eq.4) THEN    ! down z-axis
          XME =  90.
          YME =   0.
        ENDIF
        CALL XM2EYE (XME,YME,EYE)      ! convert to a unit vector
      ENDIF

C ----- make EYE,UP unit vectors -----
      CALL UNITV (EYE,3)
      CALL UNITV (UP, 3)

C ----- first appply the translation matrix 'TM' -----
      CALL UNITM (TM,4,4)
      DO I=1,3
        TM(I,4) = - COA(I) 
      ENDDO
C ----- next apply the rotation matrix 'RM' -----
      CALL UNITM (RM,4,4)
      RX(1) = UP(2)*EYE(3) - UP(3)*EYE(2)
      RX(2) = UP(3)*EYE(1) - UP(1)*EYE(3)
      RX(3) = UP(1)*EYE(2) - UP(2)*EYE(1)
      CALL UNITV (RX,3)
      DO I=1,3
        RM(1,I) = RX(I)
      ENDDO
      RM(2,1) = EYE(2)*RX(3) - EYE(3)*RX(2)  
      RM(2,2) = EYE(3)*RX(1) - EYE(1)*RX(3)  
      RM(2,3) = EYE(1)*RX(2) - EYE(2)*RX(1)  
      DO I=1,3
        RM(3,I) = EYE(I)
      ENDDO
C ----- next apply the  perspective matrix 'PM' ------
C .. this is based on the viewing distance DEYE
C ..           and the subtended eye angle FEYE

      D = ATAN( (90.-FEYE)/2. *3.14159265/180.)
      CALL UNITM (MP,4,4)
      MP(3,3)  = 0
      MP(4,3)  = 1./D
      MP(3,4)  = DEYE
      MP(4,4)  = DEYE/D

C     CALL MATMUL(PM,4,TMP,4, M ,4, 4,4,4)
C -- do a 4x4 MATMUL explicitly      
      DO I=1,4
        DO J=1,4
          MM(I,J) = RM(I,1)*TM(1,J) + RM(I,2)*TM(2,J)
     +             +RM(I,3)*TM(3,J) + RM(I,4)*TM(4,J)
        ENDDO
      ENDDO
C -- apply perspective to PM explicitly
      DO I=1,4
        DO J=1,4
          PM(I,J)= MP(I,1)*MM(1,J) + MP(I,2)*MM(2,J)
     +            +MP(I,3)*MM(3,J) + MP(I,4)*MM(4,J)
        ENDDO
      ENDDO
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE XM2EYE(XM,YM,EYE)
C
C     to convert latitude/longitude (XM,YM) to vector EYE (or LIGHT)
C     - and also resets the XM to -180 -> +180,  YM = -180 -> +180
      REAL EYE(*)
      DTR = 3.14159265/180.

      XM = MOD(XM+180.+360.,360.) -180. 
      YM = MOD(YM+180.+360.,360.) -180. 
      EYE(1) =  SIN (XM*DTR) 
      EYE(3) =  COS (XM*DTR)
      EYE(2) =  SIN (YM*DTR)
      EYE(1) = EYE(1) * COS(YM*DTR)
      EYE(3) = EYE(3) * COS(YM*DTR)
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE UNITV (V,N)
C
C     this scales a vector V into a unit vector
C
      REAL V(*)
      VL = 0.
      DO I=1,N
        VL = VL + V(I)**2
      ENDDO
      VL = SQRT(VL)
      DO I=1,N
        V(I) = V(I) / VL
      ENDDO
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE UNITM (A,IA,N)
C
C     makes the top NxN of A into a unit matrix
C
      REAL A(IA,*)
      DO I=1,N
        DO J=1,N
          A(I,J) = 0.
          A(I,I) = 1.
        ENDDO
      ENDDO
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE GOURAUD (SC,ISC,COLI, N_CONTS,IFS,IOP,ITYPE,ICMAT)
C
C     this shades a facet  -- decomposing into triangles
C
      REAL SC(ISC,*)  ,COLI(*)        !- given facet and nodal colours
     +  ,XP(3) ,YP(3) ,COL(3)         !- the abstracted triangles

      INTEGER T8NQ(3,6)
      DATA T8NQ/1,2,8, 7,6,8, 6,8,2, 6,2,4, 2,3,4, 4,5,6/  !- 8nq :-)


      IF (IOP.eq.-1) THEN      !- if 'off' skip
      ELSEIF (IOP.eq.-2) THEN  !- if 'average' : colour as the maean value
        CC = 0.
        DO I=1,IFS
          CC = CC + COL(I)
        ENDDO
        CALL DRAW_POLY (SC,ISC,IFS,INTS(CC/IFS),2)
      ELSE                     !-  ie '0'= 'black' etc. etc.
        NTRI = IFS - 2               !-   # triangles  
C--------------------------- loop the sub-triangles --------------------
        DO IT = 0,IFS-3

          IF (IFS.eq.8) THEN   ! if 8nq then dissect explicitely (9nq?)
            DO I=1,3                  !- 5nq for 14nb's ? <- sub_facet
              II = T8NQ(I,IT+1)       !- etc. ??
              XP (I) =   SC (II,1)
              YP (I) =   SC (II,2)   
              COL(I) = COLI (II)
            ENDDO

          ELSE                    !-- default is triangle at node '1'
            XP (1) =   SC (1,1)
            YP (1) =   SC (1,2)
            COL(1) = COLI (1)
            DO I=2,3
             XP (I)  =   SC (IT+I,1)
             YP (I)  =   SC (IT+I,2)
             COL(I)  = COLI (IT+I)
            ENDDO
          ENDIF  !- of 8nq/others

          CALL GOURAUD_TRI (XP,YP,COL, N_CONTS,IOP,ITYPE,ICMAT) 
        ENDDO    !- next triangle
      ENDIF   !- if not 'skipped' or 'averaged'
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE GOURAUD_TRI (XP,YP,COL, N_CONTS,IOP,ITYPE,ICMAT)
C
C     This contours the given triangle  (filled/lines)
C
      REAL XP(3), YP(3), COL(3)  !- given triangle and nodal values
      REAL X(6),Y(6)             !- the polygon to draw
      INTEGER A,B,C              ! labels of the 3 points

      BIS   (VAL,E,F) = MIN(MAX(0., (VAL-E)/(F-E+t)), 1.)  
      RATIO (VAL,G,H) = G + (H-G) * VAL            !- statement funcs
      T = 1.E-4                     !-   = a 'small' number

      A = 1                         !-  find :  'A' the largest value
      B = 2                         !           'B' the middle value
      C = 3                         !           'C' the smallest value
      IF (COL(A).LT.COL(C) ) CALL ISWAP (C,A)
      IF (COL(B).LT.COL(C) ) CALL ISWAP (B,C)
      IF (COL(A).LT.COL(B) ) CALL ISWAP (A,B)

      N = 0
      I_LO = MAX (INT(COL(C)+t), 0)
      I_HI = MIN (INT(COL(A)+1+t), N_CONTS)

C--------------- loop from the MIN contour to the MAX contour ----------
      DO I= I_LO, I_HI     
        RI = REAL (I)

        IF (I.EQ.INT(COL(B)+1.) ) THEN    !- add middle node if necessary
          N = N + 1
          X(N) = XP(B)
          Y(N) = YP(B)
        ENDIF
c                                   (  be careful if COL(B)=COL(C) !)
        N = N + 1
        IF (I.LT.INT(COL(B)+1.)) THEN            !- edge B -> C
          FAC  = BIS   (RI, COL(C),COL(B))
          X(N) = RATIO (FAC, XP(C), XP(B))
          Y(N) = RATIO (FAC, YP(C), YP(B))
        ELSE                                    !- edge A -> B 
          FAC  = BIS   (RI, COL(B),COL(A))
          X(N) = RATIO (FAC, XP(B), XP(A))
          Y(N) = RATIO (FAC, YP(B), YP(A))
        ENDIF

        N = N + 1                               !- 'long' edge A->C
        FAC  = BIS   (RI, COL(C),COL(A))
        X(N) = RATIO (FAC, XP(C), XP(A))
        Y(N) = RATIO (FAC, YP(C), YP(A))

        IF (ITYPE.EQ.1 .AND. I.GT.I_LO)                  !- not the
     +    CALL GOURAUD_FILL (X,Y,N,I,IOP, ICMAT)           !- 'first'

        IF (ITYPE.EQ.2 .AND. I.GE.COL(C).AND.I.LE.COL(A))  !- not if 
     +    CALL GOURAUD_LINE (X,Y,N,I,IOP, ICMAT)           !- 'clipped off'

        X(1) = X(N)       !- copy last 2 nodes to list base
        Y(1) = Y(N)       !- ready for the next contour band
        X(2) = X(N-1)
        Y(2) = Y(N-1)
        N = 2
      ENDDO
      END
c-----------------------------------------------------------------------
      SUBROUTINE GOURAUD_FILL (X,Y,N,ICOL,IOP, ICMAT)
C
C     This is a sub-part of 'GOURAUD' which simply draws the shapes
C     -- FILLED
C     IOP is the 'style' : black/coloured/zebra etc.
C     ITYPE =1 for filled and =2 for just the lines

      REAL X(*),Y(*)

      J = ICOL/2
      K = ICOL-2*J
      ICOL2 = 2 + MOD (ICOL,14)    !- offset into the colour table


      IF(IOP.eq.0) THEN
        CALL DR_PRIM (X,Y,N,0,2)               !- black
      ELSEIF (IOP.eq.1) THEN
        CALL DR_PRIM (X,Y,N,1,2)               !- white
      ELSEIF (IOP.eq.2) THEN
        CALL DR_PRIM (X,Y,N,ICOL2,2)            !- coloured
      ELSEIF (IOP.eq.3) THEN
        CALL DR_PRIM (X,Y,N,MOD(ICOL2,2),2)     !- zebra evens
      ELSEIF (IOP.eq.4) THEN
        CALL DR_PRIM (X,Y,N,MOD(ICOL2+1,2),2)   !- zebra odds
      ELSEIF (IOP.eq.5) THEN
        IF (K.EQ.0) CALL DR_PRIM (X,Y,N,J+1,2) !- stripe with invis
      ELSEIF (IOP.eq.6) THEN
        IF (K.EQ.0) CALL DR_PRIM (X,Y,N,J+1,2)
        IF (K.EQ.1) CALL DR_PRIM (X,Y,N,0,2)   !- stripe with 0
      ELSEIF (IOP.eq.7) THEN
        IF (K.EQ.0) CALL DR_PRIM (X,Y,N,J+1,2)
        IF (K.EQ.1) CALL DR_PRIM (X,Y,N,1,2)   !- stripe with 1
      ELSEIF (IOP.eq.8) THEN
        CALL DR_PRIM (X,Y,N,ICMAT,2)           !- mat colour
      ELSEIF (IOP.eq.9) THEN
        IF (K.EQ.0) CALL DR_PRIM (X,Y,N,ICMAT,2) !- zebra with mat col
        IF (K.EQ.1) CALL DR_PRIM (X,Y,N,1,2) 
      ELSEIF (IOP.eq.10) THEN
        ICOL3= (15-ICOL)*16 +ICMAT+1               
        CALL DR_PRIM (X,Y,N,ICOL3,2)             !- shade material ?
      ENDIF
      END
c-----------------------------------------------------------------------
      SUBROUTINE GOURAUD_LINE (X,Y,N,ICOL,IOP, ICMAT)
C
C     This is a sub-part of 'GOURAUD' which simply draws the shapes
C     -- LINES
C     IOP is the 'style' : black/coloured/zebra etc.

      REAL X(*),Y(*)

      J = ICOL/2
      K = ICOL-2*J
      IF (IOP.eq.0) THEN
         CALL DR_PRIM (X(N-1),Y(N-1),2,0,3)              !- black
       ELSEIF (IOP.eq.1) THEN                    
         CALL DR_PRIM (X(N-1),Y(N-1),2,1,3)              !- white
       ELSEIF (IOP.eq.2) THEN                        
         CALL DR_PRIM (X(N-1),Y(N-1),2,ICOL,3)           !- coloured
       ELSEIF (IOP.eq.3) THEN
         CALL DR_PRIM (X(N-1),Y(N-1),2,MOD(ICOL,2),3)    !- zebra:1
       ELSEIF (IOP.eq.4) THEN
         CALL DR_PRIM (X(N-1),Y(N-1),2,MOD(ICOL+1,2),3)  !- zebra:2
       ENDIF
      END
C-----------------------------------------------------------------------
      SUBROUTINE FILL_POLY (XI,YI,N,ICOL)
C                                                            DJK 19-3-92
C     this colours-in a polygon (integer co-ords)            
C     ... with co-incident vertex elimination :-)
C
      INTEGER    XI(*), YI(*),N,ICOL
      INTEGER*2  X(30), Y(30),IFAIL,NN,IPOLY,ICOL2
      IF (ICOL.NE.-1) then  ! skip 'invisible colour = -1'
C--------- patch to remove 'zero-height' polygons 14-07-92 -------------
        ICOL2 = ICOL
        NN = N
        iymin = xi(1)
        iymax = yi(1)
        DO I=1,N
          X(I) = XI(I)
          Y(I) = YI(I)
          iymin = min(iymin,y(i))
          iymax = max(iymax,y(i))
        ENDDO
        IF (iymin.ne.iymax) THEN
C--------------------------------------------------------------------

        CALL CREATE_POLYGON@      (X,Y,NN,IPOLY,  IFAIL)
        CALL DOSERR@(IFAIL)
        CALL FILL  _POLYGON@            (IPOLY,ICOL2,IFAIL)
        CALL DOSERR@(IFAIL)
        CALL DELETE_POLYGON_DEFINITION@ (IPOLY,  IFAIL)
        CALL DOSERR@(IFAIL)

c.. debug..info
c         write(19,'(50(''-''))')
c        write(19,'(5(''  <>'',2i4))') (x(i),y(i),i=1,nn),icol

        ENDIF
      ENDIF
      RETURN
      END
C----------------------------------------------------------------------
        SUBROUTINE DRAW_POLY (SC,ISC,IFS,ICOL,IOP)
C
C       this calls FILL_POLY to draw the a facet
C       IOP = 1 for edges, = 2 for filled
C       --->if ICOL .lt. 0  then skipped
C
C       2-5-92 CCOLS() contains the colour options (filled only ?)

C...... it would be much better to do this in GOURAUD !
      COMMON /DEST/IOUT_D     !-- eg. 5=postscript

c     COMMON /CCOLS/CCOL    !-- passed from the main program
C     INTEGER CCOL(6)

c  ccol(1) = (14)  # colour contour bands (not really used here)
c  ccol(2) = (1)   colour 'from' ie. lower limit of 'wrap'
c  ccol(3) = (14)  colour 'to'       upper limit of 'wrap' 
c  ccol(4) = (1)   offset into 'real' colours. ie. all +1
c  ccol(5) = (-1)  'infill' code = -1 =none
c                                   0 = alternate 'background' bands
c                                   1 = alternate 'invisble' bands
c  ccol(6) = (-1)  'wrap' code = -1 = clip to the max/min values
c                                 0 = clip to 'background
c                                 1 = clip to 'invisible'
c                                 2 = wrap around (so 15->1,16->2, etc.)
c                                 3 = oscilate    (so 13->1,16->12, etc.)

      REAL SC(ISC,*)    ! x,y,(z) coords of the line/polygon
      INTEGER*2 XPI(20), YPI(20), IFSI, IFSP1,ICOLI

      IF (ICOL.lt.0) RETURN      ! ie. if color = -1 exit !
      ICOLI = ICOL
      IFSI  = IFS
      DO I=1,IFS         ! copy real coords to integer x,y lists 
        XPI(I) = NINT (SC(I,1))
        YPI(I) = NINT (SC(I,2))
      ENDDO

      IF (IOP.EQ.3) THEN                      ! simple line draw 
        IF (IOUT_D.EQ.5) THEN
          CALL DR_PS (XPI,YPI,IFSI,ICOLI, 1)
        ELSE
          CALL POLYLINE@ (XPI,YPI, IFSI, ICOLI)
        ENDIF
      ELSEIF (IOP.EQ.1) THEN          ! closed polygon line draw 

        IFSP1 = IFS+1
        XPI (IFSP1) = XPI (1)
        YPI (IFSP1) = YPI (1)
        IF (IOUT_D.EQ.5) THEN
          CALL DR_PS (XPI,YPI,IFSI,ICOLI, 2)
        ELSE
          CALL POLYLINE@       (XPI,YPI, IFSP1, ICOLI)
        ENDIF
      ELSEIF (IOP.EQ.2) THEN
        IF (IOUT_D.EQ.5) THEN
          CALL DR_PS (XPI,YPI,IFSI,ICOLI, 5)
        ELSE
          CALL FILL_POLY (XPI,YPI,IFSI,ICOLI)
        ENDIF
      ELSE
        PRINT*,'IOP=',IOP,' unknown'
      ENDIF
      END

C----------------------------------------------------------------------
C----------------------------------------------------------------------
        SUBROUTINE DR_PRIM (X,Y,N,ICOL,IOP)
C
C       This is the Main primitive drawing interface
C
C       IOP = 1 for closed edges, = 2 for filled, 3=unclosed edge ,etc.
C
C       --->if ICOL .lt. 0  then skipped

      COMMON /DEST/IOUT_D     !-- eg. 5=postscript

c     COMMON /CCOLS/CCOL    !-- passed from the main program
C     INTEGER CCOL(6)

c  ccol(1) = (14)  # colour contour bands (not really used here)
c  ccol(2) = (1)   colour 'from' ie. lower limit of 'wrap'
c  ccol(3) = (14)  colour 'to'       upper limit of 'wrap' 
c  ccol(4) = (1)   offset into 'real' colours. ie. all +1
c  ccol(5) = (-1)  'infill' code = -1 =none
c                                   0 = alternate 'background' bands
c                                   1 = alternate 'invisble' bands
c  ccol(6) = (-1)  'wrap' code = -1 = clip to the max/min values
c                                 0 = clip to 'background
c                                 1 = clip to 'invisible'
c                                 2 = wrap around (so 15->1,16->2, etc.)
c                                 3 = oscilate    (so 13->1,16->12, etc.)

      REAL X(*),Y(*)        !--- x,y,(z) coords of the line/polygon
      INTEGER*2 XPI(20), YPI(20), IFSI, IFSP1,ICOLI

      IF (ICOL.lt.0) RETURN      ! ie. if color = -1 exit !
      ICOLI = ICOL
      IFSI  = N
      DO I=1,N                 ! copy real coords to integer*2 x,y lists 
        XPI(I) = NINT (X(I))   !-- (no need for PostScript :-)
        YPI(I) = NINT (Y(I))
      ENDDO

      IF (IOP.EQ.3) THEN       !----------simple line draw -------------
        IF (IOUT_D.EQ.5) THEN
          CALL DR_PS (XPI,YPI,N,ICOL, 1)
        ELSE
          CALL POLYLINE@ (XPI,YPI, IFSI, ICOLI)
        ENDIF
      ELSEIF (IOP.EQ.1) THEN   !--------closed polygon line draw--------
        IF (IOUT_D.EQ.5) THEN
          CALL DR_PS (XPI,YPI,N,ICOL, 2)
        ELSE
          IFSP1 = IFSI+1
          XPI (IFSP1) = XPI (1)
          YPI (IFSP1) = YPI (1)
          CALL POLYLINE@       (XPI,YPI, IFSP1, ICOLI)
        ENDIF
      ELSEIF (IOP.EQ.2) THEN
        IF (IOUT_D.EQ.5) THEN
          CALL DR_PS (XPI,YPI,N,ICOL, 5)
        ELSE
          CALL FILL_POLY (XPI,YPI,IFSI,ICOLI)
        ENDIF
      ELSE
        PRINT*,'IOP=',IOP,' unknown'
      ENDIF
      END

C----------------------------------------------------------------------
      SUBROUTINE DR_ARROW (pt_fr,pt_to,iangle,iscale,icol)
C
C      to draw an arrow vector (from L2.for)    c. Dan Kidger
C       eg. angle=10. (deg) , scale = 10 (%)
C
c.... nicer as a sub-function of DR_PRIM ..extends X and Y for the head! 

c     INTEGER*2 IXT,IYT 
      REAL PT_FR(*),PT_TO(*) ,  x(6),y(6)

c      IXT = PT_TO(1)    
c      IYT = PT_TO(2)    

      XL=pt_to(1) - pt_fr(1)
      YL=pt_to(2) - pt_fr(2)

c     CALL DRAW_LINE@ (INTS(PT_FR(1)),INTS(PT_FR(2)), IXT,IYT,ICOL)

      IF (ISCALE.GT.0) THEN           ! skip if no arrow head
        A  =  3.14159295/180  * IANGLE 
        SA = -SIN(A) * ISCALE/100.
        CA = -COS(A) * ISCALE/100.

c        CALL DRAW_LINE@ (IXT,IYT, INTS(IXT + XL*CA+YL*SA),
c     +                            INTS(IYT - XL*SA+YL*CA), ICOL)
c        CALL DRAW_LINE@ (IXT,IYT, INTS(IXT + XL*CA-YL*SA),
c     +                            INTS(IYT + XL*SA+YL*CA), ICOL)

        x(1) = pt_fr(1)    !- base
        y(1) = pt_fr(2)
        x(2) = pt_to(1)    !- top
        y(2) = pt_to(2)

        x(3) = x(2) + XL*CA + YL*SA     !- LH tip
        y(3) = y(2) - XL*SA + YL*CA
        x(4) = x(2)                 !- top again
        y(4) = y(2)
        x(5) = x(2) + XL*CA - YL*SA     !- RH tip
        y(5) = y(2) + XL*SA + YL*CA

        CALL DR_PRIM (X,Y,5,ICOL,3)    !----- edges

      ENDIF
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE GRATIC (IRESX,IRESY,NX,NY,ICOL)
C
C     draws a graticule of RESX by RESY boxes
C .. should really draw within the current window ??
C

      INTEGER*2 IRESX,IRESY
      REAL X(2),Y(2)
      JRESX = IRESX-1         !-- so that the far side is visible
      JRESY = IRESY-1         !--  "    "   " "   "
      DO I=0,NX                        !--- horizontal rulings
        X(1) = 0.
        X(2) = JRESX
        Y(1) = REAL(JRESY)*I / REAL(NX)
        Y(2) = Y(1)
        CALL DR_PRIM (X,Y,2,ICOL,3)    !----- edges
      ENDDO

      DO I=0,NY                         !--- vertical rulings       
        X(1) = REAL(JRESX)*I / REAL(NY)
        X(2) = X(1)
        Y(1) = 0.
        Y(2) = JRESY
        CALL DR_PRIM (X,Y,2,ICOL,3)    !----- edges
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE READ_NFF (U3,IDT,NN,NDIM,NDIME,ITYPE,
     +     GCOORD,IGCRD,INF,NEL,NUMS,INUMS,MEL,NLDS,CV,PAL)   
c
c     this reads in data from the NFF= 'Neutral FILE Format'
c     currently only polygons are read
c     -- however NFF also contains the viewing params, colours, etc.
c


      REAL GCOORD(IGCRD,INF)
      INTEGER NUMS(INUMS,MEL),U3
      INTEGER CV(*),PAL(3,0:255)
      CHARACTER code*2,t*1, line*80

      NDIM = 3
      NDIME= 2
      NLDS = 0
      ITYPE= 9

      NN   = 0
      NEL  = 0
      IMAT = 1
 1111 CONTINUE
      READ (U3,'(a)',ERR=999) CODE
      IF (CODE.eq.'f ') THEN   ! a material colour
       BACKSPACE (U3)
        IMAT = IMAT + 1
       READ (U3,'(A)') LINE
       LINE(1:2) ='  '   !-- zap the 'key-letter'
       READ (LINE,*) R,G,B
       PAL(1,IMAT) = R *256 
       PAL(2,IMAT) = G *256 
       PAL(3,IMAT) = B *256 
       DO J=1,3
         IF (PAL(J,IMAT).EQ.256) PAL(J,IMAT) = 255
       ENDDO

C-----------  a polygon / polygon-with-normals ('p','pp') --------------
      ELSEIF (CODE.eq.'p '.or.CODE.eq.'pp') THEN 
      BACKSPACE (U3)
          IF (CODE.eq.'p ') READ (U3,'( a,i2)')   t,NOD
          IF (CODE.eq.'pp') READ (U3,'(2a,i2)') t,t,NOD
        NEL = NEL+1
        NUMS(NEL,1) = NOD
        DO I=1,NOD
          NN=NN+1
          READ (U3,*) (GCOORD(J,NN),J=1,3)
          NUMS (NEL,I+1) = NN
        ENDDO
        NUMS (NEL,NOD+2) = IMAT
      ENDIF
      GOTO 1111
  999 CONTINUE
      CALL SET_ALL_DACS ()

      RETURN
      END

C------------------- Dan's .PL format ----------------------------------
C------------------ 'Object File Format' -------------------------------
C---------------------- 'David Ho's Format -----------------------------
      SUBROUTINE READ_NODES (IO,IDT,NN,NDIM,GCOORD,IGCRD,INF,NEL)
c
c  this reads in the nodal co-ordinates for various file formats
c     IDT=1 (my .PL), =2 (OFF) , =3 (D.Ho's)

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
        READ (IO,*) NN, NEL, NEDGES    ! NEL needs storing !
        READ (IO,*) ((GCOORD(J,I),J=1,3), I=1,NN)
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
      SUBROUTINE READ_ELEMS (IO,IDT,NEL,NDIME,ITYPE,NUMS,INUMS,MEL)
c
c     this reads in the 'elements' into NUMS for different file formats
c     IDT=1 (my .PL), =2 (OFF) , =3 (D.Ho's)
c
      INTEGER NUMS(INUMS,MEL)
C------------------- Dan's .PL format ----------------------------------
      IF (IDT.eq.1) THEN
        ITYPE = 1            ! default = 1
        READ (IO,*) NEL_
        NEL = 0
        DO I = 1,NEL_
          READ(IO,*,ERR=99) II,NUMS(II,1), (NUMS(II,J),J=2,2+NUMS(II,1))
          NEL = MAX (II,NEL)
        ENDDO
  99    CONTINUE
C------------------ 'Object File Format' -------------------------------
      ELSEIF (IDT.eq.2) THEN 
        NDIME = 2       !      with 2D 'polygon' elements
        ITYPE = 9       !      a 'polygonal' element
        DO II=1,NEL
          READ(IO,*) NUMS(II,1), (NUMS(II,J),J=2,1+NUMS(II,1) )
          NUMS (II,NUMS(II,1)+2) = NUMS(II,1)  ! all mat.types= NEN! :-)
        ENDDO
C---------------------- 'David Ho's Format -----------------------------
      ELSEIF (IDT.eq.3) THEN
        ITYPE = 1            ! default = 1
        READ (IO,*) NEL
        DO I = 1,NEL
          NEN = 14
          READ(IO,*) II, (NUMS(II,J+1),J=1,NEN)
          NUMS(II,1)     = NEN   ! #nodes per element
c         NUMS(II,2+NEN) = 1     ! material types
        ENDDO
        READ(IO,*) (NUMS(II,2+NEN),II=1,NEL)
      ELSE
        PRINT*,'*** WARNING: data_type=',IDT,' unknown (READ_ELEMS)'
        NN = 0
      ENDIF
      RETURN
      END
C ----------------------------------------------------------------------
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
      DO NLDS=1,99
        IBASE = NN * (NLDS-1) 
        IF (IBASE+NN.GT.IGDISPS)THEN
          PRINT*,'*** TOO Many Load steps.. only',NLDS-1,' used'
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
      ELSEIF (IDT.eq.2.or.IDT.eq.4) THEN
        print*,'** ignoring disps for ''OFF'' format'
        NLDS = 0
C---------------------- 'David Ho's Format -----------------------------
      ELSEIF (IDT.eq.3) THEN
        NODOF = 3      !--- always in 3D
        DO NLDS=1,99
        IBASE = NN * (NLDS-1) 
        IF (IBASE+NN.GT.IGDISPS)THEN
          PRINT*,'*** TOO Many Load steps.. only',NLDS-1,' used'
          GOTO 99
        ENDIF
        READ(IO,'(A)',ERR=99) LINE
        WRITE(*,'(a,i3,a,a)') 'load step >',nlds,' <', LINE(1:50)
        DO I=1,NN
          DO J=1,3
            GDISPS(J,IBASE+I) = 123.e-30  ! set the disps to zero
          ENDDO
        ENDDO
        DO I=1,NN
          READ (IO,*) II,(GDISPS(J,IBASE+II),J=1,NODOF)
        ENDDO
      ENDDO
   99 NLDS=NLDS-1
      print*, ' Total No. of Load steps = ',NLDS,'           '
      LD = NLDS
      ELSE
        PRINT*,'*** WARNING: data_type=',IDT,' unknown (READ_ELEMS)'
        NLDS = 0
      ENDIF
      RETURN
      END
C ----------------------------------------------------------------------
C----------------------------------------------------------------------
      SUBROUTINE RANGE (GCOORD,IGCRD, NN,NDIM, COA,DIAG)
C
C     this calculates range of the data in GCOORD
C      also the diagonal of the bounding box DIAG
C     (NN+1)  x,y,z minima,    
C     (NN+2)  x,y,z lengths,   
C     (NN+3)  x,y,z (centroid)
C     COA = 'middle' of the range of co-ordinats in GCOORD
C
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
      DIAG = SQRT (( GCOORD(1,NN+1)-GCOORD(1,NN+2) )**2
     +        +    ( GCOORD(2,NN+1)-GCOORD(2,NN+2) )**2
     +        +    ( GCOORD(3,NN+1)-GCOORD(3,NN+2) )**2 )

C--------- set the CENTRE_OF_ATTENTION to the 'middle' of the object ---
      DO I=1,3
        COA(I)= (GCOORD(I,NN+1) + GCOORD(I,NN+2)) /2.
      ENDDO
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE FSTRIP2(NEL,NUMS,INUMS,FACETS,PL2,IFCETS,NFCETS,NDIME)
C
C     THIS turns NUMS into FACETS-a list of element face pointers
C     .. the sequel to FSTRIP.. but all FACETS are recorded
C
      INTEGER NUMS(INUMS,*), FACETS(3,IFCETS),PL2(*),FS1(19),FS2(19)
      ITYPE= 1   !--should be an arguement (eg. ITYPE=9)

C------------- patch for speed-up if NDIME=2d --------------------------
      IF (NDIME.LT.3) THEN
        NFCETS = NEL
        DO I=1,NFCETS
          FACETS(1,I) = I   !--> = the element number
          FACETS(2,I) = 1   !--> only one face
          FACETS(3,I) = 0   !--> all facets on the 'outside'
        ENDDO
        RETURN
      ENDIF
C------------------- find the MIN node # of each facet -----------------
      IFC = 0
      DO IEL1=1,NEL
        NOD1  = NUMS(IEL1,1)
        IMAT1 = NUMS(IEL1,NOD1+2)
C.. in the next *only* NN_F is returned for: option=1
        CALL GFACE (NOD1,NDIME,ITYPE,NFACES,FS1,NN_F1,NN_FT1, 1 )

        DO IFACE1=1,NFACES         ! loop this element's faces
          CALL GFACE (NOD1,NDIME,ITYPE,IFACE1,FS1,NN_F1,NN_FT1, 2 ) 
          N_MIN= NUMS (IEL1,1+FS1(1))
          N_MAX= NUMS (IEL1,1+FS1(1))
          DO I=2,NN_F1
            N_MIN= MIN (N_MIN,NUMS (IEL1,1+FS1(I)))
            N_MAX= MAX (N_MAX,NUMS (IEL1,1+FS1(I)))
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
        ITMP = PL2(IFC1)
        DO IFC2 = IFC1+1,NFCETS
          IF (PL2(IFC2).NE.ITMP) GOTO 1       !--- skip  (no touch)
          IF (FACETS(3,IFC2).NE.FACETS(3,IFC1)) GOTO 1   !--- skip too
          IF (FACETS(3,IFC2).GT.0) GOTO 1   !--- skip (already done)

C         ----- a possible 'hit' so test all -----------------
            IFACE1= FACETS(2,IFC1)
            IFACE2= FACETS(2,IFC2)
            IEL1  = FACETS(1,IFC1)
            IEL2  = FACETS(1,IFC2)
            NOD1  = NUMS(IEL1,1)
            NOD2  = NUMS(IEL2,1)
            CALL GFACE (NOD1,NDIME,ITYPE,IFACE1,FS1,NN_F1,NN_FT1, 2 )
            CALL GFACE (NOD2,NDIME,ITYPE,IFACE2,FS2,NN_F2,NN_FT2, 2 )

            IF (NN_FT1.NE.NN_FT2) GOTO 1   ! skip if incompatible
            IBOT = NUMS(IEL1,1+FS1(1))
            DO I2=1,NN_F1        ! loop to try and find the 'base' node
              IF (NUMS(IEL2,1+FS2(I2)) .EQ. IBOT ) GOTO 3
            ENDDO
            GOTO 1    ! no match so skip to try its next facet

    3     IBASE = I2      ! so node I2==IBOT.. the 'base node'

          DO I1=2,NN_F1
            I2 = 1+MOD (NN_F1+IBASE -I1, NN_F1)     
            IF (NUMS(IEL1,1+FS1(I1)).ne.NUMS(IEL2,1+FS2(I2)) ) GOTO 1
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

      SUBROUTINE GFACE (NEN,NODOF,ITYPE,IFACE,FACE,NN_F,NN_FT,IOP)
C
C     this returns a list of nodes making up a face of an element
C       ( of NOD nodes in NODOF dimensions, type = ITPYE )
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

C--- increment the next line whenever a new element is added :-)
      PARAMETER    (MAX_INFO=227+35+10)  !14nb as 5nq's
C     PARAMETER    (MAX_INFO=227+35+10-6)  !14nb as 4nq's
      INTEGER INFO (MAX_INFO)  , FACE(*)
      SAVE NEN_OLD, NODOF_OLD, IFACE_OLD, IBASE
      DATA NEN_OLD, NODOF_OLD, IFACE_OLD /0,0,0/

      DATA INFO/ !------------------------------------------------------
C------>  = NEN,NODOF,#faces,#nodes per face,#mid-face nodes,
C          /nodes (-ve if 'not-a-corner')/ 
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
      endif
        
      IF (NEN_OLD.NE.NEN.AND.NODOF_OLD.NE.NODOF) THEN

        NEN_OLD   = NEN
        NODOF_OLD = NODOF
        ibase = 1

    2   continue
c        print*,'trying,nen=,nodof=',info(ibase),info(ibase+1)
        if (nen.ne.INFO(ibase).or.nodof.ne.INFO(ibase+1) ) then  
          ibase = ibase + info(ibase+2) * info(ibase+3) + 5
          if (ibase.lt.MAX_INFO) goto 2
        else
          goto 1     ! found this element
        endif
c         ...... unfound so warn user and exit .........
          print*,'*** WARNING, this element not found'//
     +         'NEN=',NEN,' nodof=',nodof
        RETURN
      ENDIF

C-----------------------------------------------------------------------
    1 CONTINUE  ! OK .. extract this facets' info   
      if (iop.eq.1) then
        iface = info(ibase+2) 
        return
      endif
c.... hmm I 'may' modify FACE() so must always create it !
c     if (IFACE.eq.IFACE_old) RETURN     ! no work at all to do !
      IFACE_OLD = abs (IFACE)
      ibase_f = ibase + 4 + (abs(iface)-1) * info(ibase+3)
      nn_ft =         info(ibase+3)
      nn_f  = nn_ft - info(ibase+4)
      IF(IFACE.gt.0) THEN           ! a 'clockwise' face
        DO I=1,NN_FT
          FACE(I) = ABS(INFO(IBASE_F+I))   ! .. just use ABS for now 
        ENDDO
      ELSE                          ! an 'anticlockwise' face (mirrored)
        FACE(1) = ABS(INFO(IBASE_F+1))
        DO I=2,NN_F
          FACE(NN_F+2-I) = ABS(INFO(IBASE_F+I)) 
        ENDDO
        DO I=NN_F+1,NN_FT
          FACE(I) = ABS(INFO(IBASE_F+I))
        ENDDO
      ENDIF
      RETURN
      END
C----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE SELECT_PS_PRINTER (FILE,IRES_X,IRES_Y)
C   
C     This opens a file for post-script output on unit 55
C       
      CHARACTER FILE*(*)
c    +,TIME@(8),FDATE@(20)
c     EXTERNAL TIME@,FDATE@
      OPEN (55,FILE=FILE)
      WRITE(55,'(A)')
     + '%!Adobe-2.0 EPSF'
     +,'%%Creator: Danplot Finite Element Postprocessor'
c    +,'%%Creationdate '//TIME@()//' '//FDATE@()
     +,'%%Title: Unknown'
     +,'%%BoundingBox 0 0 576 792'

     +,' ' 
     +,'/fp { newpath 3 1 roll moveto {lineto} repeat'//
     +    ' closepath gsave setrgbcolor fill grestore } def' 
     +,'/dl { newpath 3 1 roll moveto {lineto} repeat'//
     +    ' gsave setrgbcolor stroke grestore } def' 
     +,'/dp { newpath 3 1 roll moveto {lineto} repeat'//
     +    ' closepath gsave setrgbcolor stroke grestore } def' 
     +,' ' 
c    +,'gsave clippath pathbbox grestore'  !-----use Landscape mode-----
c    +,'4 dict begin'
c    +,'/ury exch def /urx exch def /lly exch def /llx exch def'
c    +,'90 rotate      ury urx sub llx ury add neg translate'
c    +,'end'
c    +,' ' 
c    +,'newpath llx lly moveto urx lly lineto urx ury lineto llx ury '//
c    + 'moveto closepath gsave 20 setlinewidth stroke grestore'
     +,' ' 
     +,' 72 300 div dup scale '      !--- units in 300ths inch
     +,' ' 


      IRES_X =  8 * 300
      IRES_Y =  8 * 300
c     IRES_Y = 11 * 300

      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE CLOSE_PS_PRINTER ()
C   
C     This closes the post-script output on unit 55
C       
      WRITE(55,'(A)')
     +       '%%   finally show the page'  
     +      ,'showpage'
      CLOSE(55)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE DR_PS (X,Y,N,ICOL, IOP)
C
C     This writes PostScript Output to unit 55 for lines,fill_areas etc.
C     IOP = 1 un-closed line
C     IOP = 2    closed line = polygon
C     IOP = 5 filled polygon
C
      INTEGER X(*),Y(*),N, ICOL,IOP
      COMMON /PALETTE/PAL
      INTEGER PAL(3,0:255)
      CHARACTER FORMAT*20,ACTION*5

      IF (IOP.EQ.1) ACTION = ' dl'
      IF (IOP.EQ.2) ACTION = ' dp'
      IF (IOP.EQ.5) ACTION = ' fp'

      WRITE(FORMAT,'(A,I2,A)') '(3f6.3,',2*N,'I5,I3,A)' 
      WRITE(55,FORMAT) (PAL(J,ICOL)/256.,J=1,3)
     +                ,(X(I),Y(I),I=N,1,-1), N-1, ACTION
      RETURN
      END
C-----------------------------------------------------------------------
C----------------------------------------------------------------------
      include 'shape.for'
C----------------------------------------------------------------------
C-----------------------------------------------------------------------
            SUBROUTINE INVAR(STRESS,SIGM,DSBAR,THETA)
C
C      THIS SUBROUTINE FORMS THE STRESS INVARIANTS (2-D)
C
      REAL STRESS(*)
      SX=STRESS(1)
      SY=STRESS(2)
      TXY=STRESS(3)
      SZ=STRESS(4)
      SIGM=(SX+SY+SZ)/3.
      DSBAR=SQRT((SX-SY)**2+(SY-SZ)**2+(SZ-SX)**2+6.*TXY**2)/SQRT(2.)
      IF(DSBAR.EQ.0.)THEN
      THETA=0.
      ELSE
      DX=(2.*SX-SY-SZ)/3.
      DY=(2.*SY-SZ-SX)/3.
      DZ=(2.*SZ-SX-SY)/3.
      XJ3=DX*DY*DZ-DZ*TXY**2
      SINE=-13.5*XJ3/DSBAR**3
      IF(SINE.GT.1.)SINE=1.
      IF(SINE.LT.-1.)SINE=-1.
      THETA=ASIN(SINE)/3.
      END IF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE TWOBY2(JAC,IJAC,JAC1,IJAC1,DET)
C
C      THIS SUBROUTINE FORMS THE INVERSE OF A 2 BY 2 MATRIX
C
      REAL JAC(IJAC,*),JAC1(IJAC1,*)
      DET=JAC(1,1)*JAC(2,2)-JAC(1,2)*JAC(2,1)
      JAC1(1,1)=JAC(2,2)
      JAC1(1,2)=-JAC(1,2)
      JAC1(2,1)=-JAC(2,1)
      JAC1(2,2)=JAC(1,1)
      DO 1 K=1,2
      DO 1 L=1,2
    1 JAC1(K,L)=JAC1(K,L)/DET
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE MVMULT(M,IM,V,K,L,Y)
C
C      THIS SUBROUTINE MULTIPLIES A MATRIX BY A VECTOR
C
      REAL M(IM,*),V(*),Y(*)
      DO 1 I=1,K
      X=0.
      DO 2 J=1,L
    2 X=X+M(I,J)*V(J)
      Y(I)=X
    1 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE MATMUL(A,IA,B,IB,C,IC,L,M,N)
C
C      THIS SUBROUTINE FORMS THE PRODUCT OF TWO MATRICES
C
      REAL A(IA,*),B(IB,*),C(IC,*)
      DO 1 I=1,L
      DO 1 J=1,N
      X=0.0
      DO 2 K=1,M
    2 X=X+A(I,K)*B(K,J)
      C(I,J)=X
    1 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FMDEPS(DEE,IDEE,E,V)
C
C      THIS SUBROUTINE FORMS THE ELASTIC PLANE STRAIN
C      STRESS/STRAIN MATRIX
C
      REAL DEE(IDEE,*)
      V1=1.-V
      C=E/((1.+V)*(1.-2.*V))
      DEE(1,1)=V1*C
      DEE(2,2)=V1*C
      DEE(3,3)=.5*C*(1.-2.*V)
      DEE(1,2)=V*C
      DEE(2,1)=V*C
      DEE(1,3)=0.
      DEE(3,1)=0.
      DEE(2,3)=0.
      DEE(3,2)=0.
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      

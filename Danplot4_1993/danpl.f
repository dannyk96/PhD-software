c     OPTIONS (DREAL,INTL,FULLCHECK,UNDEF)
      OPTIONS (fullcheck,undef)

C--------------------------------------------------------------------
      PROGRAM     P L O T T E R _ 7 c
C--------------------------------------------------------------------

      INTEGER*4 MNN,MEL,IGDISPS,IBASE, IB1,IB2         !- yuk use INTL !

      PARAMETER ( MNN = 20 000            ! Max # of nodes
     +           ,MEL = 20 000            ! Max # of elements
     +       ,IGDISPS = 200 000           ! Max # of displaced nodes
     +        ,IFCETS = 30 000       )    ! Max # of facets

C     this is Dan Kidger's general purpose plotting package

c.... NEW history:
c
c  18- 6-93 : move command-line parsing into the 'main-block'
c   5- 6-93 : * * New menus and Makefiles * *

C--------------------- revision history ------------------------------

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
C  30-10-91 on-screen time per frame  '"'.
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
C  30- 4-92  (push eye posn, light posn.), picture # to .() ?
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
C    1-4-93 reversed the default colours to black on white (and back!)
C   10-5-93 PL9.  attempt to add *KEYWORDS ('-k' on command-line)
c     ...... now cut into multiple source files ! ......
C    2-8-93 Made all large-array sizes INTEGER*4
C           reversed NUMS subscripts,  FACETS now (7,*)
C    6-8-93 Added logo and date to PS output, added ARCs to PS

C   1-11-93? *Regridding* added :-)
C   9-11-93 'V' to put mean nodal vector normals into GDISPS
C  16-11-93 moved FSTRIP etc. to drawing section (=streamline data reading)
C  16-11-93 changed all NUMS to the new style.
C  16-11-94 hence extend to parse keyword-files from the command-line
C  22-11-94 GCOORD renamed as GC

C-------- INF= Max nodes, MEL= Max Elements, MN,MDF = max elem = 20nb
C-------  IFCETS= Max. facets



      PARAMETER (   MN = 20         !- \ Max element size
     +             ,MDF = 3         !- /
     +           ,M_CV = 60         !- max # of CV entries
     +        ,MIMAGES = 200 )      !- max # of stored images


      PARAMETER (INUMS = MN+3       !- # data iterms per element
     +          ,IGC = 3 )        !- max # coords per node


c----------------------- The *Big* Arrays ------------------------------
c.. hmm IGDISPS is only used to define and in READ_LOADS

      REAL    GC (IGC,MNN)     !- nodal coordinates
     +       ,GDISPS (MDF,IGDISPS)   !- nodal disps (every loadcase)
     +       ,SCOORD (3,MNN)         !- nodal 'screen' coords (for speed?)

      INTEGER   NUMS (INUMS,MEL)     !- element steering  (now reversed)
     +       ,FACETS (7,IFCETS)      !- table of ALL the 3D faces (rev)

      REAL*4   ZBUF (IFCETS)       !- the 'depths' of the facets
      INTEGER*4 PL1 (IFCETS)       !- RSORTed order of facets  
      INTEGER   PL2 (IFCETS)       !- order for drawing FACETs

c---------------------------- other arrays -----------------------------
      REAL    X (10) ,Y (10)         !- lines.. eg. window-frame
     +       ,XTRA_L (4,30)          !- for 'added' lines on the drawing

C---------- the following need their subscripts reversing ! ------------
c.... if reversed then can pass 'slices' to subroutines

      PARAMETER (ISC = MN)     !- max number of nodes per facet
      REAL     C_F (ISC,3)         !- model-cords of a facet
     +       ,SC_F (ISC,3)         !- screen-coords of a sub-facet
     +      ,SC_SF (ISC,3)         !- screen-coords of a sub-facet
     +       ,DI_F (ISC,3)         !- the disps of a facet
     +      ,LC_SF (ISC,3)         !- local coords of a sub-facet
     +       ,ACOL (ISC)           !- colours at the nodes ob a s-f 
     +        ,VEC (3)             !- general vector
     +        ,C_R (6)             !- contour limits actual/chosen


c... the next is really only 1 part of the table of information needed
c... by the device drivers.. also image x,y size, etc.
      COMMON /DEST/IDEST,IOUT_D    !- eg. 5=postscript
      COMMON /CV/CV                !- why in common ?
      INTEGER   CV (M_CV)          !- **all the colour 'tokens' **

c     INTEGER CCOL(6)              !- color controls (obs)

      INTEGER
     +         U0,U1,U3            !- file unit numbers
     +       ,IOP_ANIM (10)        !- animation 'tokens' (re-do in CV?)
     +            ,NUM (INUMS)     !- nodes of an element    
     +             ,FS (ISC)       !- node numbers of this facet
     +        ,IC_MATS (20)        !- colours for each material type (obs)

      LOGICAL 
     +     KEY_WAITING@            !- keyboard has been pressed (func)
     +            ,EXIT            !- to quit image-drawing
     +       ,L_OVERLAP            !- if drawn FACETS overlap the menu
     +   ,CMN_EXHAUSTED            !- any tokens left on the command-line ?
     +   ,KWF_EXHAUSTED            !- any keywords left in the input file ?
     +         ,KWFOUND            !- if a keyword is known
C-------------- viewing parameters, etc. -------------------------------
      COMMON /PALETTE/PAL
      INTEGER   PAL (3,0:255)      !- 'current' colour palette
     +     ,ID_FACE (6)            !- 0/1/2 to force draw a facet type

      REAL     EYE (3)             !- coord of 'eye'
     +        ,COA (3)             !- coord of 'centre-of-attension'
     +         ,UP (3)             !- unit vector pointing 'up'
     +      ,LIGHT (3)             !- unit vector of light direction
     +         ,PM (4,4)           !- tranformation matrix
     +         ,GPT (3)            !- global coord of a point
     +         ,SPT (3)            !- screen coord of a point
     +        ,SPT2 (3)            !- another s.c. of a point
     +         ,LD,LDD             !- the load-step and 'anim' increment

      CHARACTER  LINE *132         !- general line of text
     +          ,FILE *80          !- 'new' inported file name
     +      ,FILE_DAT *80          !- DANPLOT date file name
     +       ,PICFILE *30          !- for PCX/print output
c    +         ,CHAR1 *1           !- a single character ?
     +          ,CKEY *1           !- ** character for event-parsing **
     +          ,NTXT *4           !- 'small' string for labels
     +     ,FILE_ROOT *30          !- for animation PCX files (cf PCX)
     +           ,ESC *1           !- ascii (27)
     +           ,C_T *6           !- == 'top corner of the screen' !
     +          ,PATH *60          !- current PATH (for SELECT_FILE)
     +       ,CURDIR@ *50          !- current directory (function)
     +        ,CMNAM@ *60          !- the command-line
     +        ,KWFILE *80          !- *NEW* keyword file
     +       ,KEYWORD *80          !- *NEW* current keyword

C-----------------------------------------------------------------------
C------ explicit 'short-integers' variable types for ftn77 calls
C... I should really aim to remove all of these into 'interface'
C... subroutine to avoid 'cluttering' up the main prog

      REAL*4  TIME (20), TIM       !- for status info

      INTEGER*4 
     +          BUFFER             !- memory loc. of an image
     +         ,IMAGE (MIMAGES)    !- memory locs. of stored frames
     +         ,VSCREEN            !- memory loc. 

      INTEGER*2
     +          IXM_S,IYM_S        !- mouse position
     +         ,IBUT_S             !- mouse button status
     +         ,IFAIL              !- if a subroutine 'fails'
     +         ,KEY                !- GET_KEY@ keycode

      INTEGER   
c    +          IRESCP             !- # bit-planes of colour
     +          IBITPL             !- # of col-planes in image (PCX) (obs)
     +         ,RES (6)            !- store of current screen sizes (obs?)


C-----------------------------------------------------------------------
C---------------- some initial data values -----------------------------
      DATA U0/40/ ,U3/13/ ,EXIT/.FALSE./
      DATA N_X_L /0/          !- The number of extra-lines on the mesh

C------------------- viewing parameters --------------------------------
      DATA
     +     NEL,NN,NDIM /0,0,0/     ! initially no elements or nodes 
     +   , NEL_OLD /-1/             ! no known elements hence FSTRIP
     +   , NLDS    /0/              ! no loadsteps present      (in CV?)
     +   , LD      /0./             ! current load-step number  (in CV?)
     +   , FACT    /0./             ! disp scale factor      +
      DATA
c    +     IOUT_D/1/                ! printer type
     +     AR      /-1./            ! aspect ratio for this device
     +   , UP      /0.,-1.,0./      ! 'up' vector = y-axis   
     +   , FEYE    /20./            ! viewing angle          \
     +   , DEYE    /-3./            ! viewing distance       + Integers ?
     +   , XML,YML /-20., 20./      ! angle of light vector  +
     +   , XME,YME /-30., 15./      ! angle of the eye       / 
     +   , SSC     /5./             ! image scale factor
c    +   , XDS,YDS,ZDS /1.,1.,1./   ! anisotropic disp. scale factors
     +   , C_R /0.,0.,0.,1.,1.,-1./ ! contour scale factors
     +   , IOP_ANIM /10*0/          ! animation token list

     +   , ITEM2C  /1/              ! - token for a CV() entry :-)
      DATA
c    +     CCOL/14,1,14,1,-1,-1/    ! colours:use,from,to,offset/infill,wrap
     +     ID_FACE  /6*0/           ! facet# to draw (0='selective')
     +   , IMAGE    /MIMAGES*-1/    ! saved images (switch off)
     +   , IMAGEN   /0/             ! # of auto-save images
     +   , IPIC_N   /1/             ! (animation picture number)
     +   , LD_CE    /0/             ! auto-glass:1=built-up,2=cut-out !
      DATA
     +     CMN_EXHAUSTED /.FALSE./  ! command-line arguments ?
     +   , KWF_EXHAUSTED /.TRUE./   ! no more keywords ?
     +   , L_OVERLAP     /.false./  ! if Facets overlap the menu
C---------------- colour info vector -----------------------------------
c >>> see SET_CV 'object'

C--------------------- output options ----------------------------------
C.... should ALL be as ONE output dest. option ! 'cept bufferin'
C.... IROT is really just part of the 'Transformation matrix'
      DATA IOP_BUF/0/, IROT/0/

C---------------------- initial constants ------------------------------
c     OPEN (19,file='debug')     !- switch on ONLY when necessary !
      IOUT_D = 1              !- to the 640x480 screen ! (obs)
      IDEST = 0               !- to VGA screen 
      DTR = 3.14159265/180.   !- degrees to radians
      ESC = CHAR(27)          !- ascii 'escape'
      C_T =  esc//'[2;2H'     !- position top-left of screen 
      DO i=1,20
        IC_MATS(I) = I+2      !- set material-cols to be 'themselves'
      ENDDO                   !- (never changed !)
      PATH=CURDIR@()          !- default directory is the current one ?
c     U1 = U1                 !- (just a dummy)

      CALL SET_CV (999,999,0)                 !- reset CV

C--------------------- check command-line for HELP ---------------------
      LINE = CMNAM@()        !- a token

      IF (LINE.EQ.'-?') THEN        !- direct request for help
        CALL WRITE_HELP()
        STOP
      ELSEIF (LINE.EQ.' ') THEN     !- no command-line given
        CALL WRITE_HELP()
        CALL COUA  ( 'wait ' )
        DO I=1,20
          CALL SLEEP@ (.2)
          CALL COUA ('.')
        ENDDO
      ENDIF
      CALL   CMNAMR()        !- rewind command-line

C-------------------------- into graphics ------------------------------
c... and 'other' initialisation ..

      CALL SET_PAL (11,0,0,0,13)    !- old 'default' set ? 
      CALL SET_PAL (10,0,0,0,13)    !- default menus
      CALL INTO_GRAPHICS ()         !- into VGA
      CALL SET_ALL_DACS()

      CALL STANDARD_MENU_FONT (-1)         !- reset the font
      CALL POST_MENUS (1, 0,0,Ckey,IVAL)     

      CALL RESET_SIZES  (RES)         !- ie x,y sizes of the window
      CALL DANPLOT_LOGO (RES)

c... yuk.. use a subroutine to 'set' resolution of various modes
C.. then call 'set_current_mode' by calling this !

c      IRESC = 16
c      IRESCP = NINT (LOG(REAL(IRESC))/LOG(2.))

c... do this as a subroutine  .. 'select-vga-mode'
C----------------------- reset mouse -----------------------------------
c  >> into_graphics does this :-)


C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C        S T A R T    O F    T H E    M A I N    L O O P
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 1110 CONTINUE       !- return-point *after* a redraw/ event
c       IF (L_OVERLAP) CALL POST_MENUS (1, 0,0,Ckey,IVAL) !- if obscured

        CALL CLOCK@(TIM)        !( maybe only call if we have to ?)
        TIME0 = TIM
        TIME1 = 0    !- Initialise
 1111 CONTINUE       !- return-point if no event was found
        CALL CLOCK@(TIM)
        TIME2 = TIME1             !-------- *Pause* clock -------------
        TIME1 = TIM - TIME0         
        IF (INT(TIME1).NE.INT(TIME2)) THEN       !- only every second
          WRITE (LINE,'(F5.0)') TIME1
          CALL CLEAR_SCREEN_AREA@ (600,22,637,38,2 )  ! -grey it
          CALL SET_TEXT_ATTRIBUTE@ (106, 1.2, 0.,0. )
          CALL DRAW_TEXT@ (LINE(1:4), 598, 35, 3)
c          IF (TIME1.GT.1.AND.TIME1.LT.2) THEN     !- redraw as we have the time
c            CALL POST_MENUS (1, 0,0,Ckey,IVAL)
c            CALL SLEEP@ (1.)   !- a little rest
c          ENDIF     
        ENDIF

c..... first we should tet to see if we have an 'event' in the 
c..... currently opened file
c..... if NOT then see if there is another file in the command line
c..... if NOT 
c     CALL SLEEP@(.1)              !- just a little rest :-)


c-----------------------------------------------------------------------
c------------------------ try the keyboard -----------------------------
c-----------------------------------------------------------------------
      CALL GET_KEY1@ (KEY)
      CALL GET_KEY_CHAR (KEY,CKEY)     !- to its single 'char'
      KEY2 = -99
      IF (CKEY.NE.'z') GOTO 1112       !- handle this 'key'

c-----------------------------------------------------------------------
c------------------------- try the menus -------------------------------
c-----------------------------------------------------------------------
      CALL GET_MOUSE_BUTTON_PRESS_COUNT@ (0, IBUT_S) 
      IF (IBUT_S.NE.0) THEN
  94    CALL GET_MOUSE_POS@ (IXM_S,IYM_S) 
        IXM = IXM_S
        IYM = IYM_S
        CKEY = 'z'                               !- defualt = 'no-op'
        CALL POST_MENUS (2,IXM,IYM,CKEY,IVAL)    !- get the opcode ?
        IF (IVAL.EQ.-2) IVAL = -99               !- convert (fudge!)

      CALL GET_MOUSE_BUTTON_PRESS_COUNT@ (0, IBUT_S) 
c       call get_mouse_position@ (ixm_s,iym_s,ibut_s)     !- obs !
        IF (IBUT_S.NE.0) GOTO 94      !- until mouse released

        IF (ckey.eq.' ') goto 1111    !- a 'non-event' so loop-back
        IF (ckey.EQ.'z') goto 1111    !- no events to handle so loop-back
        KEY2 = IVAL                   !- (translate)
        GOTO 1112                     !- to event-handler
      ENDIF

c-----------------------------------------------------------------------
c----------------- try the next file on the command-line ---------------
c-----------------------------------------------------------------------

      IF (.NOT. CMN_EXHAUSTED.and.KWF_EXHAUSTED) THEN
        LINE = CMNAM@()        !- a token
        IF (LINE.EQ.' ') THEN      !----- all the command-line is done ----
          CMN_EXHAUSTED = .TRUE.
        ELSEIF (LINE(1:1).NE.'-') THEN   !----- a Keyword file --------
          KWFILE = LINE
          OPEN (U0,FILE=KWFILE,STATUS='OLD') 
          KWF_EXHAUSTED = .FALSE.
          U1 = U0                         !- Start from the base file = U0
          IKEYWORD = 0                    !- & reset the 'counter' too 
        ELSEIF (LINE.eq.'-?') THEN 
          CALL TEXT_MODE@()            !- hmm already handled above ?
          CALL WRITE_HELP ()
          STOP
        ELSEIF (LINE.eq.'-d') THEN   !---- import a fixed-format file ----
          FILE = CMNAM@()                        
          CKEY = '#'
          KEY2 = 21                  
        ELSEIF (LINE.eq.'-%') THEN   !------- direct x<->y flip ---------
          IROT = 1-IROT
        ELSEIF (LINE(1:2).Eq.'-p') THEN    !------ palette number -------
          READ (LINE,'(2X,I5)') ICOL
          CALL SET_PAL ( ICOL ,999,999,999,CV(18))   !- WHY cv(18)
        ELSEIF (LINE.eq.'-D') THEN   !------- just draw the image -------
          GOTO 901
        ELSEIF (LINE(1:2).eq.'--') THEN      !- all other options ! :-)
          READ (LINE,'(2X,A1,I5)') CKEY      ! (KEY2 ?)
        ELSE
          PRINT*,'** Unknown command-line arguement **'
          PRINT*
          CALL WRITE_HELP ()
          STOP
        ENDIF  !- end_of_possible command_line_arguements
      ENDIF   !- block of code for reading the command-line

c-----------------------------------------------------------------------
C----------------------- Parse a data file -----------------------------
c-----------------------------------------------------------------------
c.... get the next keyword: may skip out  at EOF or return keyword

c     DO WHILE (.NOT.KWF_EXHAUSTED)       !(loop forever) 
          IF (.NOT.KWF_EXHAUSTED) THEN      ! do one-at-a-time
        IKEYWORD = IKEYWORD + 1
        CALL GET_KEYWORD (U1,U0,KEYWORD)
        WRITE(*,'(I3,A,A)') IKEYWORD, ' : ',KEYWORD

c------------------------ find the '*keyword' --------------------------
      KWFOUND = .TRUE.
      IF (KEYWORD.EQ.'*EOF') THEN
          KWF_EXHAUSTED = .TRUE.
          CLOSE (U0)
          GOTO 901    !- and draw the mesh
C------------------------------------------
      ELSEIF (KEYWORD.EQ.'*DRAW') THEN           !- so cause a redraw
        GOTO 901
C------------------------------------------
      ELSEIF (KEYWORD.EQ.'*CONTROL') THEN
        CALL R_OPTIONS (U1,CV)                  !- read options
C------------------------------------------
      ELSEIF (KEYWORD.EQ.'*SLEEP') THEN
        READ (U1,*) TIM                         !- just a 'pause'
        CALL SLEEP@ (TIM)
C------------------------------------------
c      ELSEIF (KEYWORD.EQ.'*TWO_DIMENSIONAL') THEN
c        NDIM = 2
C------------------------------------------
c      ELSEIF (KEYWORD.EQ.'*THREE_DIMENSIONAL') THEN
c        NDIM = 3
C------------------------------------------
C---------------------- user-supplied keywords -------------------------
      ELSE
        KWFOUND = .FALSE.
      ENDIF
C--------------------- 'standard' keywords -----------------------------
c.. note. use of PL2 for workspace
      IF (.NOT.KWFOUND) CALL KEY_MESH_READ
     +             (KWFOUND,KEYWORD,U1,GC,IGC,NDIM,NN,NUMS,INUMS,NEL)
      IF (.NOT.KWFOUND) CALL KEY_MESH_MUNG
     +           (KWFOUND,KEYWORD,U1,GC,IGC,NDIM,NN,NUMS,INUMS,NEL,PL2)
      IF (.NOT.KWFOUND) THEN
        PRINT*,'KEYWORD was not found !'
        STOP
      ENDIF

      ENDIF    !- block-if wrapper

c-----------------------------------------------------------------------
c.. OK I may be here if I generate CKEY from the command-line ?

      IF (CKEY.EQ.' ') GOTO 1111    !--- try again   .. why ?

C-----------------------------------------------------------------------
C------------------------ key event handling ---------------------------
C-----------------------------------------------------------------------

 1112 CONTINUE      !- got an event from MOUSE or KEYBOARD

c     WRITE (*,'(A,A,i3,A,$)') '<',ckey,key2,'>'   !- debug
 
c-----------------------------------------------------------------------
c     IF (KEY.eq.0) THEN            !- no-op
      IF (1.eq.2) THEN

      ELSEIF (CKEY.eq.'z') THEN    !- no-op

      ELSEIF (CKEY.eq.'_') THEN    !- draw image

        GOTO 901         !----- Jump straight to    * Image_Draw *

C----------------------- PRINT or SVGA re-draw -------------------------
C.. 0= 'Normal VGA screen'
C.. 1= 'SVGA redraw'
C.. 2= 'Print'

      ELSEIF (CKEY.EQ.'P')THEN
        IDEST = KEY2
        GOTO 901         !- so do a redraw

C--------- choose input file name -------------------------------------
c.... this block really needs hacking out into a subroutine !
c.... (cf. with the 'R' mirror options, which also include shifting,
c.... scaling, connectivity, *hardening* , etc.
      ELSEIF (CKEY.EQ.'#')THEN
C... key2 = 0 new data, =1 re-read = 2, append data etc. ??

c    0 = reset ?
c    1 = read a  keyword file
c   20 = Import a .PL/.PL2/.G/.NFF format file
c   21 = as 20 but already have the filename in FILE_DAT
c   30 = Export a .PL/.G/.NFF/.RAY format file

c   99 = Next in a sequence ??
c========================== kill all data ==============================
      IF (KEY2.EQ.0) THEN
        NN  = 0    
        NEL = 0
c========================== keyword file ===============================
      ELSEIF (KEY2.EQ.1) THEN      
        CALL GET_SCREEN_BLOCK@ (0,0, 639,479,BUFFER)   !- save screen
        CALL TEXT_MODE@ ()
        CALL SELECT_FILE@( '.\*.*',KWFILE, IFAIL)
        CALL INTO_GRAPHICS ()      !(set DACS in here ?) (& PAL ?)
        CALL SET_ALL_DACS  ()      
        CALL RESTORE_SCREEN_BLOCK@ (0,0, BUFFER,0,ifail)

        OPEN (U0,FILE=KWFILE,STATUS='OLD') 
        KWF_EXHAUSTED = .FALSE.
        U1 = U0                         !- Start from the base file = U0
        IKEYWORD = 0                    !- & reset the 'counter' too 

C====================== import a mesh file =============================
      ELSEIF (KEY2.EQ.20.or.KEY2.eq.21) THEN      
        IF (KEY2.EQ.21) THEN      !-------- filename is given ---------
          FILE_DAT = FILE    
        ELSEIF (KEY2.EQ.20) THEN  !--------- menu-driven --------------
          CALL GET_SCREEN_BLOCK@ (0,0, 639,479,BUFFER)
          CALL TEXT_MODE@ ()
          CALL SELECT_FILE@( '.\*.*',FILE_DAT, IFAIL)
          CALL INTO_GRAPHICS ()      !(set DACS in here ?) (& PAL ?)
          CALL SET_ALL_DACS  ()      
          CALL RESTORE_SCREEN_BLOCK@ (0,0,BUFFER,0,ifail)
        ENDIF

      call clock@ (time(10))

      CALL UPCASE (FILE_DAT)             ! (only for GET_IDT ?)
      FILE = FILE_DAT                    !- remember (for a re-read)
c... also can pass to the 'Menu Title' routine
c.. IDT is *only* used to select the file-read-subroutine :-)

      CALL GET_FILE_TYPE (file_dat,idt)        ! 1=Danplot,2=OFF,4=NFF ?
      OPEN (U3,FILE=FILE_DAT,status='OLD')        !- open the file here ?

C...... next should be moved into the 'new' style surely !
      IF (INDEX(FILE_DAT,'.PL2').ne.0) THEN
        PRINT*,'reading Extra lines.....'
        READ (U3,*) N_X_L , ((XTRA_L(J,I),J=1,4),I=1,N_X_L)
      ENDIF

c-------------------------------------------------------
C.. IDT==   1=Danplot,2=OFF,3=(obs),4=NFF
      NLDS = 0         !- reset
      FACT = 0.
      LD   = NLDS

c-------------------------------------------------------
      IF (IDT.eq.1) THEN
        CALL READ_NODES (U3,IDT,NN,NDIM,GC,IGC,MNN,NEL)
        CALL READ_ELEMS (U3,IDT,NEL,NUMS,INUMS,ME,NDIM)
        NODOF = NDIM         !- for 'results' reading
        CALL READ_LOADS (U3,IDT,NLDS,GDISPS,MDF,IGDISPS,NODOF,NN)
c-------------------------------------------------------
      ELSEIF (IDT.eq.2) THEN
        CALL READ_NODES (U3,IDT,NN,NDIM,GC,IGC,MNN,NEL)
        CALL READ_ELEMS (U3,IDT,NEL,NUMS,INUMS,ME,NDIM)
        NODOF = NDIM         !- for 'results' reading
C ...... no need to read 'results'for OFF files
c-------------------------------------------------------
      ELSEIF (IDT.eq.4) THEN        ! do 'NFF' format independantly
        CALL READ_NFF (U3,NN,NDIM,
     +       GC,IGC,MNN,NEL,NUMS,INUMS,MEL,NLDS,PAL) 
        NODOF = NDIM         !- for 'results' reading
      ELSE
        STOP 'imported data file format not recognised'
      ENDIF
c-------------------------------------------------------
      CLOSE(U3)     !- so now can close the input file
c-------------------------------------------------------
C------- set a 'nice' default displacement scale factor ---------------
C ... based on the 'final' load step ? .. or all of them ?
c ... a subroutine ??
c-------------------------------------------------------
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
c.. hmm need to call RANGE here so we know what DIAG is
        CALL GET_MESH_RANGE (GC,IGC ,NN,NDIM, COA,DIAG,diag2)
        FACT = DIAG / DISM / 25.      ! should be 'device-independant'!
        LD = NLDS                     !/ show the LAST load-step /!
      ENDIF

C------------------ pre-build vector normals for OFF files -------------
c.. maybe can do for ALL elements (cf NFF,OFF,etc.)
c.. really if we have ITYPE=9 then we can build the vector-normals
c.. just for those sort of elements ??

      IF (IDT.eq.2) THEN
        NLDS = 1               !- we now have some "displacements" :-)
        CALL GET_NORMALS (GC,IGC,NDIM,NN,NUMS,INUMS,NEL,GDISPS,MDF,PL2)
      ENDIF  
c----------------------- fixed view formats ----------------------------
c... I could just 'flag' to do a fixed direction draw ??
c... or auto-base on the object shape
c... if NDIM=2 then we should auto-inhibit B_P_C

        IF (NDIM.EQ.2) THEN                 ! if 2d view straight 'on'
          XME    = 0.   !- x-eye position
          YME    = 0.   !- y-eye position
        ELSEIF (NDIM.EQ.3) THEN             ! reset for 3D ?
          XME    =-30.  !- x-eye position
          YME    = 15.  !- y-eye position
        ENDIF
c-----------------------------------------------------------------------
c       READ*                        !- pause
        CALL POST_MENUS (1, 0,0,Ckey,IVAL) !- repost (cos of 'echo')?
        GOTO 901         !----- Jump straight to    * Image_Draw * ?
      ELSE
        Print*,' unknown KEY2 option in <#> '
      ENDIF   !--- of KEY2=1 = 'read a data file'

C------------------------- clear screen --------------------------------
      ELSE IF (CKEY.EQ.'!') THEN     !- whole VGA screen (obs?)
        CALL CLEAR_SCREEN_AREA@ (0,0,639,479,0)

c  .. killed 'm' = 'post_menu'
c  .. killed '/' = 'SVGA preview
c----------------------- exit/DOS shell --------------------------------
c... I Don't need 2 of these !
      ELSEIF (CKEY.EQ.'q') THEN
c       CALL LOAD_STANDARD_COLOURS@()      !- why ?
        GOTO 999    !- exit
c      ELSEIF (CKEY.EQ.CHAR(27)) THEN      !- obsolete exit method
c        CALL LOAD_STANDARD_COLOURS@()
c       CALL TEXT_MODE@
c        STOP '--escape--'
      ELSE IF (CKEY.EQ.'d') THEN             !- DOS shell 
        PAUSE '--PAUSE--'

C-------------------------- s t a t u s --------------------------------
c.. maybe I should do the status on the bottom-line ?
      ELSEIF (CKEY.EQ.'Q') THEN
        CALL RESET_SIZES (RES)
        CALL MAKE_BOX (RES(4),RES(5),RES(1),RES(2),X,Y)     ! blank
        CALL DR_PRIM (LINE,X,Y,4,0,2)                       ! out ?
        CALL SHOW_STATUS (TIME,NN,NEL,NLDS,NFCETS_D,NFCETS,FILE_DAT)
C------------------------> animation -----------------------------------
c.. A bit cummbersome .. better to ONLY drive via macros !
c.. Or set the options in 'pieces'
c.. I can use READ_EDITED_LINE to provide preset-animations

      ELSEIF (CKEY.EQ.'A') THEN

        IF (KEY2.ge.0) THEN      ! mouse driven preset options ?
c... ie. just fill-in the appropriate positions in the ANIM table

        ELSE 
          PRINT*,'C_T'//'Enter the <animation parameters>'
          PRINT*,'eg .. 1,0,  10,5,   0,-5,  1,  1,32'
          PRINT*,'="on","no disp", eye x+10,y+5,light x+0,y+5'//
     +           '--> to 320PCX as 1:32 pics'
          READ*, (IOP_ANIM(i),i=1,9)
          IF (IOP_ANIM(7).lt.10) THEN
            print*,'enter the animation root name (eg. FREDA)'
            read(*,'(a)') FILE_ROOT
          endif
          ipic_n = iop_anim(8)        !- the first frame number ??
          IF (iop_anim(2).ne.0) THEN                
            LD  = 0.                  !- start at zero ? ... YES
c           LD  = 1.         !- try starting at the first load step ?
            LDD = REAL(NLDS) / REAL(iop_anim(9)-1)
          ELSE
            LDD = 0.
          ENDIF
        ENDIF
C------------- 'facet# to always/never draw' ---------------------------
c... again I can fill in a table off the menu  0,1,2 for each of the 6
      ELSEIF (CKEY.EQ.'E')THEN           !-- was '|' , & then was á  !!
        write(*,'(a,a,6i2)')C_T,'facet codes(1..6)?=',(ID_FACE(i),i=1,6)
        write(*,*) '0 = selective, 1=never, 2=always'
        read*,(ID_FACE(i),i=1,6)

C------------------ backface polygon cull ------------------------------
c  (Obs) .... since the menu does this directly
      ELSEIF (CKEY.EQ.'\')THEN
        IF (KEY2.ge.0) THEN
          IF (KEY2.LE.3) CV(32) = KEY2            ! b_p_c
          IF (KEY2.LE.6) CV(17) = KEY2-3          ! depth sort (obs here)
        ELSE 
          CV(32) = MOD(CV(32)+1,3) ! just toggle if keyb.
        ENDIF 

C------------- facet sub-division --------------------------------------
c.. (obs) as done directly from the mneu into CV(19)
      ELSEIF (CKEY.EQ.'@')THEN
        PRINT*,C_T,'How many sub-facets per facet (',CV(19),' ) ?'
        READ*,CV(19)
        CKEY = '_'

C------------- input macro ---------------------------------------------
c... obsolete ! ... but interesting as an example
c      ELSEIF (CKEY.eq.'k') THEN    ! (rarely used !)
c        PRINT*,C_T//'Enter macro string (_=space,~=CR)'
c        READ(*,'(a)') LINE
c        DO I=1,leng (LINE)
c          IF (LINE(I:I).eq.'~') THEN
c            CALL FEED_KEYBOARD@ (13, ifail)
c          ELSEIF (LINE(I:I).eq.'_') THEN
c            CALL FEED_KEYBOARD@ (32, ifail)
c          ELSE
c            CALL FEED_KEYBOARD@ (ints(ICHAR(LINE(I:I))), ifail)
c          ENDIF
c        ENDDO

C***********************************************************************
C------------------> set an item's color -------------------------------
C----------- (or take other appropriate action)
C***********************************************************************
C... first give the code of the entity that is going to be changed'

      ELSEIF (CKEY.EQ.'I') THEN
        IF (KEY2.eq.-99) THEN
          PRINT*,C_T,'Item Number ?'       !- via dialog
          READ*,KEY2
        ENDIF           
        CALL SET_CV (KEY2,999, 2)        !- save this number (cf ITEM2C)

      ELSEIF (CKEY.EQ.'C') THEN
        IF (KEY2.eq.-99) THEN     
          PRINT*,C_T,'colour value ?'      !- dialog
          READ*,KEY2
        ELSEIF (KEY2.eq.-10) THEN     
          PRINT*,C_T,'value ?'             !- dialog
          READ*,KEY2
        ENDIF
        CALL SET_CV (ITEM2C,999, 6)        !- get ITEM2C ??

        IF (ITEM2C.GT.0.AND.ITEM2C.LE.M_CV) THEN

          CALL SET_CV (999, KEY2, 3)     !- change a CV() entry

c.. special cases
C-------------------- glassing materials -------------------------------
c.. this reverses the MAT # of an element (should use GET/PUT_ELEMENT!)
        ELSEIF (ITEM2C.EQ.207) THEN
          DO I=1,NEL
            CALL GET_EL_IMAT (NUMS,INUMS,IEL, IMAT)
            IF (ABS(IMAT).eq.KEY2) IMAT = -IMAT  !-- flip
            CALL PUT_EL_IMAT (NUMS,INUMS,IEL, IMAT)
          ENDDO
C------------- shade palette based on an existing colour ---------------
c... does this work ??
        ELSEIF (ITEM2C.EQ.201) THEN
          IR = PAL(1,KEY2)
          IG = PAL(2,KEY2)
          IB = PAL(3,KEY2)
          CALL SET_PAL ( 201,ir,ig,ib,CV(18))
C----------------- RGB definition of a single colour -------------------
c.. keyboard only !
        ELSEIF (ITEM2C.EQ.202) THEN 
          print*,C_T,'enter the new rgb values'
     +              ,PAL(1,KEY2),PAL(2,KEY2),PAL(3,KEY2)
          read*,ir,ig,ib
          CALL SET_PAL ( KEY2 +1000 ,ir,ig,ib,CV(18))

C------------------------- overscan colour -----------------------------
      ELSEIF (ITEM2C.EQ.203) then
          CALL SET_PALETTE@ (17,KEY2)
C------------------------ load step # patch ----------------------------
        ELSEIF (ITEM2C.EQ.204) then
          LD = KEY2
C------------------------- save image to memory ------------------------
c.. hmm do we not save ALL the VGA area ?
c.. careful of '-1L'
        ELSEIF (ITEM2C.EQ.301) then
          CALL GET_SCREEN_BLOCK@ (0,0, 639,479,IMAGE(KEY2))
          IF (IMAGE(KEY2).EQ.-1) THEN
            CALL SOUND@( 256,2)
            PRINT*,'****ERROR: IMAGE',KEY2,' out of memory' 
          ENDIF

C------------------------- load iamge from memory ------------------------
c.. careful of '-1L'
        ELSEIF (ITEM2C.EQ.302) then
          IF (IMAGE(KEY2).ne.-1) THEN   !- if saved
            CALL RESET_SIZES (RES)           !- in the 'normal' window
            CALL RESTORE_SCREEN_BLOCK@   
     +     (INTS(RES(4)),INTS(RES(5)),IMAGE(KEY2),0,IFAIL)
            CALL DOSERR@(IFAIL)
          ELSE
            CALL SOUND@( 256,2)
            PRINT*,'****ERROR: IMAGE',KEY2,' was not saved' 
          ENDIF
          CALL DOSERR@ (IFAIL)
C---------------------- auto save image to memory ----------------------
c.. hmm should really just be saving the 'image' part of the screen !

        ELSEIF (ITEM2C.EQ.-9) then  !- ie. in the next available space
          IMAGEN = IMAGEN+1
          CALL GET_SCREEN_BLOCK@ (0,0, 639,479,IMAGE(IMAGEN))
          IF (IMAGE(KEY2).EQ.-1) THEN
            CALL SOUND@( 256,2)
            PRINT*,'****ERROR: IMAGE',KEY2,' out of memory' 
          ENDIF
        ELSE
          PRINT*,' ITEM2c=',item2c,' not understood'
        ENDIF



C------------------------ shade palette --------------------------------
C --> linear palatte shade range specification
c.. obsolete 'cos it is also done by SET_PAL with an opcode
      ELSEIF (CKEY.EQ.';') THEN      !- cf ('.' from the menu)
        CALL SET_SHADES (KEY2)

C------- individual color specification & pre-def palettes ? -----------
c.. eg. 'rainbow','spectrum','gray' 
      ELSEIF (CKEY.EQ.'.') THEN
        ICOL = KEY2  + 1000
        CALL SET_PAL (ICOL+1000,IR,IG,IB,CV(18))  

C------------------- palette cycling -----------------------------------
c... hmm  just spin the table ! .. (careful of the menu-colors !)
c.. hmm was 'P'
      ELSEIF (CKEY.EQ.'W') THEN
        DO I=1,5*3
          CALL SET_PAL (10,IR,IG,IB,CV(18))
        ENDDO

C----------------> 'hand' contour type  selection ----------------------
      ELSEIF (CKEY.EQ.'K') THEN
        IF (KEY2.NE.-99) THEN      ! menu pointer to contour type
          CV(13) = KEY2
       ELSE
          PRINT*,C_T//'enter facet color option (',CV(35),')'
          READ*,CV(13) 
      ENDIF

C------------------- contour range setting -----------------------------
c... an easy block to push out to a subroutine !
c... an *obvious* candiadate for form-filling !
c...... also nicer to split into ../ # conts / Min,Max / Log.power /
      ELSEIF (CKEY.EQ.'^') THEN
        IF (KEY2.eq.1) THEN    !---------- auto set range ----------
          C_R(3)= c_r(1)
          C_R(4)= c_r(2)
          C_R(5)= 1.
          C_R(6)=-1.
        ELSEIF (KEY2.eq.2) THEN
           PRINT*,'Contour range = ?  (current ='
     +       ,C_R(3),' to', C_R(4) ,' from limits of'
     +       ,C_R(1),' to', C_R(2) ,' )'
           READ*, c_r(3),C_R(4)
        ELSEIF (KEY2.eq.3) THEN
           PRINT*,'Number of contours = ? (current =', CV(18),' )'
           READ*, cv(18)
        ELSEIF (KEY2.eq.4) THEN
           PRINT*,'Log Factor = ?  (current =',C_R(5),' )' 
           READ*, c_r(5)
        ELSEIF (KEY2.eq.4) THEN
           PRINT*,'Take Absolutes (if=-1)?  (current =',C_R(6),' )' 
           READ*, c_r(6)

c        ELSEIF (KEY2.eq.2) THEN    !-- manualy set the contour range ----
c
c        write(*,'(a,e12.4,a,e12.4)')
c     +      'actual  range =',C_R(1),'<---->',C_R(2)
c     +     ,'current range =',C_R(3),'<---->',C_R(4),
c     +               ' power=',C_R(5)
c           PRINT*,'enter new Min, Max, & Power, (-ve for abs)'
c           WRITE (LINE,'(2E12.3,F5.2)') C_R(3),C_R(4),C_R(5)
c           CALL READ_EDITED_LINE@ (LINE,0s,1S,15S,IFAIL)   !- 'esc' to quit
c           read(LINE,*) C_R(3), C_R(4), C_R(5)
c           C_R(6) = SIGN (1.,C_R(5))       !- just +1 or -1
c           C_R(5) =  ABS    (C_R(5))       !- lose the sign :-)
        ENDIF

C---------- contour method (=Gouraud/lines or averaged) ----------------
c.. a bit obsolete ! .. just use the 'number' menu :-)
      ELSEIF (CKEY.EQ.'g') THEN
        IF (KEY2.ge.0.and.KEY2.le.9) THEN  !-- # conts from menu ?
          CV(18) = KEY2
        ELSEIF (KEY2.eq.10) THEN           !-- # colours
          PRINT*,'# Contours?=',CV(18) 
          READ*, CV(18) 
c.. the next was never implimented !
c        ELSEIF (KEY2.eq.11) THEN           !-- color range + offset
c          PRINT*,'palette from/to &offset?=',CCOL(2),CCOL(3),CCOL(4)
c          READ*,CCOL(2),ccol(3),ccol(4)

         ENDIF

C-------------------- screen image store/restore -----------------------
c... replay the sequence stored by animation
C... need a 'query' option to tell how many frames there are :-)
c.. also these too may easily be made into subroutines
      ELSEIF (CKEY.eq.'}') THEN    !------ forwards
        CALL RESET_SIZES (RES)           !- in the 'normal' window
        DO J=1,999
          DO I=IOP_ANIM(8),IOP_ANIM(9)
           CALL RESTORE_SCREEN_BLOCK@   
     +      (RES(4),RES(5),IMAGE(I),0,IFAIL)
           IF (KEY_WAITING@()) GOTO 1111
          ENDDO
        ENDDO
      ELSEIF (CKEY.eq.'{') THEN   !---- forwards and backwards
        CALL RESET_SIZES (RES)           !- in the 'normal' window
        DO j=1,999
          DO I=iop_anim(8),iop_anim(9)-1
         CALL RESTORE_SCREEN_BLOCK@   
     +     (RES(4),RES(5),IMAGE(I),0,IFAIL)
            IF (KEY_WAITING@()) GOTO 1111
          ENDDO                                        
          DO I=iop_anim(9),iop_anim(8)-1,-1
         CALL RESTORE_SCREEN_BLOCK@   
     +     (RES(4),RES(5),IMAGE(I),0,IFAIL)
            IF (KEY_WAITING@()) GOTO 1111
          ENDDO
        ENDDO

C------------------> glassing materials --------------------------------
c.. really need a more general method !
c.. this is easily moved into CV() :-)
        ELSEIF (CKEY.EQ.'[') then    !------ built-up models------
          LD_CE = 1
        ELSEIF (CKEY.EQ.']') then   !------ excavation models-------
          LD_CE = 2

C----------------- output destination choice ---------------------------
C ........ should really be a typed numeric option
C ---> PCX screen dump .. 'just the current drawing area ?'
C.... most of these are really printers

c... IOUT_D should really be in CV as the code for the printer 'driver'

c... Ok so 'select a printer by number' then allow the 'entries' for 
c... this to be changed (eg 75/150/300 dpi)

c... Ok just give every output possibility a code..
c---- to include 300dpi etc. if appropriate .. 
c---- so 'block' this section into a subroutine  code--> open device

      ELSEIF (CKEY.EQ.'O')THEN
        IF (KEY2.eq.-99) then
          PRINT*,C_T,'Output destination ?'
          PRINT*,'(1=VGA,2=SVGA,3=PCX,4=laser,5=PS,6=pro-p,8=HP7550)'
          READ*,key2
        ENDIF       
      IF (KEY2.eq.1) THEN               !---screen----
        IOUT_D = KEY2
c       print*,'image size (640,480) ?'
c       read*,RES(1),RES(2)
        RES(1) = 640
        RES(2) = 480
c       IRESC  =  16
c       IRESCP =   4
        AR = -1.   !.. but NOT=-1.if 320x200 mode etc.

C      ELSEIF (KEY2.eq.2) THEN
C        IOUT_D = KEY2
C        RES(1) = 1024
C        RES(2) = 768
C        RES(4) = 0
C        RES(5) = 0
C        IRESC  = 256
C        IRESCP =   8
C        AR = -1.
C        CALL GRAPHICS_MODE_SET@ (INTS(RES(1)),INTS(RES(2))
C     +                ,INTL(IRESC),IFAIL)

      ELSEIF (KEY2.EQ.3) THEN             ! redraw to a PCX
        CV(51) = 5

      ELSEIF (KEY2.EQ.4) THEN    ! ---> laser-printer output
        IOUT_D = KEY2
        PICFILE='LASER'
        PRINT*,'printer resolution ? (75,100,150,300)'
        READ*,IPRES
        CALL SELECT_PCL_PRINTER@ (0,'A4',IPRES,RES(1),RES(2))

        CALL OPEN_GPRINT_FILE@     (PICFILE,IFAIL)
        AR = -1.

      ELSEIF (KEY2.EQ.5) THEN    ! ---> Postscript output
        CV(51) = 5             ! 5= 'ordinary' postscript (a4 portrait)

      ELSEIF (KEY2.EQ.6) THEN    ! ---> pro-printer output
        IOUT_D = KEY2
        CALL OPEN_GPRINT_DEVICE@     (1,IFAIL)
        RES(1) = 960
        RES(2) = 576
        AR = -1.3

        CALL SELECT_DOT_MATRIX@ (0,960,576)

      ELSEIF (KEY2.eq.8) THEN         ! ---> HP7550 plot output
        IOUT_D = KEY2
        PRINT*,C_T//'HP7550 File name:'
        READ(*,'(A)') PICFILE
        CALL OPEN_PLOT_FILE@ (PICFILE,IFAIL)
        RES(1) = 10870                        !-store the sizes in a 
        RES(2) =  7600                        ! table? so can modify
        AR = 1.                        !- 'border widths ??'

        ENDIF  ! end of ouput devices
        CKEY='_'

C--------------------- dump the full screen to a PCX -------------------
      ELSEIF (CKEY.eq.'"') THEN           !---(PCX screen dump) --
        CALL GET_SCREEN_BLOCK@ (0,0,639,479,BUFFER)
        PRINT*,C_T//'PCX screen dump File name: ?'    
        READ(*,'(A)')PICFILE
        CALL SCREEN_BLOCK_TO_PCX@ (PICFILE,BUFFER,IFAIL)
        CALL RETURN_STORAGE@ (BUFFER)
        CALL PCX_FIX (PICFILE,PAL)

C---------------- screen window size and position ----------------------
      ELSEIF (CKEY.EQ.'w') THEN
        PRINT*,C_T//'enter window size, x,y (460,460)'
        READ*,RES(1),RES(2)
        RES(4) = (479 -RES(1)) / 2     !- force it to the 
        RES(5) = (479 -RES(2)) / 2     !- centre_of_the_screen

C-------------------- buffered screen output ---------------------------
      ELSEIF (CKEY.EQ.'b') THEN       ! is this a CV option ?
        IOP_BUF = MOD(IOP_BUF+1,2)

C------------------------- s c a l i n g -------------------------------
C --> menu-driven  eye-angle/eye-pan/light-source chenaging etc. ?
C.. can make a subroutine .. nicer if we hold ALL these real numbers in
C.. a 'viewing-array' :
C.. EYE angles & vector
C.. LIGHT angles and vector
C.. COA on the object
C.. SSC scale factor
C.. DEYE perspective distance
C.. UP vector
C.. PAN screen x-y shifts (cf .. moving COA itself)
C.. really a BIG fill-in form (some affect others ?) EYE->unit vector

      ELSEIF (CKEY.EQ.'V') THEN
        IF (KEY2.EQ. 1) YME = YME - 5.     !- rotate the eye
        IF (KEY2.EQ. 2) XME = XME + 5.     !- by 5 degress
        IF (KEY2.EQ. 4) XME = XME - 5.
        IF (KEY2.EQ. 5) YME = YME + 5.
        IF (KEY2.EQ. 3) THEN
          PRINT*,C_T//'Enter the EYE longitude and latitude'
          PRINT*,'XME,YME=',XME,YME
          READ*,  XME,YME
        ELSEIF (KEY2.GE.6.AND.KEY2.LE.9) THEN
           CV(38) = KEY2 - 5         !- a 'preset' code for CAMERA x/y/z/p
        ENDIF

        IF (KEY2.EQ.11) CV(47) = CV(47) + 10.  !- shift the image
        IF (KEY2.EQ.12) CV(46) = CV(46) - 10.  !- by 10% screen width
        IF (KEY2.EQ.14) CV(46) = CV(46) + 10.
        IF (KEY2.EQ.15) CV(47) = CV(47) - 10.
        IF (KEY2.EQ.13) THEN
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
          PRINT*,'XML,YML=',XML,YML
          READ*,XML,YML
        ENDIF

        IF (KEY2.eq.31) THEN     ! change the perspective !
          DEYE = DEYE / SQRT(2.)   ! wide-angle
          SSC  = SSC / SQRT(2.)
        ELSEIF (KEY2.EQ.32) THEN  
          DEYE = DEYE * SQRT(2.)   ! narrow-angle
          SSC  = SSC * SQRT(2.)
        ENDIF

C------------------------ typed eye-position ---------------------------
      ELSEIF (CKEY.EQ.'t')THEN 
        PRINT*,C_T//'Enter the EYE longitude and latitude'
        PRINT*,'XME,YME=',xme,yme
        READ*,XME,YME

C------------------ cursor key eye-position move -----------------------
      ELSEIF (CKEY.eq.CHAR(128)) THEN
        XME = XME + 5.
      ELSEIF (CKEY.eq.CHAR(129)) THEN
        XME = XME - 5.
      ELSEIF (CKEY.eq.CHAR(130)) THEN
        YME = YME - 5.
      ELSEIF (CKEY.eq.CHAR(131)) THEN
        YME = YME + 5.
C--------------------------- image scaling -----------------------------
      ELSEIF (CKEY.EQ.'+') THEN      !-- do with 'V' above ??
         SSC = SSC * SQRT(2.)
      ELSEIF (CKEY.EQ.'-') THEN  
         SSC = SSC / SQRT(2.)
      ELSEIF (CKEY.EQ.'*') THEN  
        WRITE(*,'(A,E8.2,A)') 'image scale ?(',SSC,' )'
        READ*,SSC
C------------------------- displacemnt scaling -------------------------
      ELSEIF (CKEY.EQ.'<') THEN
        FACT = FACT / SQRT(2.)
      ELSEIF (CKEY.EQ.'>') THEN 
        FACT = FACT * SQRT(2.)
      ELSEIF  (CKEY.EQ.'s') THEN
        WRITE(*,'(A,E8.2,A)') 'displacement scale ?(',FACT,' )'
        READ*,FACT
C--------------- typed anisotropic displacemnt scaling -----------------
c      ELSEIF (CKEY.EQ.'d') THEN
c        WRITE(*,'(A,3F8.2,A)') 
c     +    '(',XDS,YDS,ZDS,' ) new x,y,z disp. scale ?'
c        READ*,XDS,YDS,ZDS
C-------------------- typed centre -of-attention -----------------------
      ELSEIF (CKEY.EQ.'') THEN  
        WRITE(*,'(A,3F10.4)') 'COA ?',COA(1), COA(2), COA(3)
        READ*,COA(1), COA(2), COA(3)
        CKEY = '_'
C---------------------- (screen x,y flipping) --------------------------
      ELSEIF (CKEY.EQ.'%') THEN 
        IROT = 1-IROT

C-------------------- 'Mirroring the model' ----------------------------
C---> push the work to a subroutine to save space ? :-)
C... also allow the option to 'remove' a block to 'un-mirror' 
c.. SUBROUTINE would need NUMS,GC,GDISPS,IDIR,VAL (xyz limits?)
c.. maybe *redo* Facet-stripping to save 'messing about'    <- Yes!

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

c-------------------------------------------------------------
        IF (IDIR.EQ.7) THEN     !  un-mirror.. = 'halve the model'
          NEL = NEL * 1./2.     !**** careful with the disps !
          NN  = NN  * 1./2.
        ELSEIF (IDIR.EQ.8) THEN !  'quarter the model' eg. (spud-can)
          NEL = NEL * 3./4.
          NN  = NN  * 3./4.
        ELSEIF (IDIR.EQ.8) THEN !  'eighth the model' eg. (spud-can)
          NEL = NEL * 7./8.
          NN  = NN  * 7./8.

c-------------------------------------------------------------
        ELSEIF (IDIR.EQ.50) THEN !  *** Connectivity ***
          CALL CONNECTIVITY_IMAT (GC,IGC,NN,NDIM
     +                               ,NUMS,INUMS,NEL,PL2)

c-------------------------------------------------------------
        ELSEIF (IDIR.EQ.51) THEN !  * HARDENING* : GC = GC + FACT*DISPS
          DO I=1,NN
            DO J=1,3
              GC(J,I) = GC(J,I) + FACT* (FAC1 * GDISPS(J,IB1+I) 
     +                                 + FAC2 * GDISPS(J,IB2+I) )
            ENDDO
          ENDDO
          NLDS = 0     !- 'wipe-out'
c-------------------------------------------------------------
        ELSEIF (IDIR.GE.1.AND.IDIR.LE.3) THEN
          IF (II.lt.0) VAL = GC(IDIR,NN+1)  ! min x,y,z is stored at NN+1
          IF (II.gt.0) VAL = GC(IDIR,NN+2)  ! max x,y,z is stored at NN+2

C--------------------- copy the displacement info ----------------------
C ----> assume the NN is constant, so copy from the top down 
 
        DO K=NLDS-1,0,-1             ! loop load-steps
          IB1 = INTL(K) * NN * 2
          IB2 = INTL(K) * NN
          DO I=NN,1,-1               ! loop nodes
            DO J=1,3                 ! loop freedoms
                GDISPS (J,IB1 +I   )  =  GDISPS(J,IB2 +I)  ! move + 
              IF (J.NE.IDIR) THEN
                GDISPS (J,IB1 +I+NN)  =  GDISPS(J,IB2 +I)  ! copy
              ELSE                         
                GDISPS (J,IB1 +I+NN ) = -GDISPS(J,IB2 +I)  ! reflect
              ENDIF
            ENDDO     
          ENDDO
        ENDDO
c.....................
        CALL MIRROR_ELEMS (GC,IGC,NN,NDIM,NUMS,INUMS,NEL,IDIR,VAL)

        ENDIF  !--- block for IDIR =1,=2, or =3

C----------------- linkage of load step # to picture #------------------
c... maybe  more possibilities
      ELSEIF (CKEY.EQ.'F') THEN 
        IF (KEY2.eq.1) THEN         !--  link load step # --
          IF (CV(37).eq.1) THEN
            CV(37) = -1         ! switch-off
          ELSE
            CV(37) = 1          ! switch-on 
          ENDIF
        ELSEIF (KEY2.eq.2) THEN     !--  link view direction --
          IF (CV(37).eq.2) THEN
            CV(37) = -1         ! switch-off
            CV(38) = -1         ! view = non-auto (ie. can rotate)
          ELSE
            CV(37) = 2          ! switch-on 
            CV(38) = 1          ! start at view #1  ...to  #4
          ENDIF
        ENDIF
C-------------------- sum/difference load steps-------------------------
c.. quite easy as a subroutine :-)
c.. just need the GDISPS table and NN,NODOF,NLDS 
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
c.. (obs)?) either +1, -1, or typed .. so do from the 'number-menu' only ?
c.. hence 'all' typed values would be prompted by their title :-)
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
c.. (obs?) again this can be done direct from the menu (also +/-1??) 
      ELSEIF (CKEY.EQ.'S') THEN
        print*,C_T//'ishrink=', cv(31),' ?'
        read*,cv(31)
C------------------- x,y,z up-direction toppling -----------------------
C..   better as a subroutine to save space ??
C      ( for Ian Williams et al. who think lop-sided :-)
C...  we can also do this by simply adjusting the mapping of the model 
c..   'xyz' onto the 3d screen model xyz ?
      ELSEIF (CKEY.EQ.'a') THEN
        T      = COA(1)
        COA(1) = COA(2)
        COA(2) = COA(3)
        COA(3) = T
        DO I=1,NN  + 2   ! ie. all the nodes AND the x,y,z min/max
                   T  = GC(1,I)  
          GC(1,I) = GC(2,I)
          GC(2,I) = GC(3,I)
          GC(3,I) = T
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
      ELSE
        WRITE (*,'(A,$)') '?'
      ENDIF
C-----------------> if CKEY='_' then we generate a refresh  

c     IF ( IBUT2.ne.0) CKEY='_'                    ! mouse-redraw
c     IF ( IBUT1.ne.0.or.KEY_WAITING@()) CKEY=' '  ! don't redraw ??
c     IF ( .not.(CKEY.eq.'_'.or.KEY.eq.13))  GOTO 1111   ! go back

      GOTO 1111               !-- back again !

C (note. KEY=13 was a HARD redraw (no interupts alllowed) )
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C----- OK either we have moved or a new load case .. so draw -----------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

  901 CONTINUE


C--------- turn the NUMS into a list of external facets 'FACETS' -------
c... only need to do if the table of elements has changed
      IF (NEL.NE.NEL_OLD) THEN
        NEL_OLD = NEL
c       print*,' Re-sorting the mesh...'
        call clock@ (time(11))
        CALL GET_MESH_RANGE (GC,IGC, NN,NDIM, COA,DIAG,diag2)  ! get DIAG,etc.
        CALL FSTRIP2 (NEL,NUMS,INUMS,FACETS,PL2,IFCETS,NFCETS,NDIME)
        call clock@ (time(12))
        call sound@ (256,1)
      ENDIF
C-----------------------------------------------------------------------

      CALL CLOCK@ (TIME(1))       !--- start of drawing --

c----------------- set the picture size parameters --------------------
      IF (IDEST.EQ.0) THEN        !- ordinary image draw
c.. set the screen 'window' sizes ?
        CALL RESET_SIZES (RES)           !- in the 'normal' window
        IRESC  = 16
        IRESCP =  4
        AR     = -1.       !--- noraml 'Aspect Ratio'

      ELSEIF (IDEST.EQ.1) THEN    !- SVGA preview mode
c.. set 1024x768 etc. the go into this mode

c       CALL RESET_SIZES (RES)           !- in the 'normal' window
        RES(1) = 1024     !-\ Hack for now
        RES(2) = 768      !-/ 
        RES(4) =   0            !-\ No window offsets ?
        RES(5) =   0            !-/
        IRESC  = 256
        IRESCP =  8
        AR     = -1. 

      ELSEIF (IDEST.EQ.2) THEN    !- output to PRINTER
        IOUT_D = CV(51)
        IF (IOUT_D.EQ.5) THEN      !- Postscript -
c.. hmm we can set the # of colour to infinity really (IRESC) ??
          IRESC  = 256
          IRESCP =  8
          AR = 1.
          PRINT*,C_T//'PostScript File name: ?'
          READ (*,'(A)') PICFILE
          IF (PICFILE.EQ.' ') PICFILE = 'PRINT'   !- so will print it
          CALL UPCASE (PICFILE)
          CALL SELECT_PS_PRINTER (PICFILE,RES(1),RES(2),
     +          RES(4),RES(5), FILE_DAT) 
        ELSEIF (IOUT_D.EQ.3) THEN      !- PCX -
          PRINT*,'X,Y,(D) RESOLUTION ?  (1200,1200,4)'  !- set a default
          READ*,RES(1),RES(2),IBITPL
          PRINT*,C_T//'PCX File name: ?'
          READ(*,'(A)')PICFILE
        AR = -1.
        ENDIF


      ENDIF    !- output destinations

C----------------- create transformation matrix ------------------------
c... is cv(38) a cheat to fix some view directions
c... if so, it is better to mung EYE() beforehand
c... I should also create the screen scale/shift matrix then product :-)

      IOP = CV(38)
      CALL XM2EYE (XME,YME,EYE,iop)      ! convert to a unit vector
      CALL XM2EYE (XML,YML,LIGHT,0)      !    "    "  "  "     "
       
      T = SIGN (1.,ABS(YME)-90.001)      ! 'fly -over-the-top' patch
      DO I=1,3
        UP(I) = SIGN( UP(I),T )
      ENDDO
      CALL CAMERA (EYE,COA,UP, DEYE*DIAG,FEYE, PM,cv(38))  ! create PM()

c------- set the aspect ratio for different ouput devices --------------
c ...?

C------------- set up the windowing parameters -------------------------
c.. ie multi-picture # --> x/y --> hilo of the window
c. rem. cv(41) = #, 42=nxpics,43=nypics, res = image width/height
c.. note that this uses RES(1..5) so.. need to call-up
c.. the info for the current device .. eg if SVGA
c.. * IF AR < 0 then the pic_origin moves AND iwin_yw inverts
c.. * also 'tied' numbers move too ? (eg. pic #)

      iwinx =  mod(cv(41)-1  , cv(42) )
      iwiny =     (cv(41)-1) / cv(42)
      iwin_xw = res(1)/ cv(42)
      iwin_yw = res(2)/ cv(43)
      iwin_xo = iwin_xw * iwinx  + res(4)  !.. the picture offsets
      iwin_yo = iwin_yw * iwiny  + res(5)  !..
      IF (AR.GT.0) THEN
        iwin_yw = -iwin_yw
        iwin_yo = res(2) - iwin_yo    !- invert if origin is bottom-left ?
      ENDIF
c.. SC_X/Y/Z are (ugly #) screen scale factors 'hacked' onto 
c.. subr. TRANSFORM to adjust image szie AND aspect ratio
c.. X/Ycen are also hacks for the centre of the 'viewable' window
c.. CV(46,47) are the % shift of the centre (-50% to 50% ?)

      SC_X = SSC * iwin_xw         ! x,y,z scaling factors (screen)
      SC_Y = SSC * iwin_xw * ar
      SC_Z = SSC * iwin_xw
      XCEN = iwin_xw/2. + iwin_xw/2.*CV(46)/100.     ! image-centre
      YCEN = iwin_yw/2. + iwin_xw/2.*CV(47)/100.     ! + panning
      IF (IOP_BUF.ne.1) then   
        XCEN = XCEN + iwin_xo  ! add picture offsets if NOT buffered
        YCEN = YCEN + iwin_yo  ! 
      ENDIF

C-------------- set up the factors of each loadstep to add -------------
c.. ILD1,2 are the lo,hi 'adjacent' load steps  
c   (maybe do externaly ?) <- but would 'break' Eigenmode plots !
c  << unless I hack eigenmodes to a 'set' of elements !! >>

      ILD2 = MAX(1,(MIN(INT(LD+1),NLDS)))
      ILD1 = ILD2 - 1
      FAC2 = LD - REAL(ILD1)
      FAC1 =  1. - FAC2
      IF (ILD1.EQ.0) THEN        !- patch
        ILD1 = 1
        FAC1 = 0.
      ENDIF
      IB1  = INTL(NN) * (ILD1-1)     !- pointers to the table
      IB2  = INTL(NN) * (ILD2-1)     
c      XDS1 = XDS * FACT * FAC1     !- anisotropic scale factors
c      YDS1 = YDS * FACT * FAC1     !- >> maybe do this AND 'juggling'
c      ZDS1 = ZDS * FACT * FAC1     !- via a 'fudge' matrix (3x3) !
c      XDS2 = XDS * FACT * FAC2     !- therefore do on the 'outside'
c      YDS2 = YDS * FACT * FAC2     !- (cf the 'contour' vector)
c      ZDS2 = ZDS * FACT * FAC2

C------------ auto-glassing of built-up/excavated models ---------------
c.. again this should maybe be done 'before' we call the 'plot' routine
c.. >> be carefull of 'auto'-animation
c... using GET/PUT element is easy :-)
      IF (LD_CE.eq.1) THEN  !--built-up models  (ie. mat. 'appears')
        DO I=1,NEL
          CALL GET_EL_IMAT (NUMS,INUMS,IEL, IMAT)
          IF (IMAT.NE.1.and.IMAT.GT.LD+1.5) IMAT = -IMAT  !<-- mend this!
          CALL PUT_EL_IMAT (NUMS,INUMS,IEL, IMAT)
        ENDDO
      ELSEIF (LD_CE.eq.2) THEN  !--excavated models  (mat. 'disappears')
        DO I=1,NEL
          CALL GET_EL_IMAT (NUMS,INUMS,IEL, IMAT)
          IF (IMAT.NE.1.and.IMAT.LT.LD+0.5+1.01) IMAT = -IMAT
          CALL PUT_EL_IMAT (NUMS,INUMS,IEL, IMAT)
        ENDDO
      ENDIF

C------------ loop the facets to get the 'draw-list' -------------------
c .. I could mark the 'touch' column as -ve if the adjoining material 
c    is of a different material type?--> so can 'always' draw?
c--> could be a subroutine : from FACETS and GC to 'list-to-draw'
c  >> note the effects of B_P_C and 'Force_draw'

c.. can I kill SCOORD and just 'recreate' ??

      NFCETS_D = 0
      L_OVERLAP = .false.      !/* Assume no facets are partly visible */
      DO IFC=1,NFCETS
        IEL   = FACETS(1,IFC)
        IFACE = ABS(FACETS(2,IFC))
        ITOUCH= FACETS(3,IFC)

        CALL GET_EL_IMAT (NUMS,INUMS,IEL, IMAT)
        IF (IMAT.LE.0) GOTO 701    ! 'skip' (as 'invis')

        IF (ID_FACE(IFACE).eq.1) GOTO 701   ! never 'draw' this side ?
        IF (ID_FACE(IFACE).eq.2) GOTO 777   ! always 'draw' this side ?

        IF (ITOUCH.eq.0) GOTO 777   !  'draw' (as it is a 'boundary')


        IEL2 = FACETS(1,ITOUCH)  ! get info on its toucher'
        CALL GET_EL_IMAT (NUMS,INUMS,IEL2, IMAT2)
        IF (IMAT2.gt.0) GOTO 701  ! 'skip' (cos 'touch' is present !)

c---------- now get & transform this facets' nodes --------------------
  777 CONTINUE
      CALL GET_ELEMENT 
     +   (NUMS,INUMS,IEL, NUM, NOD,NDIME,ITYPE, IMAT,IU1,IU2,IU3)
      CALL GFACE (NOD,NDIME,ITYPE,IFACE,FS,NN_F,NN_FT, 2 )

c     print*,'iel=',iel,' nn_f=',nn_f,' nn_ft=',nn_ft

C .. ok turn FS into the 'global node numbers .. not just the 'local' !
C      Z_MIN,Z_MAX,Z_MEAN ???
      N_VISIBLE = 0           ! # of nodes of this facet that are on-screen
      DO II=1,NN_FT
        I = NUM(FS(II))       !- the node number
        FS(II) = I            !***** ARGGH! this MUNGs FS **** (obs)
        GPT(1) = GC(1,I)
        GPT(2) = GC(2,I)
        GPT(3) = GC(3,I)
        IF (NLDS.GT.0) THEN
          DO J=1,3
            GPT(J) = GPT(J) + FACT* 
     +        ( FAC1 *GDISPS(J,IB1+I) +FAC2 *GDISPS(J,IB2+I))
          ENDDO
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
          IF (SPT(2).gt.iwin_yo.and.SPT(2).lt.iwin_yo+iwin_yw) THEN
            N_VISIBLE = N_VISIBLE + 1
          ENDIF
        ENDIF
      ELSE  !-- clip to the 'buffered' image (0->iwin_xw)  ..
        IF (SPT(1).gt.0..and.SPT(1).lt.iwin_xw) THEN
          IF (SPT(2).gt.0..and.SPT(2).lt.iwin_yw) THEN
            N_VISIBLE = N_VISIBLE + 1
          ENDIF
        ENDIF
      ENDIF

      ENDDO            !-- end-of-loop-a-facets-nodes
      IF (N_VISIBLE.eq.0) GOTO 701   !-- skip as no_nodes are on_screen
      L_OVERLAP = L_OVERLAP.OR.(N_VISIBLE.LT.NN_FT) ! record 'menu overlaps'
C------------------ back-face polygon culling --------------------------
      C3 = 0.
      DO I=1,NN_F
        IP1 = MOD (I,NN_F) + 1
        C3 = C3 + (SC_F(I,2)+SC_F(IP1,2)) * (SC_F(IP1,1)-SC_F(I,1))/2.
      ENDDO
      C3 = -C3 * AR         ! weight wrt. the screen topology
      IF (NDIM.EQ.3) THEN   ! only cull for '3D' meshes
        IF (CV(32).EQ.1.AND.C3.GT.0) GOTO 701   ! cull 'fronts'
        IF (CV(32).EQ.2.AND.C3.LT.0) GOTO 701   ! cull 'backs'
      ENDIF
C--------- OK this facet IS to be drawn.. so get its 'z'-value ---------

      NFCETS_D = NFCETS_D + 1       ! = '# of Facets to draw'
      PL2(NFCETS_D) = IFC           !-- 'index' of this facet in PL1 

       IF (NDIM.eq.3.and.CV(17).NE.0 ) THEN      !- only if we ARE depth-sorting
         ZVAL = SC_F(1,3)
         IF (CV(17).EQ.2) THEN                    ! furthest point
           DO I=1,NN_FT
             ZVAL = MIN(ZVAL,SC_F (I,3))
           ENDDO
         ELSEIF (CV(17).EQ.1) THEN                ! nearest point   ?
           DO I=1,NN_FT
             ZVAL = MAX(ZVAL,SC_F (I,3))
           ENDDO
         ELSEIF (CV(17).EQ.3) THEN                ! centroid point 
           DO I=1,NN_FT
             ZVAL = ZVAL + SC_F (I,3)
           ENDDO
           ZVAL = ZVAL / REAL(NN_FT)
         ENDIF
         ZBUF(NFCETS_D) = ZVAL
       ENDIF                       !-- only if we are sorting facets

  701 CONTINUE !.......... end-of 'skip-to' (='ignore this facet')
      ENDDO  !............ end of the facet loop (IFC)

      CALL CLOCK@ (TIME(2))  !--  time to end of transformation

C---------------- depth sort the 'drawable' facets ---------------------
      IF (NDIM.eq.3.and.CV(17).NE.0) THEN
        CALL RSORT@ (PL1,ZBUF,INTL(NFCETS_D)) 
      ELSE                      !- no depth-sort 
        DO I=1,NFCETS_D         !- ..so just number sequentially
          PL1 (I) = I
        ENDDO
      ENDIF

      CALL CLOCK@ (TIME(3))  !-- time to (after) depth-sorting the facets 

C-----------------------------------------------------------------------

c... Ok so now we have an ordered list of facets to draw
c... so open any graphics devices.. draw any 'furniture'
C... then lopp facets and draw .. and finally draw any other info (titles...)

C=======================================================================
C========================  draw the new image ==========================
C=======================================================================
      CALL SET_TEXT_ATTRIBUTE@ (100,1.,1.,1.)    !-reset font
      IOUT_D = CV(51)        !- save in COMMON for device-drivers

      IF (IDEST.EQ.0) THEN        !- ordinary image draw

      ELSEIF (IDEST.EQ.1) THEN    !- SVGA preview mode
        CALL GRAPHICS_MODE_SET@ (INTS(RES(1)),INTS(RES(2))
     +                ,INTL(IRESC),IFAIL)
        CALL SET_ALL_DACS()   !- need to reset my DAC colours

      ELSEIF (IDEST.EQ.2) THEN            !- output to PRINTER
        IOUT_D = CV(51)
        IF (IOUT_D.EQ.5) THEN             !- Postscript
C          .. postscript already set-up ?
        ELSEIF(IOUT_D.EQ.3) THEN           !- PCX file
c.. hmm where do we read in the size and the # colours ?
          CALL CREATE_SCREEN_BLOCK@ (ints(RES(1)),INTS(RES(2))
     +                             ,INTS(IBITPL),VSCREEN)
          CALL OPEN_VSCREEN@ (VSCREEN,IFAIL)        !- defer until later ?
        ENDIF
c.. open the appropriate output file
c.. and select the correct printer driver

      ENDIF

C----------- open virtual screen if buffered ---------------------------
c.. only if 'normal' VGA SCREEN ?
c.. hmm need to abstract the # of colour planes ?

c... This could also be handled within the 'Ouput Destination choice'
c... block of code above

      IF (IDEST.EQ.0.AND.IOP_BUF.NE.0) THEN
        CALL CREATE_SCREEN_BLOCK@ 
     +     (INTS(IWIN_XW),INTS(ABS(IWIN_YW)),INTS(IRESCP),VSCREEN)
        CALL OPEN_VSCREEN@ (VSCREEN,IFAIL) 
      ENDIF
C------------------------ blank window ---------------------------------
C... an alternative to this is to use a given PCX file as a backdrop
C... or a 'function' (eg. variable tints from top to bottom)
C... This MUST use the 'primitive' draw option !
c... hmmm .. really I want to use my *OWN* fill_rectangle routine !

      IF (CV(20).NE.-1) THEN
        IF (.NOT.(IDEST.EQ.2.AND.CV(51).EQ.8)) THEN ! NOT pen-plotters
          IF (CV(20).LT.1000) THEN       !- solid colour
            CALL MAKE_BOX (IWIN_XO,IWIN_YO,IWIN_XW,IWIN_YW,X,Y)
            CALL DR_PRIM (LINE,X,Y,4,CV(20),2)     !- as a filled box
          ELSE
C          .....  include code here for patterned backdrops
c          .....  eg. graded, grided,PCX-file, fractal? 
          ENDIF
        ENDIF
      ENDIF
C------------------------ draw the picture frame -----------------------
C.. this is WRONG as it does the whole page not just the current picture!
C.. carfull of DR_PRIM as well
c
c      IF (CV(21).NE.-1)   !- colour may also be the 'picture-# !':-)
c    +  CALL GRATIC (RES(1),RES(2),cv(22),cv(22), CV(21))

C--------------------- now plot the axes -------------------------------
C ... nice to draw the 'back-planes too'                             :-)
C ... ie. loop the '6' sides of the surrounding cuboid
C ... then b_p_c each, then loop each's 'squares' & draw wire-frame/fill
C ... to the nearest 10**size.. ie. square size =1mm/5m/20km etc.
c---> push this to a subroutine too !

      IF (CV(24).NE.-1) THEN
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
c**** need to use DR_PRIM to show the axes labels
c            CALL DRAW_TEXT@ (NTXT,SPT2(1),
c     +                            SPT2(2)),CV(27))
          ENDIF
        ENDDO
      ENDIF

C---------------- set min/max contour values etc. ----------------------
c.. only so we can 'spot' the min/max values found
c.. (so just initiate if this is the very first sub-facet ?
      C_R(1) =  1.e20  ! set the max/min contour values
      C_R(2) = -1.e20  ! to very high/low numbers

C------------------- loop facets and plot them -------------------------

      DO I_F = 1, NFCETS_D  ! =just the NFCETS that are to be drawn!

        IFC  = PL2 (PL1(I_F))         !- the facet's number
        IEL  = FACETS (1,IFC)         ! hence : el #
        IFACE= FACETS (2,IFC)         !       : &face dirn.
c.. also edge codes that are held in FACETS ?
            CALL GET_ELEMENT 
     +      (NUMS,INUMS,IEL, NUM, NOD,NDIME,ITYPE, IMAT,IU1,IU2,IU3)

        ICMAT= IC_MATS(MOD(IMAT-1,20)+1)  !  colour ?? (defer until later)
c       IF (IMAT.EQ.0) GOTO 81        !- skip this facet (.LE.0) ?
c                                     ! this should have happened earlier

        CALL GFACE (NOD,NDIME,ITYPE,IFACE,FS,NN_F,NN_FT, 2)  ! a facet's
        DO I=1,NN_FT                                         ! polygon
          FS(I) = NUM(FS(I))   ! **** Hmm this MUNGs FS() !
        ENDDO

C------------- get screen positions of the nodes of a facet ------------
c... Just accept the fact that I will have to call TRANSFORM :-)
      DO I=1,NN_FT    !--  hmm in a sub call could just pass 
        DO J=1,3      !--               SCOORD(J,FS(i) ?
          SC_F(I,J) = SCOORD (J,FS(I))
        ENDDO
      ENDDO


C-------------------- 'shrink' algorithm -------------------------------
c.. ignore the SC_F above and 'do_it_again' with scaling ?
c.. not quite the shrinking happens in screen-coords towards the element
c.. centroid (so screen Z shrinks too)
      IF (CV(31).ne.0) THEN
      DO J=1,3    ! no need for screen 'z' ?
        GPT(J) = 0.
        DO I=1,NOD
          GPT(J) = GPT(J) + GC (J,NUM(I))    !(last ref to NUMS)
        ENDDO
        GPT(J) = GPT(J) / REAL (NOD)     ! the element 'centroid'
      ENDDO
      CALL TRANSFORM (GPT,SPT, PM, IROT, XCEN,YCEN,SC_X,SC_Y,SC_Z)

      DO J=1,2                    ! no need for screen 'z' ?
        IF (CV(31).GT.0) then     ! as %-age shrinking
          DO I=1,NN_FT
            SC_F(I,J) = SC_F(I,J) - (SC_F(I,J)-SPT(J)) * cv(31)/100.
          ENDDO
        ELSEIF (CV(31).lt.0) THEN  ! as # of pixels to shrink by
          DO I=1,NN_FT
            SC_F(I,J) = SC_F(I,J) - SIGN (REAL(CV(31)),SC_F(I,J)-SPT(J))
          ENDDO
        ENDIF
      ENDDO
      ENDIF

C-------------------get disps of this facet-----------------------------
c.. I don't always need to do this ? .. maybe just create DISPS for all
      DO I=1,NN_FT  
        II  = FS(I)
        IB1 = NN * (ILD1-1)
        IB2 = NN * (ILD2-1)
        DO J=1,3
           C_F(I,J) = GC(J,II)
           IF (NLDS.GT.0) THEN  ! also should 'null' DI_F at the top :-)
             DI_F(I,J) = FAC1 *GDISPS(J,IB1+II) + FAC2 *GDISPS(J,IB2+II)
           ENDIF
        ENDDO
      ENDDO

C-------------------- 'undeformed' facet -------------------------------
c.. really the 'other' one .. if BASE is un-deformed then do deformed 
c.. here and vice-versa
c------ careful where I put this block of code
      IF (CV(16).NE.-1) THEN       !- this is the edges.. fill ?
        DO I=1,NN_F
          II = FS(I)
          DO J=1,3
            GPT(J) = GC(J,II)   ! start-point
          ENDDO
          CALL TRANSFORM (GPT,SPT, PM, IROT, XCEN,YCEN,SC_X,SC_Y,SC_Z)
          DO J=1,3
            SC_SF(I,J) = SPT(J)
          ENDDO
        ENDDO
        CALL DRAW_POLY (SC_SF,ISC,NN_F,CV(16),1) !  :-)
      ENDIF                        

C-----------------------------------------------------------------------
C------------------- select face/edge color options  -------------------
C-----------------------------------------------------------------------
c... maybe loop around this twice .. for 2 'layers' eg. conts+mats
c... careful with mean NN contours , contour NEN ?
c.. also chessboard MOD 3 or 4 ?? and nodal BW ?
                  
      DO IPASS = 1,2
        IF (IPASS.EQ.1) JJ = CV(1)
        IF (IPASS.EQ.2) JJ = CV(3)

        icol = -1      !- defualt is 'invis ???
        IF (JJ.EQ.-1) THEN                    ! See-thru
          ICOL = -1 
        ELSEIF (JJ.EQ.0) THEN                 ! Black (other colours ?) 
          IF (IPASS.EQ.1) ICOL = CV(2)
          IF (IPASS.EQ.2) ICOL = CV(4)
        ELSEIF (JJ.EQ.1) THEN                 ! Material # (indexed?)
          ICOL = ICMAT 
        ELSEIF (JJ.EQ.2) THEN                 ! Facet direction
          ICOL = ABS(IFACE) + 1
        ELSEIF (JJ.EQ.3) THEN                 ! Inside/Outside
          IF (C3.LT.0) ICOL = 10
          IF (C3.GT.0) ICOL =  3
        ELSEIF (JJ.EQ.4) THEN                 ! Glass front :-)
          IF (C3.LT.0) ICOL = -1                 ! only if no b_p_c ?
          IF (C3.GT.0) ICOL =  3
        ELSEIF (JJ.EQ.5) THEN                    ! Element #
          ICOL = INT(IEL / REAL(NEL) *CV(18)) + 2
        ELSEIF (JJ.EQ.6) THEN                    ! Mean Node #
          COL = 0
          DO I=1,NN_FT
            COL = COL + FS(I)
          ENDDO      
          ICOL = COL /REAL(NN_FT) /REAL(NN) * CV(18) + 3
        ELSEIF (JJ.EQ.7) THEN                    ! sub-facet #
          ICOL = 1+ I_SF
        ELSEIF (JJ.EQ.8) THEN                    ! Chessb.  sub-facets 
          ICOL = MOD(I_SF,2) + 2
        ELSEIF (JJ.EQ.9) THEN                    ! Chessboard elements
          ICOL = MOD (IEL,2) + 2

C---------------- 'flat' light-source shading of facet -----------------
c.. 10 is 'simple' ICOL= 2-->15, 11 is 'striped' with the material 
c.. colour , therefore only 'good' in 256 color mode so force OFF
C.. if there are insufficient colours ?

C---> should really compute c1,c2,c3 for EVERY facet & use above

        ELSEIF (JJ.GE.10.AND.JJ.LE.11) THEN

C... this is FLAT shading so use the screen coords in SC
C... and the 'projected-area' method  :-)
        C1 = 0. 
        C2 = 0.
        C3 = 0.
        DO I=1,NN_SF
          J  = MOD (I,NN_SF) + 1
          C1 = C1+(SC_SF(I,3)+SC_SF(J,3))*(SC_SF(J,2)-SC_SF(I,2)) /2.
          C2 = C2+(SC_SF(I,1)+SC_SF(J,1))*(SC_SF(J,3)-SC_SF(I,3)) /2.
          C3 = C3+(SC_SF(I,2)+SC_SF(J,2))*(SC_SF(J,1)-SC_SF(I,1)) /2.
        ENDDO
        C4 = ABS (C1*C1 + C2*C2 + C3*C3 )
        IF (C4.GT.1.E-12) THEN
        C4 = SQRT(C4)
          A = (LIGHT(1)*C1+LIGHT(2)*C2+LIGHT(3)*C3) /C4  ! light
          B = (                           (-1.)*C3) /C4  ! eye
          COL = SIGN(A,A*B)    ! ie. on the 'same' side as the light
c... hmm should all the colour range to be user-defined
c.. YUK use a function to to 0.->1. to a colour number
          ICOL = NINT(MIN(MAX(2.,2.+CV(18)*COL),CV(18)+1.))    ! clip
        ELSE
          ICOL = -1   ! zero-area face so skip colouring it
        ENDIF

c... use IMAT by 'striping' into 256 colour mode  (so use COL not ICOL!)
        IF (JJ.EQ.11) ICOL = (15-ICOL)*16 +IMAT+1 ! eg. palette #2

      ENDIF   ! end of facet flat-colour options

        IF (IPASS.EQ.1) ICOL_FACET = ICOL
        IF (IPASS.EQ.2) ICOL_EDGES = ICOL
      ENDDO

C------------------------ sub-faceting ---------------------------------
C... should call a subroutine to set 'up' this sub-facet
C... and another within the loop to return the local co-ords of one
C... of the points (given a 'patch' of what it looks like)

c.. surely if I am sub-faceting then I must calculate C1,C2,C3 for each !
c.. also colouring sub-facets' # etc. too.

      N_SF = CV(19)                   !- # along each edge (maybe= -1)
c     DO I_SF = 1, N_SF * N_SF        !- hence total #
      DO I_SF =    N_SF * N_SF,1,-1   !- hence total #

c...................................................................
C ... skip to main menu if a key pressed  ** mouse press ?? **
      EXIT = (KEY.ne.13.and.KEY_WAITING@()) 
      IF (EXIT) GOTO 1119
c...................................................................

C---------------------- no sub-faceting --------------------------------
c----- get the local co-ords of this facets' nodes using WTHATN
c.. ! note I am forcing 2D ! (am I really ??)
c... avoid WTHATN if 'OFF' style polygons ?? !
      IOP_CONT = CV(13)            !- ie. Are we contouring?
      IF (N_SF.eq.-1) THEN   
c----- Skip if OFF/NFF or if no contouring /disp. vectors/ regridding  -----
        IF (ITYPE.ne.9.AND.  
     +   (IOP_CONT.GE.1.or.CV(12).ge.0.or.cv(33).ge.0)) 
     +     CALL WTHATN (NN_F,2,1,LC_SF,ISC)
           NN_SF = NN_F
           DO I=1, NN_SF
           DO J=1,3                           ! (ie. get LC_SF and SC_SF)
             SC_SF(I,J) = SC_F(I,J)   
           ENDDO
         ENDDO

C--------------------------- sub-faceting ------------------------------
c.. so get local coords of this sub-facet (respect triangles)
c.. hence screen coords as FUN * SC_F -> SC_SF
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

      DO J=1,3                            !- 'close' the polygon
        SC_SF (NN_SF+1,J) = SC_SF (1,J)
      ENDDO


C----------------'flat'-colour in the (sub-) facet ---------------------

      IF (ICOL_FACET.NE.-1) 
     +     CALL DRAW_POLY (SC_SF,ISC,NN_SF,ICOL_FACET,2) 

C-----------------------------------------------------------------------
C----------- contouring of co-ord/disp/strain etc. ---------------------
C-----------------------------------------------------------------------

c... hmm light-source-shading should use the 'displaced' co-ords
c...  coord/disp/strains will use the 'original' geometry !

      N_CONTS = CV(18)
      IF (IOP_CONT.GE.14) THEN
C------------------------ light-source shading -------------------------
        IF (IOP_CONT.lt.20) THEN    !--- light-source shading

        DO I=1,NN_SF   !----- loop the nodes of a sub-facet (NN_SFT ?)

c................ from the record of the normals held in GDISPS .......
C... should never be 'sub-faceted'
         IF (IOP_CONT.EQ.18.or.iop_cont.eq.19) THEN
           VEC(1) = GDISPS(1,FS(I))
           VEC(2) = GDISPS(2,FS(I))
           VEC(3) = GDISPS(3,FS(I))
c          print*,vec(1),vec(2),vec(3),
c    +           sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
           T = 1.    !- (Already normalised)
c................ or use the shape functions at each node ............
         ELSE
           CALL SAMPLE (SC_F,ISC,DI_F,ISC,NN_FT,2,3,
     +                VAL,VEC,LC_SF,ISC,I,IOP_CONT)
           T = SQRT(VEC(1)**2+VEC(2)**2+VEC(3)**2) ! noramalise face normal
         ENDIF
c......................................................................
         IF (T.GT.0.) THEN
           A = (LIGHT(1)*VEC(1)+LIGHT(2)*VEC(2)+LIGHT(3)*VEC(3)) /T
           B = (                                   (-1.)*VEC(3)) /T  !eye
           B = 1.     !(try this ?)
           COL = SIGN(A,A*B)        ! ie. on the 'same' side as the light
C          ACOL(I) = MIN (MAX(2.,(2.+CV(18)*COL)),cv(18)+2.) ! clip to 2.--> 15 ?
           ACOL(I) = MIN (MAX(0.,(0.+CV(18)*COL)),cv(18)+0.) ! clip to 2.--> 15 ?
        ELSE
          ACOL(I) = 9   ! if zero-area draw in grey ?(col=9) NO -1 !
        ENDIF
        IF (IOP_CONT.EQ.15)  ACOL(I) = ACOL(I) + (IMAT-1) * 16 
        IF (IOP_CONT.EQ.19)  ACOL(I) = ACOL(I) + (IMAT-1) * 16 
      ENDDO
C------------------ contouring of strains, disps, etc. -----------------
c.. hmm contour Node Number also ?
c... can't this call to SAMPLE be built into the previous
c... I can still find the min/max but if shading SKIP the scale to0.->1. !

      ELSE
        DO I=1,NN_SF   !----- loop the nodes of a sub-facet (NN_SFT ?)
      
        CALL SAMPLE (C_F,ISC,DI_F,ISC,NN_FT,2,3,
     +               VAL,VEC,LC_SF,ISC,I,IOP_CONT)
        DIS = VAL
        C_R (1) = MIN (C_R(1),DIS)         !  min contour value found
        C_R (2) = MAX (C_R(2),DIS)         !  max contour value found

        IF (C_R(6).GT.0) DIS = ABS(DIS)    ! take absolute if desired
        DIS = (DIS-C_R(3))/(C_R(4)-C_R(3)) ! scale to 0.--> 1.

c....logarithic scaling  .... (skip if 'log factor' = 1 for speed)
        IF (INT(C_R(5)*10).NE.10)  DIS = DIS ** (1./c_r(5)) 


c.. so valid values to contour are 0. --> N_CONTS
C.. but values (not to contour) will lie outside these limits
c.. so where to they get 'clipped' off ?

c       ACOL(I) = 2 + DIS * (N_CONTS)     !- don't forget the offset !
        ACOL(I) = 0 + DIS * (N_CONTS)     !- don't forget the offset !
      ENDDO       
      ENDIF   ! end of contouring / lighting

      CALL GOURAUD (SC_SF,ISC,ACOL, N_CONTS,NN_SF,CV(15),1,ICMAT)  !-- filled
      CALL GOURAUD (SC_SF,ISC,ACOL, N_CONTS,NN_SF,CV(14),2,ICMAT)  !-- +lines

      ENDIF   ! endif of an 'interpolated' colour scheme

C-------------------- * REGRIDDING * -----------------------------------
      IF (CV(33).GE.0) THEN
        DO JJ=1,NDIM
          DO I=1,NN_SF   !----- loop the nodes of a sub-facet (NN_SFT ?)
            CALL SAMPLE (C_F,ISC, DI_F,ISC, NN_FT,2,3,     !- get coord
     +               VAL,VEC, LC_SF,ISC, I,21)  !(or use 20)
            VAL = VEC(JJ)
            VAL = VAL - GC(JJ,NN+1)     !- push to 0.--> tot.width
            VAL = VAL / DIAG2               !- scale  to 0. -> 1.
            n_c = 2**cv(34)
            VAL = VAL * n_c                 !- expand to 0.--> nlines

            ACOL(I) = VAL
          ENDDO       
          ICOL = CV(33)      !- hmm .. will act as a 'code' to GOURAUD !

          CALL GOURAUD (SC_SF,ISC,acol,n_c,NN_SF,icol,2,ICMAT)  !(as lines)

        ENDDO    ! loop x->y->z
      ENDIF  ! end of regridding
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-------------------- displacement vectors -----------------------------
c.. easy ! just tranform BOTH start and end-points
c.. careful of these anisotropic scale factors (ignore!)

c.. I should also be able to draw disp vector with sub-facets
c.. It will be the same as this .. just move the code into the facets' 
c     loop .. (but need to sample for the disps themselves)

c   Here we transform BOTH the start and end-points ..
c   however by (sub-)faceting we already know the end-point :-)
c
c

      IF (CV(12).NE.-1) THEN
        DO I=1,NN_SF
c. . .. . . . . .. . . .. . . 
          CALL SAMPLE (C_F,ISC, DI_F,ISC, NN_FT,2,3,     !- get coord
     +               VAL,VEC, LC_SF,ISC, I,20)
          DO J=1,3
            GPT(J) = VEC(J)
          ENDDO
          CALL TRANSFORM (GPT,SPT, PM, IROT, XCEN,YCEN,SC_X,SC_Y,SC_Z)
c...................
          CALL SAMPLE (C_F,ISC, DI_F,ISC, NN_FT,2,3,     !- get disp
     +               VAL,VEC, LC_SF,ISC, I,30)
          DO J=1,3
            GPT(J) = GPT(J) + VEC(J) * FACT
          ENDDO
          CALL TRANSFORM (GPT,SPT2, PM, IROT, XCEN,YCEN,SC_X,SC_Y,SC_Z)
c. . .. . . . . .. . . .. . . 
          CALL DR_ARROW  (SPT,SPT2,10,20,CV(12))
        ENDDO
      ENDIF
c... This is the end of the face-filling so maybe do all together here ?      


C-------------------- edge the sub-facet -------------------------------
      ICOL = cv(5)
      IF (ICOL.NE.-1) CALL DRAW_POLY (SC_SF,ISC,NN_SF,ICOL,1) ! edges

C.. OK have a CV(6) (say) for edging the 'exterior' only ?
c.. This would have the effect of drawing 'smooth' element edges
c.. if quads need the NXE/NYE of the sub facet
c.. so draw LH edge if NXE=1.. etc., etc. (easy)
c.. in triangle it is a little harder BOT if 1-> NSF, RH if 1,3,5 ???  

C------------------ end of the sub-faceting ----------------------------
c  82 CONTINUE    ! skip-to index for backface sub-facet cull (icol=-1)
      ENDDO
C---------------------- edge facet -------------------------------------
C... need to sub-sample for 'smooth' edges 
C...  -also identify the geometry/material 'edges'
c... If I have done the 'smooth' edges of a sub-facet then there is no
c... need to do any more ?
c... BUT I could have a second N_SUB_FAC & loop and sample along the edges
C... to produce a l-o-n-g polyline to draw :-)

c.. hmm this is 'easy' ! .. just draw edge #1 if LHS of element etc. :-)
c.. so need a CV code to edge all of a sub-facet or just the boundary
c.. or better: a second sub-facet edge colour code for 'boundary' only

      IF (ICOL_EDGES.NE.-1) 
     +     CALL DRAW_POLY (SC_F,ISC,NN_F,ICOL_EDGES,1) ! edges

c--------------------- 'coded' facet edgeing ---------------------------
c.. Loop facet edges (3 or 4)
c.. Extract its code from FACETS (4-5-6-7)
c.. If code = 1 (material boundary) etc.
c.. Extract its colour code from CV(50+) then draw ths edge in
c--( Also there will be a code to use sub-sampling to get smoother edges)
c--( This will use the shape funs of a sinle edge (3nl element ?)



C----------------------- draw nodes ------------------------------------
c.. should use my primitive draw-er (hence PS)
c.. at least for now .. skip this if we are doing PS :-) ?
      IF (CV(9).NE.-1) THEN
        RADIUS = CV(9)
        II = NINT (CV(10)*RES(1)/640.)  !-- for 'large pcx/plt's etc.
        DO I=1,NN_FT
          X(1) = SC_F(I,1)
          Y(1) = SC_F(I,2)
          X(2) = NINT (CV(10)*RES(1)/640.)        !- radius
          Y(2) = X(2)                  !- let y-radius=x-radius
          CALL DR_PRIM (LINE,X,Y,2,CV(9),10)          !- (code 10=circle)
        ENDDO
      ENDIF
C-------------------- node numbering -----------------------------------
C.... need an option to draw 'centred' on the node?

      IF (CV(11).NE.-1) THEN
        DO I=1,NN_FT       ! ie. 'all' the nodes on this facet
          WRITE (LINE,'(I4)') FS(I)
c         CALL TRIM (LINE)           ! ie. left-justify ?
          X(1) = SC_F(I,1)
          Y(1) = SC_F(I,2)
          CALL DR_PRIM (LINE,X,Y,1,CV(11),20)
c         CALL DRAW_TEXT@ (NTXT,INTS(SC_F(I,1))
c    +                         ,INTS(SC_F(I,2)),ints(CV(11)))
        ENDDO
      ENDIF

C------------------element numbering -----------------------------------
      IF (CV(8).NE.-1) THEN
        WRITE (LINE,'(I4)') IEL
c       CALL TRIM (LINE)            !       right-justify (or centre?)
        X(1) = 0.
        Y(1) = 0.                   ! find the centroid of the face
        DO J=1,NN_F                 ! (or NN_FT ?)
          X(1) = X(1) + SC_F(J,1)
          Y(1) = Y(1) + SC_F(J,2)
        ENDDO
        X(1) = X(1) /REAL(NN_F)
        Y(1) = Y(1) /REAL(NN_F)
        CALL DR_PRIM (LINE,X,Y,1,CV(8),20)

c        CALL DRAW_TEXT@ (NTXT,INTS(IXP/NN_F),INTS(IYP/NN_F)
c     +                  ,INTS(CV(8)))
      ENDIF
C------------------ material numbering -----------------------------------
      IF (CV(7).NE.-1) THEN
        WRITE (LINE,'(I4)') IMAT
c       CALL TRIM (LINE)            !       right-justify (or centre?)
        X(1) = 0.
        Y(1) = 0.                   ! find the centroid of the face
        DO J=1,NN_F                 ! (or NN_FT ?)
          X(1) = X(1) + SC_F(J,1)
          Y(1) = Y(1) + SC_F(J,2)
        ENDDO
        X(1) = X(1) /REAL(NN_F)
        Y(1) = Y(1) /REAL(NN_F)
        CALL DR_PRIM (LINE,X,Y,1,CV(7),20)
      ENDIF
C-------------------- displacement vectors -----------------------------
C (MOVED TO INSIDE SUB-FACET LOOP ) 29-10-93
C---------------------- end of this facet ------------------------------
c  81 CONTINUE
      ENDDO    

c--------------------- if no mesh to draw  -----------------------------
      IF (NEL.EQ.0) THEN
        LINE = 'NO MESH'
        X(1) = IWIN_XO - IWIN_XW /2
        Y(1) = IWIN_YO - IWIN_YW /2
        CALL DR_PRIM (LINE,X,Y,1,13,20)
      ENDIF
C--------- Draw the 'extra-lines (eg. model boundaries) ----------------
c... really a subroutine .. maybe can do as a list of nodes as well
c... 'cos' this is easier to store (hence flag end-of chains with -1 ?
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

        CALL DR_PRIM (LINE,X,Y,2,1,1)       !--- draw_a_line
      ENDDO
C--------------- finished model so draw border etc.----------------------
C... need to invert for the plotter (its origin is the bottom left)
c.. I don't think so .. just gets drawn upside-down which is OK
c--> also if buffered then there is NO IWIN_XO to add !
c  (( or buffer-> screen   before drawing the frames etc. ? ))
c
c...  Should really have *Two* borders : one around each image
c...  and one (maybe thicker ?) around the whole picture
c

C... ok loop and draw around the current window
      IF (CV(28).GE.0) THEN
        CALL MAKE_BOX (IWIN_XO,IWIN_YO,IWIN_XW,IWIN_YW,X,Y)
        CALL DR_PRIM (LINE,X,Y,4,CV(28),1)     !--- as a closed polygon
      ENDIF

c-------------------- show the load step number ------------------------
      IF (CV(39).GT.-1) THEN
        WRITE (LINE,'(A,F8.2)') 'LS =',LD
        X(1) = RES(4)               !.. in top LH corner
        Y(1) = RES(5) + 20.  
        CALL DR_PRIM (LINE,X,Y,1,CV(39),20)
      ENDIF
c-------------- show the (animation) image number ----------------------
      IF (CV(40).GT.-1) THEN
        WRITE (LINE,'(A,I4)') 'PIC=',IPIC_N
        X(1) = RES(4)               !.. in top LH corner
        Y(1) = RES(5) + 40.  
        CALL DR_PRIM (LINE,X,Y,1,CV(40),20)
      ENDIF
C-----------------------------------------------------------------------
c... jump-to point if we are 'aborting' a picture-draw

 1119 IF (IOP_BUF.EQ.1)THEN        !-- only if buffered ?
        IF (.NOT.EXIT)             ! inhibit if only partly drawn ?
     +     CALL VSCREEN_TO_SCREEN@ 
     +       (INTS(iwin_xo),INTS(iwin_yo),INTS(0),IFAIL) 
        CALL CLOSE_VSCREEN@()      ! do I need to do this ? ,_ yes ! -why? 
        CALL RETURN_STORAGE@ (VSCREEN)  ! remove it completely -why ?
        IF (IFAIL.NE.0) PRINT*,'Ifail=',IFAIL,' VS_2_S' 
      ENDIF
      IF (EXIT) GOTO 1110   ! jump back to the top of the program

C-----------------------------------------------------------------------
C--------------------- auto-'draw-next image' options ------------------
C-----------------------------------------------------------------------
c.. can do at the begining of 'draw_page' as if 1st pass .. rest valuss
c.. else update values 
c.. (AT end .. if multi.. jump back if not yet last)  

      IF (CV(37).GT.0) THEN
C       IPIC  = CV(41)
C       NPICS = CV(42) * CV(43)
        CV(41) = CV(41) + 1
      IF (CV(37).EQ.1) THEN  !-------------- auto_load_step ------------
        LD = LD+1                              ! load step
        IF (CV(41).LE.CV(42)*CV(43) ) GOTO 901 ! redraw
        CV(41) = 1                             !- reset pic #
        LD = 1                  !  reset to beginning ready for a re-draw
      ELSEIF (CV(37).eq.2) THEN !----------- auto_multi_view -----------
        CV(38) = cv(38) +1                     ! view direction
        IF (CV(41).le.CV(42)*CV(43) ) GOTO 901 ! redraw
        CV(41) = 1                             ! reset pic #
        CV(38) = 1  ! reset to the first picture ready for a re-draw
      ELSE
        PRINT*,' Multi-view type', CV(37) ,' is unknown'
      ENDIF
      ENDIF

C--------------------- end of image drawing ----------------------------
c.. now close any output devices and rest to VGA mode

      CONTINUE   ! 'exit'- label (if a key/mouse-button pressed) ?? 
      CALL CLOCK@ (TIME(4))

C-------------  dump image to printer/plotter file etc.-----------------
      IF (IDEST.EQ.0) THEN    !---- if VGA ------
c... Ok redraw pictures frames if obscured

        IF (L_OVERLAP) CALL POST_MENUS (1, 0,0,Ckey,IVAL) !- if obscured

      ELSEIF (IDEST.EQ.1) THEN    !-------- SVGA prevew ----------
        READ*                    !- wait for a keypress
        CALL INTO_GRAPHICS ()      !- re-initialise VGA ?
        CALL SET_ALL_DACS()        !- and use the old palette
        L_OVERLAP = .TRUE.       !- will cause menus to be redrawn

      ELSEIF (IDEST.EQ.2) THEN    !- Printer output

      IOUT_D = CV(51)
      CALL SOUND@ (1024,1)
      IF (IOUT_D.EQ.3) THEN         !------ PCX output to file ------
        CALL VSCREEN_TO_PCX@ (PICFILE,IFAIL)
        CALL PCX_FIX (PICFILE,PAL) 
c.. check this PCX_PALETTE fixing .. is it not yet automatic
c.. and update to allow for 256 colour palettes
        CALL CLOSE_VSCREEN@ ()
        CALL RETURN_STORAGE@ (VSCREEN)  !- why ?

      ELSEIF (IOUT_D.EQ.4) THEN        !------- HP-LJ to file ----------
        CALL CLOSE_GRAPHICS_PRINTER@ () 
        IF (PICFILE.eq.'LASER') THEN      ! auto print to HP-UX :-)
          CALL CISSUE@ ('pr_hpux',ifail)    !( sends & deletes 'LASER')
        ENDIF         
        IF (PICFILE.eq.'PRINT') THEN     
          CALL CISSUE@ ('copy/b print prn',ifail) ! prints
        ENDIF

      ELSEIF (IOUT_D.EQ.5) THEN        !--------- PostScript ------------
        IF (INDEX(PICFILE,'.EPS').GT.0) THEN
           CALL CLOSE_PS_PRINTER () 
C...       nothing else to do ??
        ELSEIF (PICFILE.eq.'PRINT') THEN            ! (auto print to HP-UX ??)
          CALL CLOSE_PS_PRINTER () 
          CALL CISSUE@ ('copy print prn',ifail)     ! copy to printer
          CALL CISSUE@ ('echo showpage >prn',ifail) ! and show the page
          CALL CISSUE@ ('erase print',ifail)         ! then delete
        ELSE
          WRITE (55,'(A)')  'showpage'             ! add the showpage 
          CALL CLOSE_PS_PRINTER () 
        ENDIF

      ELSEIF (IOUT_D.EQ.6) THEN        !------ pro-printer to LPT1 -----
        CALL CLOSE_GRAPHICS_PRINTER@() 

      ELSEIF (IOUT_D.EQ.8) THEN        !------- HP7550 to file ---------
        CALL CLOSE_PLOTTER@ ()            ! HP7550 to plotter ?

      ENDIF   !- end of printer choices

      CALL SOUND@ (1024,1)
c     CALL RESET_SIZES (RES)           !- in the 'normal' window
c     IRESC  = 16
c     IRESCP =  4
c     AR     = -1.       !--- back to noraml

      ENDIF              !- end of screen/SVGA/printer outputs
      IDEST = 0          !- always revert back to the normal VGA mode
      CALL CLOCK@ (TIME(5))   !- time at the end of ALL drawing

C-------------------------- animation ----------------------------------
c.. this is whereby the curent picture is saved (PCX/memory)
c.. then some viewing values are updated and A Force is
c.. made to loop-back to redraw the image
c.. maybe ONLY allow as Macros.. SPIN is easy but 'load-steps' ?

      IF (IOP_ANIM(1).GT.0) THEN
                        
c---------- save the image if desired (memory/PCX) ---------------------
        IF (IOP_ANIM(7).EQ.0) THEN
          CONTINUE                           ! don't save
        ELSE
        IOP = MOD (IOP_ANIM(7),10)           ! the 'bit-left over'

c         IF (IPIC_N.GT.IOP_ANIM(9)) IOP_ANIM(1) = 0 ! kill at the end :-)
c         WRITE (*,'(F4.2)') LD              !- load-step !
c          IF (IOP.EQ.1) CALL GET_SCREEN_BLOCK@      !-- 320x200 'window' 
c     +       (IWIN_XO+0,IWIN_YO+0,IWIN_XO+320-1,IWIN_YO+200-1,BUFFER)
c          IF (IOP.EQ.3) CALL GET_SCREEN_BLOCK@      !-- 'square'
c     +       (0,0,480-1,480-1,BUFFER)

c-----------------------------------------------------------------------
          IF (IOP.EQ.2) CALL GET_SCREEN_BLOCK@      !-- all VGA
     +   (0,0,640-1,480-1,BUFFER)
          IF(IOP.EQ.4) CALL GET_SCREEN_BLOCK@       !-- just window
     +   (IWIN_XO,IWIN_YO,
     +    IWIN_XO+IWIN_XW-1,IWIN_YO+IWIN_YW-1,BUFFER)
c-----------------------------------------------------------------------

c       PRINT*,C_T//'picture #=',IPIC_N        !- frame number (status bar?)

        IF (IOP_ANIM(7).LT.10) THEN  !--- dump to a PCX file
          WRITE (PICFILE,'(a,i3.3,a)') !- generate the file_name
     +           FILE_ROOT(1:leng(FILE_ROOT)),ipic_n,'.pcx'
c         print*,' saving image #', ipic_n
          CALL SCREEN_BLOCK_TO_PCX@ (PICFILE,BUFFER,IFAIL)
          CALL DOSERR@(IFAIL)    !- If it doesn't work
          CALL PCX_FIX (PICFILE,PAL)
          CALL RETURN_STORAGE@ (BUFFER)     !- is this necessary ? <- YES !
        ELSE
          IF (IMAGE(IPIC_N).NE.-1) 
     +    CALL RETURN_STORAGE@ (IMAGE(IPIC_N))   !- zap any old image
          IMAGE (IPIC_N) = BUFFER                !- remember this frame
        ENDIF
        ENDIF             ! end of 'save-image-to-pcx-or-memory'

C------------------- Update the viewing parameters ---------------------

        I = IOP_ANIM(2)     ! = # of frames to reach final load increment
        IF (I.NE.0) THEN
          LD = LD + LDD                              ! increment the load
          IF (I.GT.0.AND.LD.GT.NLDS) LD = 1.         !---- wrap
          IF (I.LT.0.AND.LD.GE.NLDS) LDD = -ABS(LDD) !---- bounce down! 
          IF (I.LT.0.AND.LD.LE.1   ) LDD =  ABS(LDD) !---- bounce up! 
        endif
c.... also sinusoidal animation of a single load step for eigenmodes ....
        XME = XME + IOP_ANIM(3)  !------- move the EYE/LIGHT -----------
        YME = YME + IOP_ANIM(4)
        XML = XML + IOP_ANIM(5)  ! why were these once halved ??
        YML = YML + IOP_ANIM(6)

        IPIC_N = IPIC_N +1             ! increment pic # for next time
        IF (IPIC_N.GT.IOP_ANIM(9)) IOP_ANIM(1) = 0  !- so kill
        IF (IOP_ANIM(1).NE.0) GOTO 901  !-- straight to re-draw ?
      ENDIF

C-------> goto the start of the loop --> --> --> --> --> --> --> --> -->
      GOTO 1110
C--------------------> E X I T   P R O G R A M <------------------------
c....... eg. from various errors or if event = 'q' 
  999 CONTINUE  
      CALL TEXT_MODE@()
      END                         !- The end of the Main program

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c-----------------------------------------------------------------------
      SUBROUTINE SET_CV (ITEM,IVAL,IOP)
C
C     This database holds the values of CV .. the 'ITEM COLOURS' etc.
C     if IOP = 0  : all values are initialised
C        IOP = 1  : CV(ITEM) is set to IVAL
C        IOP = 2  : the current ITEM is rememberd
C        IOP = 3  : CV(rem.ITEM) is set to IVAL
C  (     IOP = 4  : CV(rem.ITEM  is re-initialised ? )
C
C        IOP = 5  : CV(ITEM) is returned in IVAL
C        IOP = 6  :  The remembered item is returned
C        IOP = .. :  any other options ?
C
C---------------- colour info vector -----------------------------------
      SAVE
      COMMON  /CV/CV

      PARAMETER (   M_CV = 60)
      INTEGER  CV  (M_CV)         !- the data table itself
     +        ,CV1 (M_CV)         !- the initial copy of CV

      DATA ITEM_STORE/1/          !- default stored 'item'
      DATA (CV1(I),I=1,10) /
     +   1    !  1.   face colour option 1=mats, 2=facet dirn, etc.
     +,  0    !  2.   actual face colour if CV(1) option = 0 ?Mat table!

     +,  0    !  3.   edge colour option
     +,  1    !  4.   actual edge colour if CV(3) option = 0 ?Mat table!

     +, -1    !  5.   col: sub-edges
     +, -1    !  6.   col:   ( material edges )
     +, -1    !  7.   col: material number
     +, -1    !  8.   col: element number
     +, -1    !  9.  col:  nodes
     +,  2 /  ! 10.    node size (pixels or 1/640ths. of image x ?)
      DATA (CV1(I),I=11,20) /
     +  -1    ! 11.  col:  node # 
     +, -1    ! 12.  col:  disp vector 
     +, -1    ! 13.  *** contour code ***\
     +, -1    ! 14.  *  contour line type -1,0,1,(2=coloured)
     +,  2    ! 15.  *  contour fill type -1,0,1,(2=col,7/8=zebra)
     +, -1    ! 16.   undeformed facet colour     ! (obs) ?
     +,  2    ! 17.   depth sort (0=off,1=nearest,2=furthest,3=central)
     +, 13    ! 18.  ** number of colours to use for contouring **
     +, -1    ! 19.  * # sub-facets per element (per side)
     +,  0 /  ! 20.  * screen background colour
      DATA (CV1(I),I=21,30) /
     +  -1    ! 21. col: graticule
     +,  1    ! 22.      # x-graticule lines
     +,  1    ! 23.      # y-graticule lines (?) 
     +, -1    ! 24. col: axes
     +, -1    ! 25.       ( axes lengths )
     +, -1    ! 26.       ( flag fo -ve axes too )
     +, -1    ! 27. col:  axes labels (x,y,z)
     +,  1    ! 28. col:  window border
     +, -1    ! 29.       ( redraw . stamp colour)  <--- what ??
     +, -1 /  ! 30.
      DATA (CV1(I),I=31,40) /
     +   0    ! 31.    ! shrink factor (+ve=%-age,-ve = in pixels)
     +,  1    ! 32.    ! backface-poly cull 0=none,1=back,2=front
     +, -1    ! 33.    ! colour of regridded lines     ***
     +,  4    ! 34.    ! # regridded lines = 2**n      *** (eg = 16)
     +, -1    ! 35.    !
     +, -1    ! 36.    ! 
     +, -1    ! 37.    ! multi-view type (-1=none,1=load#,2=view)
     +, -1    ! 38.    ! preset view number (-1 = via XME,YME) (obs?)
     +, -1    ! 39.    !  col:  load step #  
     +, -1 /  ! 40.    !  col: image # (when animating)
      DATA (CV1(I),I=41,50) /       ! windowing etc.
     +   1    ! 41.    ! picture #
     +,  1    ! 42.    ! # pics. in x-direction
     +,  1    ! 43.    ! # pics. in y-direction           
     +,  0    ! 44.    ! image size - x              = res(1)    ?
     +,  0    ! 45.    !            - y              = res(2)    ?
     +,  0    ! 46.    !            - offset in x    = res(4)    ?
     +,  0    ! 47.    !            - offset in y    = res(5)    ?
     +,  0    ! 48.    !%image centre offset - x  (ie. 'panning' )
     +,  0    ! 49.    !                     - y 
     +,  5 /  ! 50.    ! SVGA mode number (1024x768x256?) to use

      DATA (CV1(I),I=51,60) /       ! etc.
     +   5    ! 51.    ! Printer device to use (eg. 5=PS etc.)
     +,  1    ! 52.    !
     +,  1    ! 53.    !
     +,  0    ! 54.    !
     +,  0    ! 55.    !
     +,  0    ! 56.    !
     +,  0    ! 57.    !
     +,  0    ! 58.    !
     +,  0    ! 59.    !
     +,  5 /  ! 60.    !

c-----------------------------------------------------------------------

      IF (IOP.EQ.0) THEN           !-------- reset data ---------
        DO I=1,M_CV
          CV(I) = CV1(I)
        ENDDO
      ELSEIF (IOP.EQ.1) THEN       !-------- set data -----------
        CV (ITEM) = IVAL
      ELSEIF (IOP.EQ.2) THEN
        ITEM_STORE = ITEM
      ELSEIF (IOP.EQ.3) THEN
        CV (ITEM_STORE) = IVAL
      ELSEIF (IOP.EQ.4) THEN
        CV (ITEM_STORE) = CV1(ITEM_STORE)

      ELSEIF (IOP.EQ.5) THEN       !------- extract data --------
        IVAL = CV (ITEM)

      ELSEIF (IOP.EQ.6) THEN       !--- return 'remembered'
        ITEM = ITEM_STORE
      ELSE                         !---------- error ------------
        CALL MYERROR (2,'Unknown SET_CV option')
      ENDIF

      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE CONNECTIVITY_IMAT (GC,IGC,NN,NDIM,NUMS,INUMS,NEL,P)
C
C     This will find all elements that are interconnected and give them
C     the same material property number.
C     on return the IMAT column of NUMS is altered, P is workspace but 
C     also returns the node group numbers (hence orphans = -1)
C
      REAL       GC (IGC,NN)         !-
      INTEGER  NUMS (INUMS,NEL)      !-
     +         ,NUM (27)             !-
     +           ,P (NN)             !- workspace
                                  
C.............. initialise ................
          DO I=1,NN
            P(I) = -1       !- zap all nodes' codes
          ENDDO
          DO IEL=1,NEL
            CALL PUT_EL_IMAT (NUMS,INUMS,IEL, -1)  !- zap all Mats
          ENDDO
c.......... loop and impliment
          IMAT = 0
          DO IEL2=1,NEL
            CALL GET_ELEMENT  (NUMS,INUMS,IEL2, NUM, 
     +         NOD_A,NDIME,ITYPE, IMAT_a,IU1,IU2,IU3)
            IF (IMAT_a.EQ.-1) THEN    !- only if this not yet done
              imat = imat + 1                  
              print*, 'doing material #', imat
              CALL PUT_EL_IMAT (NUMS,INUMS,IEL2, IMAT)
              DO J=1,NOD_a               !- set all these nodes as 'done'
                P(NUM(J)) = IMAT       !-
              ENDDO                      !-
c..................
c.. hmm should iterate this until nothing further changes
 1024         idone = 0
              DO IEL=IEL2+1,NEL              !- Loop the rest
                CALL GET_ELEMENT  (NUMS,INUMS,IEL, NUM, 
     +            NOD,NDIME,ITYPE, IMAT_B,IU1,IU2,IU3)
c               NOD = NUMS(1,IEL)
                IF (IMAT_B.eq.-1) then     !- only if not yet done ?
                  DO J=1,NOD               !- check if 'any' hit
                    IF (P(NUM(J)).eq.IMAT) then  
                      CALL PUT_EL_IMAT (NUMS,INUMS,IEL, IMAT)
                      DO JJ=1,NOD              !- and all its nodes
                        P(NUM(JJ)) = IMAT    !- 'on' too
                      ENDDO                    !-
                      idone = 1      !- flag that an 'hit' has happened
                      GOTO 1023
                    ENDIF
                  ENDDO
 1023             CONTINUE    !- jump out
                ENDIF
              ENDDO     !- loop the rest
              IF (idone.eq.1) goto 1024       !- re-iterate
c..................

            ENDIF   !- skip already 'done' elements
          ENDDO
          print*,' There are ',imat,' different zones in the mesh!'
          RETURN
          END
C-----------------------------------------------------------------------
        SUBROUTINE MIRROR_ELEMS (GC,IGC,NN,NDIM,
     +                              NUMS,INUMS,NEL,IDIR,VAL)
C
C       This mirrors a mesh about a given plane
C        ie. where IVAL (x/y/z) = VAL
C
      REAL       GC (IGC,*)          !-
      INTEGER  NUMS (INUMS,*)        !-
     +         ,NUM (27)             !-
c    +           ,P (NN)             !- workspace
c..... be careful about the size of the mesh as we mirror it ...

c.......................... first copy the nodes ....................... 
        DO I=1,NN  
          DO J=1,NDIM
            IF (J.NE.IDIR) THEN
              GC (J,NN+I) = GC(J,I)        ! copy
            ELSE
              GC (J,NN+I) = 2*VAL - GC(J,I)  ! reflect 
            ENDIF
          ENDDO     
        ENDDO

c....................... then copy the elements ........................
        DO IEL=1,NEL 
          CALL GET_ELEMENT 
     +      (NUMS,INUMS,IEL, NUM, NOD,NDIME,ITYPE, IMAT,IU1,IU2,IU3)
          DO J=1,NOD
            NUM(J) = NUM(J) + NN 
          ENDDO
          CALL PUT_ELEMENT 
     +      (NUMS,INUMS,IEL+NEL, NUM, NOD,NDIME,ITYPE, IMAT,IU1,IU2,IU3)
        ENDDO

        NN     = NN  * 2
        NEL    = NEL * 2
        RETURN
        END
C-----------------------------------------------------------------------
      SUBROUTINE GET_NORMALS (GC,IGC,NDIM,NN,NUMS,INUMS,NEL, 
     +                           GDISPS,MDF,P)
c
c     This find the mean normal direction at each node
c     and stores it in the GDISPS table (as the first load-case)
c     This is only really useful for 3D polygonal meshes (eg. OFF format)
c         Dan Kidger   18-11-93
c
C ... I *could* avoid P() by using the 3rd column of GDISPS ??
      INTEGER NUMS (INUMS,*)      ! elements
     +          ,P (*)            ! # of elems per node (workspace)
     +        ,NUM (27)           ! nodes of an element
     +         ,FS (27)           ! nodes on a facet

      REAL      GC (IGC,*)        ! nodal coords
     +     ,GDISPS (MDF,*)        ! BIG table of displacements
c------------------------- null the arrays ----------------------------
      DO I=1,NN
        DO J=1,NDIM
          GDISPS(J,I) = 0.
        ENDDO
      ENDDO
      DO I=1,NN
        P(I) = 0       !- to count the number of facets at a node
      ENDDO
c----------------- loop the elements and build the facets --------------
      DO IEL=1,NEL
        IFACE = 1          !.. cos always only 1 facet per elem
        CALL GET_ELEMENT   (NUMS,INUMS,IEL, NUM, NOD,NDIME,ITYPE,
     +                                            IMAT,IU1,IU2,IU3)
        CALL GFACE (NOD,NDIME,ITYPE,IFACE,FS,NN_F,NN_FT, 2)
c........ here FS() is always just 1,2,3,4..NOD !
      DO I=1,NN_FT                    
        FS(I) = NUM(FS(I))       !- ** careful if we MUNG FS() **
c        DO J=1,3
c          SC_SF(I,J) = GC(J,FS(I))
c        ENDDO
      ENDDO
C ( Note. here SC_SF is *NOT* screen coords but just model coords )
        C1 = 0.
        C2 = 0.
        C3 = 0.
        DO I=1,NN_FT
          K  = FS (I)                  !- node #
          J  = FS ( MOD(I,NN_F) + 1)   !- node #+1
          C1 = C1 + (GC(3,K)+GC(3,J)) * (GC(2,J)-GC(2,K)) /2.
          C2 = C2 + (GC(1,K)+GC(1,J)) * (GC(3,J)-GC(3,K)) /2.
          C3 = C3 + (GC(2,K)+GC(2,J)) * (GC(1,J)-GC(1,K)) /2.
        ENDDO
        C4 = ABS (C1*C1 + C2*C2 + C3*C3 )
        IF (C4.GT.1.E-12) THEN      !- skip Zero-area-elems :-)
          C4 = SQRT(C4)
          DO I=1,NN_FT
            K = FS(I)               !- the actual node number
            GDISPS(1,K) = GDISPS(1,K) + C1 /C4
            GDISPS(2,K) = GDISPS(2,K) + C2 /C4
            GDISPS(3,K) = GDISPS(3,K) + C3 /C4
            P(K) = P(K) + 1      !- incement # of faces at this node
          ENDDO
        ENDIF
      ENDDO !- loop elements
C-------------------- now average the values --------------------------
      DO I=1,NN
        IF (P(I).GT.0) THEN        !- skip 'orphan' nodes
          DO J=1,3
            GDISPS(J,I) = GDISPS(J,I) / real (P(I))
          ENDDO
        ENDIF
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE ESTRIP (NEL,NUMS,INUMS,FACETS,PL2,IFCETS,NFCETS,NDIME
     &      ,RINGS,IRINGS)
C
C     This strips the facets into a set of edge-polygons
C
C     FACETS Data as:  Element #, Facet # (1-6), 'touching' facet #
C
      INTEGER
     +     NUMS (INUMS,*)        !- The element steerings
     +  ,FACETS (7,IFCETS)       !- 3->7 so can store edge info too?
     +   ,RINGS (IRINGS)         !- edge polygons (output)

     +   ,LINES (2,9999)         !- workspace for the edges :-(

     +     ,PL2 (IFCETS)         !- workspace only
     +     ,NUM (27)             !- steering for one element
     +     ,FS  (19)             !- nodes on a facet


      DO I=1,NFCETS
        PL2(I) = 0      ! mark all facets as 'not yet done'
      ENDDO
c---------------------- loop the 'materials' ---------------------------
      DO IF_BASE = 1,NFCETS
c.. hmm kill IMAT=0 at this stage ?
        IF (PL2(IF_BASE).NE.0) GOTO 100     !- 'cycle' to next

c------------------- loop the facets (of this material) ----------------
      NLINES = 0         !- no lines stored yet
      DO IFC = IF_BASE, NFCETS        !- loop the list

        IEL   = FACETS(1,IFC)
        IFACE = ABS(FACETS(2,IFC))
        ITOUCH= FACETS(3,IFC)
          CALL GET_ELEMENT (NUMS,INUMS,IEL, NUM, NOD,NDIME,ITYPE
     +     ,IMAT,IU1,IU2,IU3)
c       IF (IMAT.EQ.0) GOTO 81        !- skip this facet (.LE.0) ?
c                                     ! this should have happened earlier

        CALL GFACE (NOD,NDIME,ITYPE,IFACE,FS,NN_F,NN_FT, 2)  ! a facet's
        DO I=1,NN_F                                          ! polygon
          FS(I) = NUM(FS(I))  
        ENDDO
c.. got the nodes on this face .. now loop them as 'edges'
c.. better to abstract them as whole edges (eg. an 8nq has 4 edges)

      DO IBIT = 1,NN_F
        I1 = FS(I)
        I2 = FS(MOD(I,NN_F)+1)
        IFROM = MIN (I1,I2)           !- lo node
        ITO   = MAX (I1,I2)           !- hi node
        DO ILINES = 1,NLINES
          IF (IFROM.EQ.LINES(1,ILINES)) THEN
            IF (ITO.EQ.LINES(2,ILINES)) THEN        !- got a 'hit'
              LINES (1,ILINES) =  LINES (1,NLINES)  !-    so 
              LINES (2,ILINES) =  LINES (2,NLINES)  !- collapse out
              NLINES = NLINES - 1                   !-  this line
              GOTO 123   !- jump out to the next edge-piece
            ENDIF
          ENDIF
        ENDDO
        NLINES = NLINES + 1         !- no 'hit' so add 
        LINES (1,NLINES) = IFROM    !- to the end of the table
        LINES (1,NLINES) = ITO
  123   CONTINUE  !- 'jump out' after a 'hit'

      ENDDO     !-- loop the edge-bits of each facet

C--------------- Now daisy-chain the pieces together -------------------
C... so get the 'start' node. then build a chain
c.. hmm NLINES should drop by one each time we use a piece
      RINGS(1) = 0   !- No rings yet found
      IC = 2         !- current pointer in RINGS()

C----------------------- loop the polygons -----------------------------
      IL_BASE = 0
      DO WHILE (IL_BASE.LT.NLINES)
        IL_BASE = IL_BASE + 1
        ISTART = LINES(1,NLINES)   
        IFIND = ISTART
        IRING_START = IC            !- where to store the polygon length
c.. skip if already done?

c---------------------- build up each polygon --------------------------
        ILINE = IL_BASE -1
        DO WHILE (ILINE.LT.NLINES)
          ILINE = ILINE + 1
         IF (LINES(1,ILINE).EQ.IFIND ) THEN
           IEND = 1
         ELSEIF (LINES(2,ILINE).EQ.IFIND ) THEN
           IEND = 2
         ELSE
           IEND = 0
         ENDIF
c----------------------- found an edge segment -------------------------
c.. hmm do I re-store the first piont at the end ?
         IF (IEND.NE.0) THEN
            IC = IC + 1                            !- # of points 
            RINGS(iC) = LINES(IEND,ILINE)          !- store first
            IFIND = LINES(3-IEND,ILINE)            !- remember last
            LINES(1,ILINE) = LINES(1,NLINES)      ! / contract
            LINES(2,ILINE) = LINES(2,NLINES)      !|    the 
            NLINES = NLINES - 1                   ! \  list

            IF (IFIND.EQ.ISTART) THEN  !------- loop has been closed ---------
              IC = IC + 1
              RINGS(IC) = ISTART          !- store the first-point again
              RINGS(IRING_START) = IC     !- store polygon-length
              RINGS(1) = RINGS(1) + 1     !- #rings++
              PRINT*, RINGS(1),' : # Points = ', IC-IRING_START
c             ... jump to next starting point ?
              GOTO 124
            ENDIF       !- polygon now complete
          ENDIF       !- found an edge segment
        ENDDO       !- loop edges
        STOP 'Shouldn''t ever get to this line (incomplete polygon)'
 124    CONTINUE  !- completed polygon
      ENDDO     !- loop start of each polygon
C-----------------------------------------------------------------------
        ENDDO       !- loop facets
 100    CONTINUE       !- 'skip-to' if this facet already done
      ENDDO      !- loop 'materials'

      PRINT*, RINGS(1),' Polygons were found'
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------


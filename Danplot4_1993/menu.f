C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
      SUBROUTINE POST_MENUS (IOP,ixm,iym,char,ival)
C
C    This maintains the Menus
C    if IOP = 1 then all the 'current' menues are re-drawn
C    if IOP = 2 then the menus are searched for a given button-press
C
      parameter (n_menus=17)
      character char*1
      integer post(n_menus)
      data (post(i),i=1,n_menus)     !- initial menus
     + /2       !1= window-title-line
     + ,2       !2= pull-down menu-line
     + ,2       !3= bottom colour-line
     + ,1       !4= side logo pane

     + ,0       !5= - 'File' 
     + ,0       !6= - 'Col-mesh'
     + ,0       !7= - 'Palette'

     + ,0       !8= - 'Edit'
     + ,0       !9= - 'View'
     + ,0       !10=- 'Contour'
     + ,0       !11=- 'SVGA setup'

     + ,0       !12=  'Faces & edges' colours'
     + ,2       !13=  'screen NXE x NYE pos.' - Image itself
     + ,0       !14=  'image attributes' - background colour etc.
     + ,0       !15=  'General Program (menu) config options'

     + ,0       !16=  'Animation'
     + ,0       !17=  'Help'
     + /

c.........also 'printer setup', 'SVGA mode' , 'animation', etc.

c------------------------- draw menus ----------------------------------
      IF (IOP.EQ.1) THEN
        DO I=1,N_MENUS
          IF (POST(I).GE.1) CALL POST_A_MENU 
     +                 (I,1,IXM,IYM,CHAR,IVAL)
        ENDDO
c------------------------- search menus --------------------------------
      ELSEIF (IOP.EQ.2) THEN
        ival = -2                       !- flag as 'unfound'
        char = ' '            !- unfound ?
        DO I=N_MENUS,1,-1
          IF (POST(I).GE.1) CALL POST_A_MENU 
     +                 (I,2,IXM,IYM,CHAR,IVAL)
          IF (char.ne.' ') GOTO 2         !- found the option
        ENDDO
c... not found so exit
c       CALL SOUND@ (1,1000)
        RETURN

c------------------- handle some of the 'events' here ------------------
    2   IMENU = I
        IF (CHAR.EQ.'m') THEN                !- handle 'post-menu'
          IF (IVAL.GE.1.and.ival.le.n_menus) THEN   !- Post a new menu
            DO I=1,N_MENUS
              IF (POST(I).EQ.1) POST(I) = 0  !- KILL 'transients'
            ENDDO
            POST (IVAL) = 1                  !- = 'ON'
            CALL POST_A_MENU (IVAL,1,IXM,IYM,CHAR,IVAL)
          ELSEIF (IVAL.EQ.-1) THEN           !- repost all
            DO I=1,N_MENUS
              IF (POST(I).GE.1) 
     +        CALL POST_A_MENU (I,1,IXM,IYM,CHAR,IVAL)
            ENDDO
          ELSE
            stop 'unknown menu number' 
          ENDIF
          CHAR = ' '       !-  done opyion so 'switch off'
        ENDIF
      ELSE
        STOP '*** unknown option (POST_MENUS)' 
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE STANDARD_MENU_FONT(IOP)
C
C     Just a common entry-point for the 'standard' menu-font
C     .. can also change the standard font (and its size ?)
      SAVE
      COMMON /MENU_FONT/ IFONT
      IF (IOP.EQ.-1) THEN 
        IFONT = 101
        SIZE  = 1.
        IOFF  = 0.
      ELSEIF (IOP.EQ.0) THEN 
        CALL SET_TEXT_ATTRIBUTE@ ( IFONT,  SIZE, 0.,0.)
      ELSEIF (IOP.EQ.1) THEN
        PRINT*,'Which font do you want (1,100-112) ?'
        READ*,IFONT
        IF (IFONT.EQ.1) THEN
          IOFF = -8
        ELSE         ! set 'origin' for the text
          IOFF = 0
        ENDIF
      ELSEIF (IOP.EQ.2) THEN
        PRINT*,'Which size of font do you want (1.) ?'
        READ*,SIZE
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE POST_A_MENU (I,IOP,IXM,IYM,CHAR,IVAL)
C
C     This simply posts a menu by number (post or inquire!)
C
      character char*1
      call hide_mouse_cursor@()
      call standard_menu_font(0)
      if (i.eq.1) then
         call post_menu_1 (IOP,ixm,iym,char,ival)      
      elseif (i.eq.2) then
         call post_menu_2 (IOP,ixm,iym,char,ival)
      elseif (i.eq.3) then
         call post_menu_3 (IOP,ixm,iym,char,ival)
      elseif (i.eq.4) then
         call post_menu_4 (IOP,ixm,iym,char,ival)
      elseif (i.eq.5) then
         call post_menu_5 (IOP,ixm,iym,char,ival)
      elseif (i.eq.6) then
         call post_menu_6 (IOP,ixm,iym,char,ival)
      elseif (i.eq.7) then
         call post_menu_7 (IOP,ixm,iym,char,ival)
      elseif (i.eq.8) then
         call post_menu_8 (IOP,ixm,iym,char,ival)
      elseif (i.eq.9) then
         call post_menu_9 (IOP,ixm,iym,char,ival)
      elseif (i.eq.10) then
         call post_menu_10 (IOP,ixm,iym,char,ival)
      elseif (i.eq.11) then
         call post_menu_11 (IOP,ixm,iym,char,ival)
      elseif (i.eq.12) then
         call post_menu_12 (IOP,ixm,iym,char,ival)
      elseif (i.eq.13) then
         call post_menu_13 (IOP,ixm,iym,char,ival)
      elseif (i.eq.14) then
         call post_menu_14 (IOP,ixm,iym,char,ival)
      elseif (i.eq.15) then
         call post_menu_15 (IOP,ixm,iym,char,ival)
      elseif (i.eq.16) then
         call post_menu_16 (IOP,ixm,iym,char,ival)
      elseif (i.eq.17) then
         call post_menu_17 (IOP,ixm,iym,char,ival)
      endif
      call display_mouse_cursor@()
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE POST_MENU_1 (iop,ixm,iym,char,ival)
c
c     This is menu #1 ... the title-bar (plus 'exit' button etc)
c     if IOP = 1 then the menu is posted
c     if IOP = 2 then the op-codes for a 'hit' on one of this menus 
c     items is returned (maybe 'nul' if no 'hit')
c
      parameter (nbut=3)
      integer q(4),q1(4),p(4,nbut), opc(nbut)
      character  text(nbut)*90, opt(nbut)*1, CHAR*1,char1*1
      data ((p(j,i),j=1,4),text(i),opt(i),opc(i),i=1,nbut)  
     + / 
     +    20,  0,  600, 20,    'title', 'm' ,-1     !- all menus ? 
     +  ,  0,  0,   20, 20,    '-',     'q' ,0      !- exit
     +  ,620,  0,   20, 20,    'x',     'P' ,1      !- SVGA redraw
     + /      
      DATA Q/0,0,640,20/, Q1/0,0,0,0/

      IF (IOP.EQ.1) THEN
c     CALL POST_TEXT_WIDGET (Q1,Q, 1, 1, 1,0, 2,' ' )    !- menu pane

c.. hmm better to be able to set the title from the program !
      CALL SET_TEXT_ATTRIBUTE@ ( 107,  1.2, 0.,0.)       !- for title
      WRITE (TEXT(1),'(t11,a)') 
     +      '- Danplot v9 - ''\MODELS\KIDGER_6.PL'' -'
      i = 1
      CALL POST_TEXT_WIDGET   (Q,P(1,I), 5,1,1,0, 2,TEXT(I) ) 
      DO I = 2,3
        CALL POST_TEXT_WIDGET (Q,P(1,I), 2,1,1,0, 2,TEXT(I) )
      ENDDO
      ELSEIF (IOP.EQ.2) THEN
        CALL Q_MENU (IXM,IYM,Q,P,NBUT,ITEM,XM,YM)

        IF (ITEM.GT.0) THEN    !- if button found
          CHAR1 = OPT(ITEM)
          IVAL = OPC(ITEM)
            CHAR = CHAR1
        ENDIF

      ENDIF
      END

c-----------------------------------------------------------------------
      SUBROUTINE POST_MENU_2 (IOP,ixm,iym,char,ival)
c
c     This is menu #2 ...   the    TITLES across the top
c
c     if IOP = 1 then the menu is posted
c     if IOP = 2 then the op-codes for a 'hit' on one of this menus 
c     items is returned (maybe 'nul' if no 'hit')
c
      parameter (nbut=7)
      integer q(4),q1(4),p(4,nbut) ,opc(nbut)
      character  text(nbut)*12, opt(nbut)*1, char*1,CHAR1*1
      data ((p(j,i),j=1,4),text(i),opt(i),opc(i),i=1,nbut)  
     +  /
     +    10, 2, 45,16,    'File'         ,'m', 5
     +  , 60, 2, 45,16,    'Edit'         ,'m', 8
     +  ,110, 2, 65,16,    'Viewing'      ,'m', 9
c    +  ,180, 2, 65,16,    '       '      ,'m', 4     ! Mesh-generate ?
     +  ,180, 2, 65,16,    'Contour'      ,'m',10
     +  ,250, 2, 65,16,    'Mesh'         ,'m', 6
c    +  ,320, 2, 65,16,    'Faces+E'      ,'m',12   
     +  ,320, 2, 65,16,    'Animation'    ,'m',16
     +  ,400, 2, 65,16,    'Help'         ,'m',17
     +  /
      DATA Q/0,20,639,20/, Q1/0,0,0,0/

c------------------------ post menus -----------------------------------
      IF (IOP.EQ.1) THEN 
        CALL POST_TEXT_WIDGET (Q1,Q, 2, 1, 1,0, 2,' ' )    !- menu pane
      DO I=1,NBUT   
c       CALL POST_TEXT_WIDGET (Q,P(1,I), -1, 0, 0,2, 0,TEXT(I) )  
        CALL POST_TEXT_WIDGET (Q,P(1,I), -1, 1, 1,0, 0,TEXT(I) )
      ENDDO

c------------------------- search menus --------------------------------
      ELSEIF (IOP.EQ.2) THEN
        CALL Q_MENU (IXM,IYM,q,p,NBUT,ITEM,XM,YM)
        IF (ITEM.GT.0) THEN
          CHAR1 = OPT(ITEM)
          IVAL = OPC(ITEM)
          CHAR = CHAR1
        ENDIF
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE POST_MENU_3 (IOP,IXM,IYM,char,ival)
C
C     This is Menu # 3 . . . . . . . . . . the  NUMBERS
C
C     if IOP = 1 then the menu is posted
C     if IOP = 2 then the op-codes for a 'hit' on one of this menus 
C     items is returned (maybe 'nul' if no 'hit')
C
      COMMON /PALETTE/PAL
      INTEGER   PAL (3,0:255)      !- 'current' colour palette

      save
      parameter (nbut=26)
      integer q(4),q1(4),p(4,nbut) , ibc(nbut),   opc(nbut)
      character  text(nbut)*12, opt(nbut)*1, CHAR*1,CHAR1*1
      DATA Q/0,455,639,25/, Q1/0,0,0,0/

      data ((p(j,i),j=1,4),text(i),ibc(i),opt(i),opc(i),i=1,22)  
     +  /
     +    60, 4, 30,16,   'Off' ,1   ,'C', -1
     +  , 90, 4, 20,16,   ' 0'  ,0   ,'C',  0
     +  ,110, 4, 20,16,   ' 1'  ,1   ,'C',  1
     +  ,130, 4, 20,16,   ' 2'  ,2   ,'C',  2

     +  ,155, 4, 20,16,   ' 3'  ,3   ,'C',  3
     +  ,175, 4, 20,16,   ' 4'  ,4   ,'C',  4
     +  ,195, 4, 20,16,   ' 5'  ,5   ,'C',  5
                                            
     +  ,220, 4, 20,16,   ' 6'  ,6   ,'C',  6
     +  ,240, 4, 20,16,   ' 7'  ,7   ,'C',  7
     +  ,260, 4, 20,16,   ' 8'  ,8   ,'C',  8
     +  ,280, 4, 20,16,   ' 9'  ,9   ,'C',  9
     +  ,300, 4, 20,16,   '10'  ,10  ,'C', 10

     +  ,325, 4, 20,16,   '11'  ,11  ,'C', 11
     +  ,345, 4, 20,16,   '12'  ,12  ,'C', 12
     +  ,365, 4, 20,16,   '13'  ,13  ,'C', 13
     +  ,385, 4, 20,16,   '14'  ,14  ,'C', 14
     +  ,405, 4, 20,16,   '15'  ,15  ,'C', 15

     +  ,430, 4, 20,16,   '16'  ,2   ,'C', 16
     +  ,450, 4, 20,16,   '17'  ,2   ,'C', 17
     +  ,470, 4, 20,16,   '18'  ,2   ,'C', 18
     +  ,490, 4, 20,16,   '19'  ,2   ,'C', 19
     +  ,510, 4, 20,16,   '20'  ,2   ,'C', 20
     +  /

      data ((p(j,i),j=1,4),text(i),ibc(i),opt(i),opc(i),i=22+1,nbut)  
     +  /
     +   535, 4, 20,16,   'us'  ,0   ,'C', -10    !- typed value

     +  , 10, 4, 20,16,   'P'   ,2   ,'m',  7     !- post-palette
     +  ,580, 4, 20,16,   '<'   ,2   ,'C', -11    !-  val -1
     +  ,610, 4, 20,16,   '>'   ,2   ,'C', -12    !-  val +1
     +  /

c--------------------------- post menu ---------------------------------
      IF (IOP.EQ.1) THEN
         call post_text_widget (Q1,Q, 2, 1, 1,0, 2,' ' )  !- pane
         call draw_line@ (ints(q(1)), ints(q(2)),
     +     ints(q(1)+q(3)),ints(q(2)),  1)      !- edging in white

      call set_text_attribute@ ( 101,  .6, 0.,0.)   !- 'small' text

      DO I=1,23 
        colb = (pal(1,ibc(i)) +pal(2,ibc(i)) +pal(3,ibc(i))) /3./256.
c       print*,'colb = ',colb
        if (colb.gt.0.5) icolt = 0    ! = a 'contrasting' colour
        if (colb.le.0.5) icolt = 1
        CALL POST_TEXT_WIDGET (Q,P(1,I), ibc(i),icolt, 0,0, 1,TEXT(I))
      ENDDO

      CALL STANDARD_MENU_FONT(0)
      DO I=24,NBUT 
        colb = (pal(1,ibc(i)) +pal(2,ibc(i)) +pal(3,ibc(i))) /3./256.
        if (colb.gt.0.5) icolt = 0    ! = a 'contrasting' colour
        if (colb.le.0.5) icolt = 1
        CALL POST_TEXT_WIDGET (Q,P(1,I), ibc(i),icolt,  1,0, 1,TEXT(I) )
      ENDDO

c-------------------------- search menu --------------------------------
      ELSEIF (IOP.EQ.2) THEN
        CALL Q_MENU (IXM,IYM,q,p,NBUT,ITEM,XM,YM)

        IF (ITEM.GT.0) THEN
          CHAR = OPT (ITEM)
          IVAL = OPC (ITEM)
c          IF (CHAR1.EQ.'m') THEN                  !- menu post
c            CHAR = CHAR1     !- pass-back
c            ITEM  = IVAL
c          ELSE
c            CHAR = CHAR1     !- pass-back
c          ENDIF
        ENDIF

c       IF (ITEM.NE.-1) CALL SET_PAL (ITEM,0,0,0,13)

      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE POST_MENU_4 (iop,ixm,iym,char,ival)
c
c     This is Menu # 4 . . . . . . . . . . the  DANPLOT logo
c
c     if IOP = 1 then the menu is posted
c     if IOP = 2 then the op-codes for a 'hit' on one of this menus 
c     items is returned (maybe 'nul' if no 'hit')
c
      save
      parameter (nbut=1,  ncols = 19)
      integer q(4),q1(4),p(4,nbut), opc(nbut),icols(ncols)
      character  text(nbut)*12, opt(nbut)*1, char*1,CHAR1*1
      data ((p(j,i),j=1,4),text(i),opt(i),opc(i),i=1,nbut)  
     + / 
     +    10,  10,  80, 16,    '    Logo    ', ' ' , 0
     + /      
      data q/0,40,120,415/,q1/0,0,0,0/          !- menu pane itself
c     DATA ICOLS / 0,0,0, 3,4,5, 6,7,8, 9,10,11, 12,13,14, 15,  0,0,0/
c     DATA ICOLS / 0,3,0, 4,0,5, 0,6,0, 7,0,8  ,  0, 9, 0, 10, 0,11,0/
      DATA ICOLS / 0,1,1, 1,1,1, 1,1,0, 1,1,1  ,  1, 1, 1,  1, 0, 0,0/

      IF (IOP.EQ.1) THEN    !----- post the menu ----------
        call post_text_widget (Q1,Q, 2,1, 1,0, 2,' ' )     !- background

        CALL FILL_RECTANGLE@
     +    (Q(1)+10,Q(2)+ 33,Q(1)+Q(3)-10,Q(2)+ 33+1, 0S)
      DO I=1,NBUT 
        CALL POST_TEXT_WIDGET (q,p(1,i), 2, 1, 0,1, 0,TEXT(I) )
      ENDDO
c................ the logo bit .....................
      CALL SET_TEXT_ATTRIBUTE@ ( 107,  6., 90.,0.)   !- logo font !
      DO I=1,NCOLS
        ICOL = ICOLS(I)
        CALL DRAW_TEXT@ ('Danplot',q(1)+q(3)-45+I,q(2)+q(4)-30+I, icol)
      ENDDO
c...................................................
      ELSEIF (IOP.EQ.2) THEN

        CALL Q_MENU (IXM,IYM,q,p,NBUT,ITEM,XM,YM)
        IF (ITEM.GT.0) THEN
          CHAR1 = OPT(IVAL)
          IVAL = OPC(ITEM)
          CHAR = CHAR1
        ENDIF

      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE POST_MENU_5 (IOP,ixm,iym,char,ival)
c
c     This is menu #5 ...   the   FILE menu 
c
c     if IOP = 1 then the menu is posted
c     if IOP = 2 then the op-codes for a 'hit' on one of this menus 
c     items is returned (maybe 'nul' if no 'hit')
c
      parameter (nbut= 8)
      integer q(4),q1(4),p(4,nbut) ,opc(nbut)
      character  text(nbut)*12, opt(nbut)*1, char*1,CHAR1*1
      data ((p(j,i),j=1,4),text(i),opt(i),opc(i),i=1,nbut)  
     +  /
     +   10, 10, 80,16,    'New  '          ,'#',  0
     +  ,10, 30, 80,16,    'Open..'         ,'#',  1
     +  ,10, 50, 80,16,    'Import..'       ,'#', 20
     +  ,10, 70, 80,16,    'Export..'       ,'#', 30

c    +  ,10, 50, 80,16,    'Re-Open'        ,'#',  3

     +  ,10, 95, 80,16,    'Print'          ,'P',  2      !- ??
     +  ,10,160, 80,16,    'Configure..'    ,'m', 15
     +  ,10,200, 80,16,    'Status..'       ,'Q',  0
     +  ,10,225, 80,16,    'Exit'           ,'q',  0
     +  /
      data q/0,40,120,415/,q1/0,0,0,0/          !- menu pane itself

      IF (IOP.EQ.1) THEN     !--------- post menus ------------
        CALL POST_TEXT_WIDGET (q1,q, 2,1, 1,0, 2,' ' )
        CALL FILL_RECTANGLE@
     +    (Q(1)+10,Q(2)+ 92,Q(1)+Q(3)-10,Q(2)+ 92+1, 0s)
        CALL FILL_RECTANGLE@
     +    (Q(1)+10,Q(2)+138,Q(1)+Q(3)-10,Q(2)+138+1, 0s)
        CALL FILL_RECTANGLE@
     +    (Q(1)+10,Q(2)+222,Q(1)+Q(3)-10,Q(2)+222+1, 0s)

      DO I=1,NBUT   
        CALL POST_TEXT_WIDGET (Q,P(1,i),
     +    2,1,0,1, 0,TEXT(I) )
      ENDDO
      ELSEIF (IOP.EQ.2) THEN  !-------- search menus ----------
        CALL Q_MENU (IXM,IYM,q,p,NBUT,ITEM,XM,YM)
        IF (ITEM.GT.0) THEN
          CHAR = OPT(ITEM)
          IVAL = OPC(ITEM)
        ENDIF
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE POST_MENU_6 (IOP,ixm,iym,char,ival)
c
c     This is menu #6 ...   the   ITEM_COLOURS 
c
c     if IOP = 1 then the menu is posted
c     if IOP = 2 then the op-codes for a 'hit' on one of this menus 
c     items is returned (maybe 'nul' if no 'hit')
c
      parameter (nbut=14)
      integer q(4),q1(4),p(4,nbut), opc(nbut)        
      character  text(nbut)*16, opt(nbut)*1, char*1,CHAR1*1
      data ((p(j,i),j=1,4),text(i),opt(i),opc(I),i=1,nbut)  
     +  /
     +   10, 10, 80,16,    'Faces+ edges'  ,'m', 12    !- as a menu

     +  ,10, 55, 80,16,    'Sub-edges'     ,'I',  5
     +  ,10, 75, 80,16,    'Nodes'         ,'I',  9
     +  ,10, 95, 80,16,    'Node size'     ,'I', 10
     +  ,10,120, 80,16,    'Node #'        ,'I', 11
     +  ,10,140, 80,16,    'Element #'     ,'I',  8    ! was 13
     +  ,10,160, 80,16,    'Material #'    ,'I',  7

     +  ,10,185, 80,16,    'Orig mesh E'   ,'I',  16
     +  ,10,205, 80,16,    'Orig mesh F'   ,'I',  16
     +  ,10,225, 80,16,    'Disp vectors'  ,'I',  12

     +  ,10,250, 80,16,    'Load step # '  ,'I',  39
     +  ,10,270, 80,16,    'Picture #   '  ,'I',  40

     +  ,10,360, 80,16,    '# sub-facets'  ,'I',  19
     +  ,10,380, 80,16,    'Shrink....'    ,'I',  31

     +  /
      data q/0,40,120,415/,q1/0,0,0,0/          !- menu pane itself

      IF (IOP.EQ.1) THEN     !--------- post menus ------------
        CALL POST_TEXT_WIDGET (Q1,Q, 2,1, 1,0, 2,' ' )
        CALL FILL_RECTANGLE@
     +    (Q(1)+10,Q(2)+ 33,Q(1)+Q(3)-10,Q(2)+ 33+1, 0S)
        CALL FILL_RECTANGLE@
     +    (Q(1)+10,Q(2)+118,Q(1)+Q(3)-10,Q(2)+118+1, 0S)
c        CALL FILL_RECTANGLE@
c     +    (Q(1)+10,Q(2)+158,Q(1)+Q(3)-10,Q(2)+138+1, 0S)
        CALL FILL_RECTANGLE@
     +    (Q(1)+10,Q(2)+181,Q(1)+Q(3)-10,Q(2)+181+1, 0S)

      DO I=1,NBUT   
        CALL POST_TEXT_WIDGET (Q,P(1,I),
     +    2,1,0,1, 0,TEXT(I) )
      ENDDO

c--------------------------- search menus ------------------------------
      ELSEIF (IOP.EQ.2) THEN  
        CALL Q_MENU (IXM,IYM,Q,P,NBUT,ITEM,XM,YM)
        IF (ITEM.GT.0) THEN
          CHAR1 = OPT(ITEM)
          IVAL = OPC(ITEM)
          IF (CHAR1.EQ.'m') THEN
            CHAR = CHAR1
            ITEM = ival                        !- eg. Face menu
c          ELSEIF (CHAR1.EQ.'I') THEN
c            CALL SET_CV (IVAL,999,  2)        !- change a value
          ELSE
            CHAR = CHAR1
          ENDIF

        ENDIF
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE POST_MENU_7 (IOP,ixm,iym,char,ival)
c
c     This is menu #7 ...   the   PALETTES
c
      save
      PARAMETER (NBUT=17)
      INTEGER Q(4),Q1(4),P(4,NBUT) , OPC(NBUT)
      CHARACTER  TEXT(NBUT)*14, OPT(NBUT)*1, char*1,CHAR1*1
      DATA ((P(J,I),J=1,4),TEXT(I),OPT(I),OPC(I),I=1,NBUT)  
     +  /
     +   10, 10,90,16,    'Presets..'     ,' ',  0
     +  ,10, 30,90,16,    'Default'       ,'.', 11
     +  ,10, 50,90,16,    'Rainbow'       ,'.', 12
     +  ,10, 70,90,16,    'Hot-iron'      ,'.', 13
     +  ,10, 90,90,16,    'Spectrum'      ,'.', 14

     +  ,10,110,100,16,    'Red -> green'  ,'.', 21
     +  ,10,130,100,16,    'Blue-> yellow' ,'.', 22
     +  ,10,150,100,16,    'Cyan->b->red'  ,'.', 23
     +  ,10,170,100,16,    'Gray'          ,'.', 25
     +  ,10,190,100,16,    'User-defined'  ,'.', 24
     +  ,10,210,100,16,    'Shade..'       ,'I', 201   !- via 'items' ?

     +  ,10,300,100,14,    'Adjust..'      ,' ',  0   ! #11 ?
     +  ,10,315,100,14,    'Reverse'       ,'.',  2
     +  ,10,330,100,14,    'Invert'        ,'.',  3
     +  ,10,345,100,14,    'B <-> W'       ,'.',  4
     +  ,10,360,100,14,    'Swap evens'    ,'.',  5
     +  ,10,380,100,14,    'define..'      ,'I', 202
     +  /
      data q/0,40,120,415/,q1/0,0,0,0/          !- menu pane itself

c---------------------------- post menus -------------------------------
      IF (IOP.EQ.1) THEN

      CALL POST_TEXT_WIDGET (Q1,Q,  2,1, 1,0, 2,' ' )
      DO I=1,NBUT   
        CALL POST_TW_2 (Q,P(1,i),text(i),opt(i))
      ENDDO

c--------------------------- search menus ------------------------------
      ELSEIF (IOP.EQ.2) THEN
        CALL Q_MENU (IXM,IYM,Q,P,NBUT,ITEM,XM,YM)
        IF (ITEM.GT.0) THEN
          CHAR1 = OPT(ITEM)
          IVAL = OPC(ITEM)
          IF (CHAR1.EQ.'.') THEN              !- do palette set directly
            CALL SET_PAL (IVAL, 0,0,0,13)
          ELSE
            CHAR = CHAR1
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE POST_MENU_8 (IOP,ixm,iym,char,ival)
C
C     This is menu #8 ...   the   EDIT menu   ??
C
C     if IOP = 1 then the menu is posted
C     if IOP = 2 then the op-codes for a 'hit' on one of this menus 
C     items is returned (maybe 'nul' if no 'hit')
C
      save
      parameter (nbut=18)
      integer q(4),q1(4),p(4,nbut) , opc(nbut)
      character  text(nbut)*12, opt(nbut)*1, char*1,CHAR1*1
      data ((p(j,i),j=1,4),text(i),opt(i),opc(i),i=1,nbut)  
     +  /

     +   10, 40, 80,16,    'mirror:'       ,' '  ,0
     +  ,10, 60, 20,16,    'x'             ,'R' ,1
     +  ,40, 60, 20,16,    'y'             ,'R' ,2
     +  ,70, 60, 20,16,    'z'             ,'R' ,3
     +  ,10, 80, 20,16,    'X'             ,'R' ,4
     +  ,40, 80, 20,16,    'Y'             ,'R' ,5
     +  ,70, 80, 20,16,    'Z'             ,'R' ,6
                                              
     +  ,10,100, 80,16,    'un-mirror:'    ,' '  ,0
     +  ,10,120, 30,22,    '1/2'           ,'R' ,7
     +  ,45,120, 30,22,    '1/4'           ,'R' ,8
     +  ,80,120, 30,22,    '1/8'           ,'R' ,9

     +  ,10,160, 80,16,    'linkage'       ,'R'  ,50
     +  ,10,180, 80,16,    '"Harden"'      ,'R'  ,51


     +  ,10,320, 80,16,    'Load steps:'   ,' ' ,0
     +  ,10,340, 40,16,    'Sum'           ,'z' ,0
     +  ,60,340, 40,16,    'Diff'          ,'z' ,0

     +  ,10,360, 80,16,    'x->y->z'       ,'a' ,0
     +  ,10,380, 80,16,    'Glass..  '     ,'I' ,207  !- ???
     +  /
      data q/0,40,120,415/,q1/0,0,0,0/          !- menu pane itself

c------------------------ post menus -----------------------------------
      IF (IOP.EQ.1) THEN     

      CALL POST_TEXT_WIDGET (Q1,Q, 2,1,1,0, 2,' ' )  !- menu pane
      DO I=1,NBUT   
        CALL POST_TW_2 (Q,P(1,i),text(i),opt(i))

c        IF (OPT(I)(1:1).EQ.' ') THEN       !- do titles differently
c          CALL POST_TEXT_WIDGET (Q,P(1,i), 2,0, 1,0, 0,text(i) )
c        ELSE
c          CALL POST_TEXT_WIDGET (Q,P(1,i), 2,1, 1,0, 1,text(i) )
c        ENDIF
      ENDDO
c------------------------- search menus --------------------------------
      ELSEIF (IOP.EQ.2) THEN
        CALL Q_MENU (IXM,IYM,Q,P,NBUT,ITEM,XM,YM)
        IF (ITEM.GT.0) THEN
          CHAR1 = OPT(ITEM)
          IVAL = OPC(ITEM)
          IF (CHAR1.EQ.'.') THEN
            CALL SET_PAL (IVAL, 0,0,0,13)
          ELSE
            CHAR = CHAR1
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE POST_MENU_9 (IOP,ixm,iym,char,ival)
C
C     This is menu #9 ...   the   VIEW menu   ??
C
      save
      parameter (nbut=33)
      integer q(4),q1(4),p(4,nbut) , opc(nbut)
      character  text(nbut)*14, opt(nbut)*1, char*1,CHAR1*1
      data ((p(j,i),j=1,4),text(i),opt(i),opc(i),i=1,nbut)  
     +  /
     +   10, 10, 90,16,    '--- view ---'  ,' '  ,0

     +  ,10, 30, 90,16,    'image size'    ,' '  ,0
     +  ,10, 50, 20,16,    '<'             ,'-'  ,0
     +  ,35, 50, 40,16,    'typed'         ,'*'  ,0
     +  ,80, 50, 20,16,    '>'             ,'+'  ,0

     +  ,10, 70, 90,16,    'disp scale'    ,' '  ,0
     +  ,10, 90, 20,16,    '<'             ,'<'  ,0
     +  ,35, 90, 40,16,    'typed'         ,'s'  ,0
     +  ,80, 90, 20,16,    '>'             ,'>'  ,0

     +  ,10,110, 90,16,    'load step#..'  ,'I'  ,204      !- via color-menu
     +  ,10,130, 20,16,    '<'             ,'z'  ,0     !- ?check
     +  ,35,130, 40,16,    'typed'         ,'l'  ,0
     +  ,80,130, 20,16,    '>'             ,'z'  ,0     !- ?check
     + ,105,130, 20,16,    'L'             ,'F'  ,1   

     +  ,10,170, 90,16,    'backface-cull' ,' '  ,0
     +  ,10,190, 33,16,    'B'             ,'\'  ,1
     +  ,45,190, 33,16,    'n'             ,'\'  ,0
     +  ,80,190, 33,16,    'F'             ,'\'  ,2

     +  ,10,210, 90,16,    'Depth-sort'    ,' '  ,0
     +  ,80,210, 33,16,    'off'           ,'|'  ,0
     +  ,10,230, 33,16,    'B'             ,'|'  ,1
     +  ,45,230, 33,16,    'M'             ,'|'  ,3
     +  ,80,230, 33,16,    'F'             ,'|'  ,2

     +  ,10,260, 90,16,    'x'' <-> y'' '  ,'%'  ,0

     +  ,10,300, 90,16,    'colours..'     ,' '  , 0
     +  ,10,320, 70,16,    'background'    ,'I'  ,20
     +  ,83,320, 20,16,    'Os'            ,'I'  ,203   !- hack
     +  ,10,340, 70,16,    'borderline'    ,'I'  ,28
                                              
     +  ,10,360, 50,16,    'pic #'         ,'I'  ,41
     +  ,60,360, 25,16,    'Nx'            ,'I'  ,42
     +  ,90,360, 25,16,    'Ny'            ,'I'  ,43

     +  ,10,380, 50,16,    'save.'         ,'I'  ,301   !<- hack
     +  ,65,380, 50,16,    'load.'         ,'I'  ,302
     +  /
      data q/0,40,120,415/,q1/0,0,0,0/          !- menu pane itself

c------------------------ post menus -----------------------------------
      IF (IOP.EQ.1) THEN     
        CALL POST_TEXT_WIDGET (Q1,Q, 2,1,1,0, 2,' ' )  !- menu pane

      DO I=1,NBUT   
        IF (OPT(I)(1:1).EQ.' ') THEN       !- do titles differently
          CALL POST_TEXT_WIDGET (Q,P(1,i), 2,0, 1,0, 0,text(i) )
        ELSE
          CALL POST_TEXT_WIDGET (Q,P(1,i), 2,1, 1,0, 1,text(i) )
        ENDIF
      ENDDO
c... marbles 
        CALL SET_CV (32, IVAL, 5)    !- get current mode
        DO I=1,3
          IBUT = I + 14   !-- offset
          CALL POST_MARBLE (opc(ibut).eq.ival,Q,P(1,IBUT) ) 
        ENDDO
        CALL SET_CV (17, IVAL, 5)    !- get current mode
        DO I=1,4
          IBUT = I + 18   !-- offset
          CALL POST_MARBLE (opc(ibut).eq.ival,Q,P(1,IBUT) )
        ENDDO

c------------------------- search menus --------------------------------
      ELSEIF (IOP.EQ.2) THEN
        CALL Q_MENU (IXM,IYM,Q,P,NBUT,ITEM,XM,YM)
        IF (ITEM.GT.0) THEN
          CHAR1 = OPT(ITEM)
          IVAL = OPC(ITEM)
          IF (CHAR1.EQ.'\') THEN        !---------- b_p_cull -----------
           CALL SET_CV (32,ival  ,1)      !- set new value
            DO I=1,3
              IBUT = I + 14   !-- offset
              CALL POST_MARBLE (opc(ibut).eq.ival,Q,P(1,IBUT) )
            ENDDO
        ELSEIF (CHAR1.EQ.'|') THEN      !--------- depth sort ----------
           CALL SET_CV (17,ival  ,1)                     !- set new value
            DO I=1,4
              IBUT = I + 18   !-- offset
              CALL POST_MARBLE (opc(ibut).eq.ival,Q,P(1,IBUT) )
            ENDDO
          ELSE
            CHAR = CHAR1
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE POST_MENU_10 (IOP,ixm,iym,char,ival)
C
C     This is menu #9 ...   the   CONTOURING   menu
C
      SAVE
      PARAMETER (NBUT=21)
      INTEGER Q(4),Q1(4),P(4,NBUT) , OPC(NBUT)
      CHARACTER  TEXT(NBUT)*16, OPT(NBUT)*1, CHAR*1,char1*1
      DATA ((P(J,I),J=1,4),TEXT(I),OPT(I),OPC(I),I=1,NBUT)  
     +  /
     +   10, 10, 90,16,    '- Contours -'      ,' ',  0

     +  ,10, 40, 90,16,    'Disps:'            ,' ' ,  0 

     +  ,10, 60, 90,16,    'ë-total'           ,'K'  ,30

     +  ,10, 80, 30,16,    'ëx'                ,'K'  ,31
     +  ,45, 80, 30,16,    'ëy'                ,'K'  ,32
     +  ,80, 80, 30,16,    'ëz'                ,'K'  ,33
     +  ,10,100, 30,16,    '46'                ,'K'  ,46
     +  ,45,100, 30,16,    '47'                ,'K'  ,47
     +  ,80,100, 30,16,    '48'                ,'K'  ,48

     +  ,10,120, 90,16,    'Others...'         ,'z' ,  0 
     +  ,10,140, 90,16,    'Shear Strain'      ,'K' ,  47 
     +  ,10,160, 90,16,    'Vol. Strain'       ,'K' ,  46 
     +  ,10,180, 90,16,    'Vec. Normals'      ,'K' ,  18 
     +  ,10,200, 90,16,    '  + IMAT'          ,'K' ,  18 
             
     +  ,10,260, 60,16,    'Rescale: '         ,' ' ,  0 
     + , 80,260, 16,16,    'A'                 ,'^' ,  1 
     + ,100,260, 16,16,    'M'                 ,'^' ,  2 

     +  ,10,300, 90,16,    '# contours'        ,'I' ,18
     +  ,10,320, 90,16,    'line type..'       ,'I' ,14
     +  ,10,340, 90,16,    'face type..'       ,'I' ,15
     +  ,10,360, 90,16,    '..as a menu'       ,'z' ,  0

     +  /
      data q/0,40,120,415/,q1/0,0,0,0/          !- menu pane itself

c------------------------ post menus -----------------------------------
      IF (IOP.EQ.1) THEN     

      CALL POST_TEXT_WIDGET (Q1,Q, 2,1,1,0, 2,' ' )  !- menu pane
      DO I=1,NBUT   
        CALL POST_TW_2 (Q,P(1,i),text(i),opt(i))
c        CALL POST_TEXT_WIDGET (Q,P(1,i),
c     +    2,1,0,1, 1,TEXT(I) )
      ENDDO
c------------------------- search menus --------------------------------
      ELSEIF (IOP.EQ.2) THEN
        CALL Q_MENU (IXM,IYM,Q,P,NBUT,ITEM,XM,YM)
        IF (ITEM.GT.0) THEN
          CHAR1 = OPT(ITEM)
          IVAL = OPC(ITEM)
c          IF (CHAR1.EQ.'.') THEN
c            CALL SET_PAL (IVAL, 0,0,0,13)     !- n/a !
c          ELSE
            CHAR = CHAR1
c          ENDIF
        ENDIF
      ENDIF
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE POST_MENU_11 (IOP,ixm,iym,char,ival)
C
C     This is menu #11 ...   the   SVGA   modes menu
C
      save
      parameter (nbut=13)
      integer q(4),q1(4),p(4,nbut) ,opc(nbut)
      character      text(nbut)*15, opt(nbut)*1, char*1,char1*1

      integer*2 igx(40), igy(40), igm(40)
      integer*4 igc(40)
      logical*2 igb(4)

      logical built 
      data built /.false./
      data ((p(j,i),j=1,4),text(i),opt(i),opc(i),i=1,nbut)  
     +  /
     +   10, 10,100,16,    '-- SVGA modes --'   ,' ',    0
     +  ,10, 50,100, 0,    'n/a'                ,'V', 1
     +  ,10, 70,100, 0,    'n/a'                ,'V', 2
     +  ,10, 90,100, 0,    'n/a'                ,'V', 3
     +  ,10,110,100, 0,    'n/a'                ,'V', 4
     +  ,10,130,100, 0,    'n/a'                ,'V', 5
     +  ,10,150,100, 0,    'n/a'                ,'V', 6
     +  ,10,170,100, 0,    'n/a'                ,'V', 7
     +  ,10,190,100, 0,    'n/a'                ,'V', 8
     +  ,10,210,100, 0,    'n/a'                ,'V', 9
     +  ,10,230,100, 0,    'n/a'                ,'V',10
     +  ,10,250,100, 0,    'n/a'                ,'V',11 
     +  ,10,270,100, 0,    'n/a'                ,'V',12
     +  /
      data q/0,40,120,415/,q1/0,0,0,0/          !- menu pane itself

c------------------------ build menus ----------------------------------
      IF (.NOT.BUILT) THEN
        CALL GET_GRAPHICS_MODES@( IGX,IGY,IGC,IGM,IGB) 
        DO I=1,40
          IBUT = I+1
          IF ( IGX(I).EQ.-1) goto 2
          P (4,IBUT) = 18
          WRITE (TEXT(IBUT),'(i4,a,i3,a,i3)')
     +   IGX(I), ' ',IGY(I), ' ', IGC(I)
        ENDDO
    2   Nmodes = I
        BUILT = .TRUE.
      ENDIF

c------------------------ post menus -----------------------------------

      IF (IOP.EQ.1) THEN     
        CALL POST_TEXT_WIDGET (Q1,Q, 2,1,1,0, 2,' ' )      !- menu pane
c       CALL STANDARD_MENU_FONT (0)    
        CALL SET_TEXT_ATTRIBUTE@ ( 101,  .8, 0.,0.)   !(slightly smaller)
        DO I=1,NBUT   

        IF (OPT(I).EQ.' ') THEN          !- do titles differently
          CALL POST_TEXT_WIDGET (Q,P(1,i), 2,0, 1,0, 0,text(i) )
        ELSE
          CALL POST_TEXT_WIDGET (Q,P(1,i), 2,1, 1,0, 1,text(i) )
        ENDIF
        ENDDO
c............. radio buttons  'marbles' .........
        CALL SET_CV (50, IVAL, 5)    !- get current mode
        DO I=1,Nmodes
          IBUT = I + 1   !-- offset
          CALL POST_MARBLE (opc(ibut).eq.ival,Q,P(1,IBUT) )   !- On
        ENDDO
c------------------------- search menus --------------------------------
      ELSEIF (IOP.EQ.2) THEN
        CALL Q_MENU (IXM,IYM,Q,P,NBUT,ITEM,XM,YM)
        IF (ITEM.GT.0) THEN
          CHAR1 = OPT(ITEM)
          IVAL = OPC(ITEM)
          IF (CHAR1.EQ.'V') THEN
            CALL SET_CV (50,ival  ,1)               !- set new value
            DO I=1,nmodes      !-- draw marbles
              IBUT = I + 1     !-- offset
              CALL POST_MARBLE (opc(ibut).eq.ival,Q,P(1,IBUT) )
            ENDDO
          ELSE
            CHAR = CHAR1
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE POST_MENU_12 (IOP,ixm,iym,char,ival)
C
C     This is menu #12 ...   the   Faces & Edges menu
C
      save
      parameter (nbut=40)
      integer q(4),q1(4),p(4,nbut) ,opc(nbut)
      character      text(nbut)*15, opt(nbut)*1, char*1, CHAR1*1

      data ((p(j,i),j=1,4),text(i),opt(i),opc(i),i=1,nbut)  
     +  /
     +   10, 10,100,16,    'face & edges'       ,' ',    0

     +  ,10, 30, 70,16,    'Invis'              ,' ', 0
     +  ,10, 50, 70,16,    'Solid'              ,' ', 0
     +  ,10, 90, 70,16,    'Material'           ,' ', 0
     +  ,10,110, 70,16,    'Facet '             ,' ', 0
     +  ,10,130, 70,16,    'Front/back'         ,' ', 0
     +  ,10,150, 70,16,    'Glass '             ,' ', 0
     +  ,10,170, 70,16,    'Node #'             ,' ', 0
     +  ,10,190, 70,16,    'Element #'          ,' ', 0
     +  ,10,210, 70,16,    'nodal b/w'          ,' ', 0
     +  ,10,230, 70,16,    'chess?'             ,' ', 0
     +  ,10,250, 70,16,    'chess'              ,' ', 0 
     +  ,10,270, 70,16,    'shaded'             ,' ', 0


     +  ,87, 30, 16,16,    ' '                  ,'c', -1
     +  ,87, 50, 16,16,    ' '                  ,'c',  0
     +  ,87, 90, 16,16,    ' '                  ,'c',  1
     +  ,87,110, 16,16,    ' '                  ,'c',  2
     +  ,87,130, 16,16,    ' '                  ,'c',  3
     +  ,87,150, 16,16,    ' '                  ,'c',  4
     +  ,87,170, 16,16,    ' '                  ,'c',  5
     +  ,87,190, 16,16,    ' '                  ,'c',  6
     +  ,87,210, 16,16,    ' '                  ,'c',  7
     +  ,87,230, 16,16,    ' '                  ,'c',  8
     +  ,87,250, 16,16,    ' '                  ,'c',  9
     +  ,87,270, 16,16,    ' '                  ,'c', 10
         
     + ,100, 30, 16,16,    ' '                  ,'|', -1
     + ,100, 50, 16,16,    ' '                  ,'|',  0
     + ,100, 90, 16,16,    ' '                  ,'|',  1
     + ,100,110, 16,16,    ' '                  ,'|',  2
     + ,100,130, 16,16,    ' '                  ,'|',  3
     + ,100,150, 16,16,    ' '                  ,'|',  4
     + ,100,170, 16,16,    ' '                  ,'|',  5
     + ,100,190, 16,16,    ' '                  ,'|',  6
     + ,100,210, 16,16,    ' '                  ,'|',  7
     + ,100,230, 16,16,    ' '                  ,'|',  8
     + ,100,250, 16,16,    ' '                  ,'|',  9
     + ,100,270, 16,16,    ' '                  ,'|', 10

     +  ,20, 70, 60,16,    ' + color..'         ,' ',  0  
     + , 87, 70, 16,16,    'F'                  ,'I',  2  ! CV() -faces
     + ,100, 70, 16,16,    'E'                  ,'I',  4  ! CV() -edges

     +  /
      data q/0,40,120,415/,q1/0,0,0,0/          !- menu pane itself


c------------------------ post menus -----------------------------------

      IF (IOP.EQ.1) THEN     
        CALL POST_TEXT_WIDGET (Q1,Q, 2,1,1,0, 2,' ' )      !- menu pane
        DO I=1,NBUT   

        IF (OPT(I).EQ.' '.or.opt(i).eq.'c'.or.opt(i).eq.'|')  THEN 
          CALL POST_TEXT_WIDGET (Q,P(1,i), 2,0, 1,0, 0,text(i) )   !- plain
        ELSE
          CALL POST_TEXT_WIDGET (Q,P(1,i), 2,1, 1,0, 1,text(i) )   !- bdr
          IF (I.EQ.5)      !- *just a test*
     +    CALL POST_TEXT_WIDGET (Q,P(1,i), 2,1, 0,1, 1,text(i) )   !- bdr
        ENDIF
        ENDDO
c............. radio buttons  'marbles' .........
        CALL SET_CV (1, IVAL, 5)    !- get current mode
        DO I=1,12
          IBUT = I + 13  !-- offset
          CALL POST_MARBLE (opc(ibut).eq.ival,Q,P(1,IBUT) )   !- On
        ENDDO
        CALL SET_CV (2, IVAL, 5)    !- get current mode
        DO I=1,12
          IBUT = I + 25  !-- offset
          CALL POST_MARBLE (opc(ibut).eq.ival,Q,P(1,IBUT) )   !- On
        ENDDO
c------------------------- search menus --------------------------------
      ELSEIF (IOP.EQ.2) THEN
        CALL Q_MENU (IXM,IYM,Q,P,NBUT,ITEM,XM,YM)
        IF (ITEM.GT.0) THEN
          CHAR1= OPT (ITEM)
          IVAL = OPC (ITEM)
c          IF (CHAR1.EQ.'I') THEN
c           CALL SET_CV (IVAL,999,  2)        !- point to CV()
          IF (CHAR1.EQ.'c') THEN
              CALL SET_CV (1, IVAL, 1)         !- set new value
              DO I=1,12      !-- draw marbles
                IBUT = I + 13     !-- offset
                CALL POST_MARBLE (opc(ibut).eq.ival,Q,P(1,IBUT) )
              ENDDO
          ELSEIF (CHAR1.EQ.'|') THEN
              CALL SET_CV (3, IVAL, 1)         !- set new value

              DO I=1,12      !-- draw marbles
                IBUT = I + 25   !-- offset
                CALL POST_MARBLE (opc(ibut).eq.ival,Q,P(1,IBUT) )
              ENDDO
          ELSE
            CHAR = CHAR1
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE POST_MENU_13 (IOP,ixm,iym,char,ival)
c
c     This is menu #13 ...   the  PICTURE itself !
c
      SAVE
      PARAMETER (NBUT=1)
      INTEGER Q(4), P(4,NBUT)
      CHARACTER CHAR*1

      DATA P/  0, 0,519,415/
      DATA Q/120,40,519,415/          !- menu pane itself

c---------------------------- post menus -------------------------------
      IF (IOP.EQ.1) THEN
c........ nothing to draw !

c--------------------------- search menus ------------------------------
      ELSEIF (IOP.EQ.2) THEN
        CALL Q_MENU (IXM,IYM,Q,P,NBUT,ITEM,XM,YM)
        IF (ITEM.GT.0) THEN    !- A 'HIT' on this menu

        call sound@( 400s,1s)   !- debug
C.. get NXE,NYE hence pic #

          CALL SET_CV (42,NXE,5)       !- NXE
          CALL SET_CV (43,NYE,5)       !- NYE
          IPIC = NXE * INT (NYE * YM) + INT (NXE * XM) + 1
          CALL SET_CV (41,IPIC,1)       !- Picture number
          CHAR = '_'             !- redraw :-)
        ENDIF
      ENDIF
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE POST_MENU_14 (IOP,ixm,iym,char,ival)
c
c     This is menu #14 ...   the   Image attributes (background etc.)
c
c     if IOP = 1 then the menu is posted
c     if IOP = 2 then the op-codes for a 'hit' on one of this menus 
c     items is returned (maybe 'nul' if no 'hit')
c
      parameter (nbut=8)
      integer q(4),q1(4),p(4,nbut), opc(nbut)        
      character  text(nbut)*16, opt(nbut)*1, char*1,CHAR1*1
      data ((p(j,i),j=1,4),text(i),opt(i),opc(I),i=1,nbut)  
     +  /
     +   10, 55, 70,16,    'Background'    ,'I',  1
     +  ,10, 80, 70,16,    'Border'        ,'I',  1
     +  ,10,100, 80,16,    '# sub-facets'  ,'I',  19

     +  ,10,120, 80,16,    'Node size'     ,'I', 10
     +  ,10,140, 80,16,    'Element #'     ,'I', 13

     +  ,10,180, 80,16,    'Orig mesh E'   ,'I',  16
     +  ,10,200, 80,16,    'Orig mesh F'   ,'I',  16
     +  ,10,220, 80,16,    'Disp vectors'  ,'I',  12

     +  /
      data q/0,40,120,415/,q1/0,0,0,0/          !- menu pane itself

      IF (IOP.EQ.1) THEN     !--------- post menus ------------
        CALL POST_TEXT_WIDGET (Q1,Q, 2,1, 1,0, 2,' ' )
        CALL FILL_RECTANGLE@
     +    (Q(1)+10,Q(2)+ 33,Q(1)+Q(3)-10,Q(2)+ 33+1, 0S)
        CALL FILL_RECTANGLE@
     +    (Q(1)+10,Q(2)+ 78,Q(1)+Q(3)-10,Q(2)+ 78+1, 0S)
        CALL FILL_RECTANGLE@
     +    (Q(1)+10,Q(2)+138,Q(1)+Q(3)-10,Q(2)+138+1, 0S)
        CALL FILL_RECTANGLE@
     +    (Q(1)+10,Q(2)+161,Q(1)+Q(3)-10,Q(2)+161+1, 0S)

      DO I=1,NBUT   
        CALL POST_TEXT_WIDGET (Q,P(1,I),
     +    2,1,0,1, 0,TEXT(I) )
      ENDDO

c--------------------------- search menus ------------------------------
      ELSEIF (IOP.EQ.2) THEN  
        CALL Q_MENU (IXM,IYM,Q,P,NBUT,ITEM,XM,YM)
        IF (ITEM.GT.0) THEN
          CHAR1= OPT(ITEM)
          IVAL = OPC(ITEM)
c          IF (CHAR1.EQ.'I') THEN
c            CALL SET_CV (IVAL,999, 2)        !- change a value
C          ELSEIF (CHAR1.EQ.'m') THEN
C            ITEM = ival                        !- eg. Face menu
c          ELSE
            CHAR = CHAR1       !- pass back
c          ENDIF

        ENDIF
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE POST_MENU_15 (IOP,ixm,iym,char,ival)
c
c     This is menu #15 ...   the DANPLOT Configuration options
c
      parameter (nbut=7)
      integer q(4),q1(4),p(4,nbut), opc(nbut)        
      character  text(nbut)*16, opt(nbut)*1, char*1,CHAR1*1
      data ((p(j,i),j=1,4),text(i),opt(i),opc(I),i=1,nbut)  
     +  /
     +   10, 50, 80,16,    'Print. setup..' ,'z',  0
     +  ,10, 70, 80,16,    'SVGA mode..'    ,'m', 11

     +  ,10,100, 70,16,    'Menu color'    ,'.',  10
     +  ,10,120, 70,16,    'Menu font..'   ,'|',  1
     +  ,10,140, 70,16,    'Font size..'   ,'|',  2

     +  ,10,200, 80,16,    'Flash..'       ,'z', 0 
     +  ,10,220, 80,16,    'Beep..'        ,'z', 0 

     +  /
      data q/0,40,120,415/,q1/0,0,0,0/          !- menu pane itself

      IF (IOP.EQ.1) THEN     !--------- post menus ------------
        CALL POST_TEXT_WIDGET (Q1,Q, 2,1, 1,0, 2,' ' )

      DO I=1,NBUT   
        CALL POST_TEXT_WIDGET (Q,P(1,I), 2,1,0,1, 1,TEXT(I) )
      ENDDO

c--------------------------- search menus ------------------------------
      ELSEIF (IOP.EQ.2) THEN  
        CALL Q_MENU (IXM,IYM,Q,P,NBUT,ITEM,XM,YM)
        IF (ITEM.GT.0) THEN
          CHAR1= OPT(ITEM)
          IVAL = OPC(ITEM)
          IF (CHAR1.EQ.'.') THEN         !- do menu colour directly
            CALL SET_PAL (IVAL, 0,0,0,13)
          ELSEIF (CHAR1.EQ.'|') THEN     !- change menu font
            CALL STANDARD_MENU_FONT(IVAL)
          ELSE
            CHAR = CHAR1
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE POST_MENU_16 (IOP,ixm,iym,char,ival)
c
c     This is menu #16 ...   * Animation *
c
      parameter (nbut=2)
      integer q(4),q1(4),p(4,nbut), opc(nbut)        
      character  text(nbut)*16, opt(nbut)*1, char*1,CHAR1*1
      data ((p(j,i),j=1,4),text(i),opt(i),opc(I),i=1,nbut)  
     +  /
     +   10, 55, 70,16,    '# frames'      ,'z', 0

     +  ,10,160, 80,16,    'Animate..'     ,'z', 0 

     +  /
      data q/0,40,120,415/,q1/0,0,0,0/          !- menu pane itself

      IF (IOP.EQ.1) THEN     !--------- post menus ------------
        CALL POST_TEXT_WIDGET (Q1,Q, 2,1, 1,0, 2,' ' )

      DO I=1,NBUT   
        CALL POST_TEXT_WIDGET (Q,P(1,I), 2,1,0,1, 1,TEXT(I) )
      ENDDO

c--------------------------- search menus ------------------------------
      ELSEIF (IOP.EQ.2) THEN  
        CALL Q_MENU (IXM,IYM,Q,P,NBUT,ITEM,XM,YM)
        IF (ITEM.GT.0) THEN
          CHAR= OPT(ITEM)
          IVAL = OPC(ITEM)
        ENDIF
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE POST_MENU_17 (IOP,ixm,iym,char,ival)
c
c     This is menu #17 ...   * Help Menu *
c
      parameter (nbut=6)
      integer q(4),q1(4),p(4,nbut), opc(nbut)        
      character  text(nbut)*16, opt(nbut)*1, char*1,CHAR1*1
      data ((p(j,i),j=1,4),text(i),opt(i),opc(I),i=1,nbut)  
     +  /
     +   10, 20, 80,16,    'About..'        ,'m',  4

     +  ,10,100, 80,16,    'Sorry  '       ,'z', 0 
     +  ,10,120, 80,16,    '  No   '       ,'z', 0 
     +  ,10,140, 80,16,    ' Help  '       ,'z', 0 
     +  ,10,160, 80,16,    ' Yet   '       ,'z', 0 
     +  ,10,180, 80,16,    'Availaible'    ,'z', 0 

     +  /
      data q/0,40,120,415/,q1/0,0,0,0/          !- menu pane itself

      IF (IOP.EQ.1) THEN     !--------- post menus ------------
        CALL POST_TEXT_WIDGET (Q1,Q, 2,1, 1,0, 2,' ' )

      DO I=1,NBUT   
        CALL POST_TEXT_WIDGET (Q,P(1,I), 2,1,0,1, 1,TEXT(I) )
      ENDDO

c--------------------------- search menus ------------------------------
      ELSEIF (IOP.EQ.2) THEN  
        CALL Q_MENU (IXM,IYM,Q,P,NBUT,ITEM,XM,YM)
        IF (ITEM.GT.0) THEN
          CHAR = OPT (ITEM)
          IVAL = OPC (ITEM)
        ENDIF
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
c.....             Menu posting/inquiring routines
c-----------------------------------------------------------------------

      SUBROUTINE POST_TW_2 (q,p,TEXT,txt2)
C
C     A header to 'POST_TEXT_WIDGET' which colour the text of the button 
C     if a title (txt2=' ') :
C           text in *black* with no button box
C     else
C           text in *white* within a button box, with b+w shadow of 1 pixel
C
C   .. I could also have options to: disable colour, hi-light sub-menus,etc.
C
      INTEGER Q(4),P(4)
      CHARACTER TEXT*(*),TXT2*1

      IF (TXT2.EQ.' ') THEN                                !* A title *!
        CALL POST_TEXT_WIDGET (Q,P, 2,0, 1,0, 0,text)
      ELSE                                                 !* An option *!
        CALL POST_TEXT_WIDGET (Q,P, 2,1, 1,0, 1,text)
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE POST_TEXT_WIDGET (q,p,
     +                   ICOLB,ICOLT,ICOL1,ICOL2,ITHK,TEXT)
c                                                        
c     this posts a menu item 'text' at position 'pos' with the given 
c     colours
c   
c                                              Dan Kidger  29-11-91
c--> need to extend to include 'control' characters eg. colours, font,
c    'dials', 'sliders', etc. !
c.. also nice to optionaly switch off the shadowing of the buttons

      COMMON /MENU_FONT/ IFONT

      INTEGER Q(4),P(4)
      CHARACTER TEXT*(*)

      IXF = q(1) + p(1)
      IYF = q(2) + p(2)
      IXT = IXF  + p(3)
      IYT = IYF  + p(4)

      IF (P(4).EQ.0) RETURN   !- dummy skip of invisible buttons

      DO J=0,ithk -1                   !- 'thick' borders
        CALL DRAW_LINE@ (ints(IXF+J),ints(IYT-J),
     +                   ints(IXF+J),ints(IYF+J),ints(ICOL1))  ! ÉÍ   top
        CALL DRAW_LINE@ (ints(IXF+J),ints(IYF+J),
     +                   ints(IXT-J),ints(IYF+J),ints(ICOL1))  !      right
        CALL DRAW_LINE@ (ints(IXT-J),ints(IYF+J),
     +                   ints(IXT-J),ints(IYT-J),ints(ICOL2))  !   Í¼ bottom
        CALL DRAW_LINE@ (ints(IXT-J),ints(IYT-J),
     +                   ints(IXF+J),ints(IYT-J),ints(ICOL2))  !      left
      ENDDO
      IF (ICOLB.GE.0) CALL CLEAR_SCREEN_AREA@    !- cf. fill_rectangle@ ??
     +    (ints(IXF+ithk),ints(IYF+ithk)
     +   , ints(IXT-ithk),ints(IYT-ithk),ints(ICOLB))

      IBLACK = 0                            !- ('black' is col #0)
      IF (LENG(TEXT).GT.0) THEN              !- no need to do for 'blanks'

        IF (IFONT.EQ.1) THEN
          IOFF = -10
        ELSE         ! set 'origin' for the text
          IOFF = 0
        ENDIF

      IBDR = (IYT - IYF - 14 ) /2  + ithk-1 !- thickness of the 'border' ??

      IF (ICOLT.NE.IBLACK.and.ifont.eq.1)      !- shadow non-black letters
     + CALL DRAW_TEXT@ (TEXT(1:LENG(TEXT)),
     +  ints(IXF+ithk+1), ints(IYT-ibdr-4+IOFF), ints(IBLACK))

       CALL DRAW_TEXT@ (TEXT(1:LENG(TEXT)),
     +  ints(IXF+ithk+2), ints(IYT-ibdr-3+IOFF)  , ints(ICOLT))

      ENDIF
      END
c-----------------------------------------------------------------------
      SUBROUTINE POST_MARBLE (YES, Q, P )
c                                                        
c     This adds a 'marble' at the *end* of a (pre-posted) button
c     to indicate its status
c                                              Dan Kidger  9-6-93

      INTEGER Q(4),P(4)
      LOGICAL YES

      IF (P(4).EQ.0) RETURN        !- skip invisible buttons
      IXF = q(1) + p(1)
      IYF = q(2) + p(2)
      IXT = IXF  + p(3)
      IYT = IYF  + p(4)

      IXC = IXT - 10
      IYC = (IYF+IYT)/2
  
      IF (.NOT.YES) THEN                                    !- Off
        CALL FILL_ELLIPSE@ (ints(IXC), ints(IYC)  , 4S,4S, 2) !-gray
        CALL FILL_ELLIPSE@ (ints(IXC), ints(IYC)  , 1S,1S, 0)
      ELSEIF (YES) THEN                                     !- On
        CALL FILL_ELLIPSE@ (ints(IXC), ints(IYC)  , 4S,4S, 1)
        CALL FILL_ELLIPSE@ (ints(IXC), ints(IYC)  , 2S,2S, 0)
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
      SUBROUTINE Q_MENU (IXM,IYM, Q,P,N_MENU, ITEM,XM,YM)
c
c     this returns the menu item (and sub-postion) of a selection
c     ih,iv = mouse position (screen)
c     xm,ym = fractional position within a menu item
c
c... the table will be just one menu so no need for 'menup()'

      integer q(4), p(4,N_MENU)

      integer*4 buffer

      item=  -1
      xm  = .5
      ym  = .5
      IH = IXM - Q(1)
      IV = IYM - Q(2)
      IF (IXM.LT.Q(1)) RETURN     !- exit *FAST* if the point 
      IF (IYM.LT.Q(2)) RETURN     !- is outside the menu
      IF ( IH.GT.Q(3)) RETURN
      IF ( IV.GT.Q(4)) RETURN
      do i = n_menu,1,-1      !- (try reverse search)
        if (iv.LT.p(2,i))        goto 1
        if (iv.GT.p(2,i)+p(4,i)) goto 1   !- skip to next button
        if (ih.LT.p(1,i))        goto 1
        if (ih.GT.p(1,i)+p(3,i)) goto 1
        IF ( P(4,i).eq.0)        goto 1   !- 'invisible' button
                item = i
                xm = real(ih-p(1,i)) / real(p(3,i))
                ym = real(iv-p(2,i)) / real(p(4,i))
       
      IF (REAL(P(4,I)-P(2,I)) * REAL(P(3,I) - P(1,I)).LT.100.*100.) then
C... only for screen areas *less* than   100 by 100 pixels :-)
c... 'flash'

          call hide_mouse_cursor@()
          call get_screen_block@(q(1)+p(1,i),q(2)+p(2,i)
     +       ,q(1)+p(1,i)+p(3,i), q(2)+p(2,i)+p(4,1), buffer)
          call clear_screen_area@(q(1)+p(1,i),q(2)+p(2,i)
     +       ,q(1)+p(1,i)+p(3,i), q(2)+p(2,i)+p(4,1), 15S)
          call sleep@ (0.1)
          call restore_screen_block@ (q(1)+p(1,i),q(2)+p(2,i)
     +       ,buffer,3,ifail)                         !- XOR
c         call sleep@ (0.2)
          call restore_screen_block@ (q(1)+p(1,i),q(2)+p(2,i)
     +       ,buffer,0,ifail)                         !- replace
          call display_mouse_cursor@()
c... end-'flash'
      ENDIF
          return
    1  continue
      enddo
      return
      end
c-----------------------------------------------------------------------


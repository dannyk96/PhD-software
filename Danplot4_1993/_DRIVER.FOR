c     options (fullcheck,undef)

c-----------------------------------------------------------------------
      PROGRAM MENU
C-----------------------------------------------------------------------
C     This program is a 'driver' to my 'new' menuing system
C     .. It merely calls POST_MENUS to draw all the current menus
C     .. note that this is also Object-Oriented as it also calls this 
C     .. code to inquire which button has been pressed
C-----------------------------------------------------------------------
C    2-7-93 MENU1A .. works !, try posting menus direct
C
C-----------------------------------------------------------------------
C--      real a
c--      integer i
      integer*2 ixm_s,iym_s,i_but_s, ifail
      integer*4 buffer
      character*17 pal
      character char*1
c-----------------------------------------------------------------------
      call into_graphics ()
      call set_cv (999,999,0)              !- reset CV
      call post_menus (1, 0,0,char,ival)

c................ dummy picture .................
      call hide_mouse_cursor@()
      call pcx_to_screen_block@('teapot.pcx', buffer,pal,ifail)
      call restore_screen_block@(121s,41s,buffer,0s,ifail)

c.............. dummy picture border ...............
      call draw_line@ (121s, 41s, 639s, 41s,  14)
      call draw_line@ (639s, 41s, 639s,454s,  14)
      call draw_line@ (639s,454s, 121s,454s,  14)
      call draw_line@ (121s,454s, 121s, 41s,  14)
      call display_mouse_cursor@()

c------------------ Loop and handle mouse-events -----------------------
    2 call get_mouse_button_press_count@ (0S, i_but_s) 
c     call debounce_buttons ()
c     call yield@()
      if (i_but_s.eq.0) goto 2

c-------- OK a mouse press
c.. loop back to here until the mouse is released

  94  call get_mouse_pos@ (ixm_s,iym_s)
      ixm = ixm_s
      iym = iym_s
      char = 'z'                              !- defualt = 'no-op'
      call post_menus(2,ixm,iym,char,ival)      !- get the opcode ?



c------------------------------------------------------
      call get_mouse_position@( ixm_s,iym_s,ibut_s)
      if (ibut_s.ne.0) goto 94     !- until mouse released

c-------------------------- handle the events --------------------------

      if (1.eq.1) goto 2    !-- back to the event handler
c------------------------- Exit ----------------------------------------
      read*
      call text_mode@()

      end



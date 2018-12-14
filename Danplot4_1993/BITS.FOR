c-----------------------------------------------------------------------
      SUBROUTINE JUST_JUNK
C .. just an empty subroutine
      END
C-----------------------------------------------------------------------
      SUBROUTINE INVARIENTS (STRESS,NODOF,SIGM,DSBAR,THETA)
C
C     This forms the stress invariants from the stress components 2D/3D
C
      REAL STRESS(*)
      SQ3=SQRT(3.)
      SX=STRESS(1)
      SY=STRESS(2)
      IF (NODOF.EQ.2) THEN
        SZ=STRESS(4)      !- sigma-z
        S4=STRESS(3)      !- shear term
        S5 = 0.
        S6 = 0.
      ELSEIF (NODOF.EQ.3) THEN
        SZ=STRESS(3)
        S4=STRESS(4)
        S5=STRESS(5)
        S6=STRESS(6)
      ENDIF
      SIGM = (SX+SY+SZ)/3.
      D2   = ((SX-SY)**2+(SY-SZ)**2+(SZ-SX)**2)/6. +S4**2 +S5**2 +S6**2
      DSX = SX-SIGM
      DSY = SY-SIGM
      DSZ = SZ-SIGM
      D3 = DSX*DSY*DSZ - DSX*S5**2 - DSY*S6**2 - DSZ*S4**2 + 2.*S4*S5*S6
      DSBAR = SQ3*SQRT(D2)
      IF (ABS(DSBAR).LT.1.E-10) THEN      !- ie no shear
        THETA=0.
      ELSE
        SINE=-3.*SQ3*D3/(2.*SQRT(D2)**3)
        SINE = MIN(MAX(-1.,SINE),1.)       !-- clip the lode angle
        THETA=ASIN(SINE)/3.          !- (in radians)
      ENDIF
      RETURN
      END

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE MAKE_BOX (IX,IY,IXW,IYW,X,Y)
C
C     simply stores xlo/xhi info as a polygon (so draw_a_box)
C
      REAL X(4),Y(4)
      X(1) = IX
      X(2) = IX + IXW  -1
      X(3) = IX + IXW  -1
      X(4) = IX

      Y(1) = IY
      Y(2) = IY
      Y(3) = IY + IYW  -1
      Y(4) = IY + IYW  -1
      RETURN
      END
C-------------------------- Show status --------------------------------
      SUBROUTINE SHOW_STATUS (TIME,NN,NEL,NLDS,NFCETS_D,NFCETS,FILE_DAT)
C
C     This shows in the main drawing window 
C     a table of status information about the mesh (+timings)  
C
C...  show other anotation too !
C..... eg. load step # , viewing angle, light angle, disp scale, etc.
C---> 29-11-93 all the INTEGER*2 removed (for FTN90) cos' it still works

      REAL*4 TIME(20), T
      CHARACTER LINE*80,FILE_DAT*(80)
      T = 1.
      CALL SET_TEXT_ATTRIBUTE@ (1,T,T,T)
c     CALL POST_WIDGET (TABLE,TXT(menup(9)),menup(9))   !- obs

      write(line,'(a)')'   < '//
     +      FILE_DAT(1:leng(file_dat)-1)//' >'
        call draw_text@ (line(1:leng(line)),300, 70,14)
      write(line,'(a,i5)')  '   NN =',NN
        call draw_text@ (line(1:leng(line)),300,115,14)
      write(line,'(a,i5)')  '   NEL=',NEL
        call draw_text@ (line(1:leng(line)),300,130,14)
      write(line,'(a,i5)')  '#loads=',NLDS
        call draw_text@ (line(1:leng(line)),300,145,14)
      write(line,'(2(a,i4))')  '#facets=',NFCETS_D,'/',NFCETS
        call draw_text@ (line(1:leng(line)),300,160,1)

      write(line,'(a,f10.3)')  '--- timings ---'
        call draw_text@(line(1:leng(line)),300,300,1)
      write(line,'(a,f10.3)')  'total =',time(5)-time(1)
        call draw_text@(line(1:leng(line)),300,315,15)
      write(line,'(a,f10.3)')  'transf=',time(2)-time(1)
        call draw_text@(line(1:leng(line)),300,330,1)
      write(line,'(a,f10.3)')  'd.sort=',time(3)-time(2)
        call draw_text@(line(1:leng(line)),300,345,1)
      write(line,'(a,f10.3)')  'draw  =',time(4)-time(3)
        call draw_text@(line(1:leng(line)),300,360,1)
      write(line,'(a,f10.3)')  'buffer=',time(5)-time(4)
        call draw_text@(line(1:leng(line)),300,375,1)

        write(line,'(a,f10.3)')  'reading=',time(11)-time(10)
          call draw_text@(line(1:leng(line)),300,390,15)
        write(line,'(a,f10.3)')  'FSTRIP =',time(12)-time(11)
          call draw_text@(line(1:leng(line)),300,405,15)

      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GET_KEY_CHAR (KEY,CKEY)
C
C     This simple translates a simple 'keycode' into an ascii character
C
      CHARACTER CKEY*1
      IF (KEY.EQ.0)  THEN    !- 'no-op'
        CKEY = 'z'
      ELSEIF (KEY.EQ.32) THEN    !- alias for 'redraw' option
        CKEY = '_'                  
      ELSEIF (KEY.EQ.13) THEN       !- alias for 'redraw' option
        CKEY = '_'                  

      ELSEIF (KEY.EQ.331) THEN      ! - cusor   : left
         CKEY = CHAR(128)           
      ELSEIF (KEY.EQ.333) THEN      !           : right
         CKEY = CHAR(129)           
      ELSEIF (KEY.EQ.328) THEN      !           : down
         CKEY = CHAR(130)
      ELSEIF (KEY.EQ.336) THEN      !           : up
         CKEY = CHAR(131)

      ELSEIF (KEY.GE.0.and.KEY.LE.255) THEN
        CKEY = CHAR(KEY)
      ELSE
        CKEY = 'Z'        !- patch to show key is >255 (eg F9) ? obs?
      ENDIF
      RETURN
      END

c-----------------------------------------------------------------------
      SUBROUTINE DANPLOT_LOGO (RES)
C
C     Ok need more work on this to get a *nice* logo
C
      INTEGER RES(*)
      REAL X(20),Y(20)  ,X1(2),Y1(2)
      CHARACTER LINE*80

      IXS = RES(4)
      IYS = RES(5)
      IXW = RES(1)
      IYW = RES(2)

      CALL MAKE_BOX (IXS,IYS,IXW,IYW,X,Y)
      CALL DR_PRIM (LINE,X,Y,4,  9, 2)     !- the 'extra space'
      IXS = IXS + 30
      IYS = IYS + 30
      IXW = IXW - 30 *2
      IYW = IYW - 30 *2
      CALL MAKE_BOX (IXS,IYS,IXW,IYW,X,Y)
      CALL DR_PRIM (LINE,X,Y,4, 15, 2)     !----- 'box edge' ------
      IXS = IXS + 3                     !>  nice to be able to do  <!
      IYS = IYS + 3                     !>  rounded corners  :-)   <!
      IXW = IXW - 3 *2
      IYW = IYW - 3 *2
      CALL MAKE_BOX (IXS,IYS,IXW,IYW,X,Y)
      CALL DR_PRIM (LINE,X,Y,4, 3, 2)     !- 'main box'

      X1(1) = IXS  + 75
      Y1(1) = IYS  + 50
      WRITE(LINE,'(A)') 'DANPLOT-2'
      CALL SET_TEXT_ATTRIBUTE@ (107, 4.,0.,0.)
      CALL DR_PRIM (LINE,X1,Y1,1, 1, 20)     !- main title
                                      
                                   
      CALL SET_TEXT_ATTRIBUTE@ (107, 2. ,0.,0.)
      X1(1) = IXS +145
      Y1(1) = IYS +75
      WRITE (LINE,'(A)') 'version 2.0'
      CALL DR_PRIM (LINE,X1,Y1,1, 15, 20)     !- version number

      CALL SET_TEXT_ATTRIBUTE@ (107, 1.5,0.,0.)

      X1(1) = IXS +140
      Y1(1) = IYS +290
      WRITE (LINE,'(A)') 'by Dr. Dan Kidger'
      CALL DR_PRIM (LINE,x1,y1,1, 1, 20) 

      X1(1) = IXS + 70
      Y1(1) = IYS +310
      WRITE (LINE,'(A)') ' Dept. of Civil Engineering'
      CALL DR_PRIM (LINE,x1,y1,1, 15, 20) 

      X1(1) = IXS + 70
      Y1(1) = IYS +330
      WRITE (LINE,'(A)') 'University of Manchester, UK'
      CALL DR_PRIM (LINE,x1,y1,1, 15, 20) 


      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE WRITE_HELP ()
C
C     just writes out the 'usage' information for DANPLOT v2
C
      WRITE(*,'(a)')
     + 'DANPLOT finite element visualization package. Version 2.0'
     +,'Copyright 1993  Dr. Dan Kidger, all rights reserved.'
     +,' '
     +,'Usage: DANPLOT [keyfile] [-d file ] [-p#] [-PS] [-?]'
     +,' '
     +,'[keyfile]  Specifies keyword based data file'
     +,' '
     +,'-d file    Imports a mesh file where file ends in :'
     +,'               .PL    for old versions of DANPLOT' 
     +,'               .PL2   for old version containing "extra lines"'
     +,'               .GEO   for "Object File Format" '
     +,'               .NFF   for "Neutral File Format" '
     +,' '
     +,'-PS        Postscript output'
     +,'-p#        The palette number to use'
     +,'-%         Swaps X and Z on the screen '
     +,'-A         makes Z the vertical axis'
     +,' '
     +,'To get a full copy of the DANPLOT package please contact: '
      WRITE(*,'(a,t40,a)')
     + 'Dr. Dan Kidger'
     +                    ,' Tel   : 061-275-4402 / 4375'
     +,'Dept. of Engineering,'
     +                    ,' FAX   : 061-274-3384'  
     +,'University of Manchester, UK'
     +                    ,' E-MAIL: Dan.Kidger@manchester.ac.uk'
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE R_OPTIONS (IO,IOPS)
C
C     This reads in various program control options as an integer array
C     .. such as Output options for DISPS,ACTIONS,STRESSES etc.
C     .. and control such as Max.# of iterations, 'No-tension' supports
C         flag >limit / echo 'old'/'new' ??
C
      INTEGER IOPS(*)
      DO I=1,9999
    1   READ(IO,*,IOSTAT=IOS) IOP,IVAL
        CALL IN_TEST(IO,IOS,*1,*999)
        IOPS(IOP) = IVAL
      ENDDO
      CALL MYERROR (2,'-- should never get to this line') !-code = 3 ?
  999 RETURN
      END

C----------------------------------------------------------------------
      SUBROUTINE FORM_DEE (DEE,IDEE,E,V,BW,NODOF,IOP)
C
C     This forms most stress/strain matrices in 2D and 3D
C        Dan Kidger   c. 1989
C
C     Set IOP = (1)NORM, (2)K, (3)G, 4(LAME), 5(G),6(PLANE-STRESS)
C     BW = modulus of the water ( =0.0 if none)
C    
      REAL DEE(IDEE,*), FACT(10)
      INTEGER M(3,6)
      DATA M/5,4,2, 3,3,1, 7,8,2, 4,4,1, 6,1,2,  9,10,2/

      DO I=1,IDEE
        DO J=1,IDEE
          DEE(I,J) = 0.
        ENDDO
      ENDDO
C     CALL NULL (DEE,IDEE,IDEE,IDEE)

C----- fact= (1)0.,(2)G,(3)K,(4)LAME,(5)L2,(6)2G,(7)4G/3,(8)=2G/3 ------

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

      DO I=1,NODOF                     !----- direct terms -----
        DEE(I,I) = FACT(M(1,IOP))
        DO J=1,NODOF
          IF (I.NE.J) DEE(I,J) = FACT(M(2,IOP))
        ENDDO
      ENDDO
      DO I=NODOF+1, NODOF*(NODOF+1)/2  !----- shear terms ------
        DEE(I,I) = FACT(M(3,IOP))
      ENDDO

      IF (NODOF.EQ.2.AND.IDEE.GE.4) THEN  !---- Sz terms if 2D -----
        DEE (1,4) = DEE (1,2)
        DEE (2,4) = DEE (1,2)
        DEE (3,4) = 0.
        DEE (4,4) = DEE (1,1)

        DO I=1,4
          DEE(4,I) = DEE(I,4)   !- symmetry
        ENDDO
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE MVMULT (M,IM,V,K,L,Y)
C
C     This multiplies a matrix by a vector
C
      REAL M(IM,*),V(*),Y(*)
      DO I=1,K
        X = 0.
        DO J=1,L
          X = X + M(I,J)*V(J)      !-- reverse order ??
        ENDDO
        Y(I) = X
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE MTVMULT (M,IM,V,K,L,Y)
C
C     This multiplies a matrix-transposed by a vector
C
      REAL M(IM,*),V(*),Y(*)
      DO I=1,K
        X = 0.
        DO J=1,L
          X = X + M(J,I)*V(J)      !-- reverse order ??
        ENDDO
        Y(I) = X
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE MAT_MULT (A,IA,B,IB,C,IC,L,M,N)
C
C     This forms the product of two matrices
C
      REAL A(IA,*),B(IB,*),C(IC,*)
      DO I=1,L
        DO J=1,N
          X=0.0
          DO K=1,M
            X=X+A(I,K)*B(K,J)
          ENDDO
          C(I,J)=X
        ENDDO
      ENDDO
      END
c-----------------------------------------------------------------------
      SUBROUTINE TWOBY2 (JAC,IJAC,JAC1,IJAC1,DET)
C
C     This forms the inverse of a 2 by 2 matrix
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

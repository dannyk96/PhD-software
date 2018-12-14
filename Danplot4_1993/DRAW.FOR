C-----------------------------------------------------------------------
c >>>  This the GRAPHICS part of the code, incluing
c       1/ all œD transformations
c       2/ all Palette adjusting
c       3/ all Primitive Drawing
c       4/ the GOURAUD interpolated shader
c
c      Dan Kidger 15-6-93
C-----------------------------------------------------------------------

C----------------------------------------------------------------------
        SUBROUTINE DR_PRIM (TEXT,X,Y,N,ICOL,IOP)
C
C  ******** This is the Main primitive drawing interface *********
C
C       ---> if ICOL .lt. 0  then skips
C
C       IOP   1 = for closed edges
c             2 = for filled
c             3 = unclosed edge ,etc.
c            10 = filled circle (eg. nodes)
C            20 = text string (eg node numbers) -- normal font. 
c
C --> here we go into INT*2 for DBOS drawing routines

      COMMON /DEST/IDEST,IOUT_D     !-- eg. 5=postscript

      REAL X(*),Y(*)        !--- x,y,(z) coords of the line/polygon
      INTEGER*2 XPI(20), YPI(20), IFSI, IFSP1,ICOLI
      CHARACTER TEXT*80

      IF (ICOL.lt.0) RETURN      ! ie. if color = -1 exit !
      ICOLI = ICOL
      IFSI  = N
      DO I=1,N                 ! copy real coords to integer*2 x,y lists 
        XPI (I) = NINT (X(I))   !-- (no need for PostScript :-)
        YPI (I) = NINT (Y(I))
      ENDDO

      IF (IOP.EQ.3) THEN       !----------simple line draw -------------
        IF (IDEST.EQ.2.AND.IOUT_D.EQ.5) THEN
          CALL DR_PS (TEXT,X,Y,N,ICOL, 1)
        ELSE
          CALL POLYLINE@ (XPI,YPI, IFSI, ICOLI)
        ENDIF
      ELSEIF (IOP.EQ.1) THEN   !--------closed polygon line draw--------
        IF (IDEST.EQ.2.AND.IOUT_D.EQ.5) THEN
          CALL DR_PS (TEXT,X,Y,N,ICOL, 2)
        ELSE
          IFSP1 = IFSI+1
          XPI (IFSP1) = XPI (1)
          YPI (IFSP1) = YPI (1)
          CALL POLYLINE@       (XPI,YPI, IFSP1, ICOLI)
        ENDIF
      ELSEIF (IOP.EQ.2) THEN   !--------- filled polygon draw ----------
        IF (IDEST.eq.2.and.IOUT_D.EQ.5) THEN
          CALL DR_PS (TEXT,X,Y,N,ICOL, 5)
        ELSE
          CALL FILL_POLY (XPI,YPI,IFSI,ICOLI)
        ENDIF
      ELSEIF (IOP.EQ.10) THEN   !------ filled circle (eg nodes) -------
        IF (IDEST.eq.2.and.IOUT_D.EQ.5) THEN
          CALL DR_PS (TEXT,X,Y,2,ICOL, 10)
        ELSE
          CALL FILL_ELLIPSE@ (XPI(1),YPI(1),XPI(2),XPI(2),ICOLI)
        ENDIF
      ELSEIF (IOP.EQ.20) THEN   !------ Text-string -------
        IF (IDEST.eq.2.and.IOUT_D.EQ.5) THEN
          CALL DR_PS (TEXT,X,Y,1,ICOL, 20)
        ELSE
          CALL DRAW_TEXT@ (TEXT,XPI(1),YPI(1),ICOLI)
        ENDIF
      ELSE
        PRINT*,'IOP=',IOP,' unknown'
      ENDIF
      END

C----------------------------------------------------------------------
      SUBROUTINE DRAW_POLY (SC,ISC,IFS,ICOL,IOP)
C
C     This is simply a header to DR_PRIM 
C     it unpacks a 2D list of coords into X() and Y() then calls DR_PRIM
C
C     IOP = 1 for edges, = 2 for filled
C
       SAVE                  !/ faster ? /!
       REAL SC(ISC,*)        ! x,y,(z) coords of the line/polygon
       REAL X(50), Y(50)
       CHARACTER TEXT*80
 
C     IF (ICOL.lt.0) RETURN      ! ie. if color = -1 exit !
c      ICOLI = ICOL
c      IFSI  = IFS

      DO I=1,IFS    
        X(I) = SC(I,1)
        Y(I) = SC(I,2)
      ENDDO

      CALL DR_PRIM (TEXT,X,Y,IFS,ICOL, IOP)   
      RETURN
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
C.... this is too slow ... should build the IROT, etc. bits
C.... directly into PM
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

C----------------------------------------------------------------------
      SUBROUTINE CAMERA(EYE,COA,UP, DEYE,FEYE,PM)
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

C ----- make sure are EYE,UP unit vectors ? -----
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
      SUBROUTINE XM2EYE (XM,YM,EYE,IOP)
C
C     to convert latitude/longitude (XM,YM) to vector EYE (or LIGHT)
C     - and also resets the XM to -180 -> +180,  YM = -180 -> +180
C      .. IF IOP /= 0 then XM, YM are reset ('ignored' ??)
C
      REAL EYE(*)
      DTR = 3.14159265/180.

      XM = MOD(XM+180.+360.,360.) -180.    !- to range +/- 180.
      YM = MOD(YM+180.+360.,360.) -180. 

      IF (IOP.NE.0) CALL GET_PRESET_VIEW (XM,YM,IOP)
      EYE(1) =  SIN (XM*DTR) 
      EYE(3) =  COS (XM*DTR)
      EYE(2) =  SIN (YM*DTR)
      EYE(1) = EYE(1) * COS(YM*DTR)
      EYE(3) = EYE(3) * COS(YM*DTR)
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE GET_PRESET_VIEW (XME,YME,IOP)
C
C     if IOP>0 then the view is set to a predefined direction
C
      IF (IOP.EQ.1) THEN         ! down y-axis
        XME =   0.
        YME =  89.999    ! (to avoid the singularity)
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
        ENDDO
        A(I,I) = 1.
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
c.... Don't use DRAW_POLY ?? (use DR_PRIM instead) cos' of PS output
c       CALL DRAW_POLY (SC,ISC,IFS,INTS(CC/IFS),2)
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
C-----------------------------------------------------------------------
      SUBROUTINE ISWAP (I,J)
c
c     this just swaps the 2 integer arguements
c 
        K = I
        I = J
        J = K
      END
c-----------------------------------------------------------------------
      SUBROUTINE GOURAUD_FILL (X,Y,N,ICOL,IOP, ICMAT)
C
C     This is a sub-part of 'GOURAUD' which simply draws the shapes
C     -- FILLED
C     IOP is the 'style' : black/coloured/zebra etc.
C     ITYPE =1 for filled and =2 for just the lines

      REAL X(*), Y(*)
      CHARACTER T*80

      IB = 2           !- offset into the colour table

      J = ICOL/2            !- Ok half the colour
      K = ICOL-2*J          !- either 1=odd, or 2=even
c     ICOL2 = IB + MOD (ICOL,13)    !- offset into the colour table
      ICOL2 = IB + ICOL             !- offset into the colour table
c... hmm so fill will never pick a colour outside 1->16

      IF(IOP.eq.0) THEN
        CALL DR_PRIM (T,X,Y,N,0,2)               !- black
      ELSEIF (IOP.eq.1) THEN
        CALL DR_PRIM (T,X,Y,N,1,2)               !- white
      ELSEIF (IOP.eq.2) THEN
        CALL DR_PRIM (T,X,Y,N,ICOL2,2)            !- coloured
      ELSEIF (IOP.eq.3) THEN
        CALL DR_PRIM (T,X,Y,N,MOD(ICOL2,2),2)     !- zebra evens
      ELSEIF (IOP.eq.4) THEN
        CALL DR_PRIM (T,X,Y,N,MOD(ICOL2+1,2),2)   !- zebra odds
      ELSEIF (IOP.eq.5) THEN
        IF (K.EQ.0) CALL DR_PRIM (T,X,Y,N,J+IB,2) !- stripe with invis
      ELSEIF (IOP.eq.6) THEN
        IF (K.EQ.0) CALL DR_PRIM (T,X,Y,N,J+1,2)
        IF (K.EQ.1) CALL DR_PRIM (T,X,Y,N,0,2)   !- stripe with 0
      ELSEIF (IOP.eq.7) THEN
        IF (K.EQ.0) CALL DR_PRIM (T,X,Y,N,J+1,2)
        IF (K.EQ.1) CALL DR_PRIM (T,X,Y,N,1,2)   !- stripe with 1
      ELSEIF (IOP.eq.8) THEN
        CALL DR_PRIM (T,X,Y,N,ICMAT,2)           !- mat colour .. why ?
      ELSEIF (IOP.eq.9) THEN
        IF (K.EQ.0) CALL DR_PRIM (T,X,Y,N,ICMAT,2) !- zebra with mat col
        IF (K.EQ.1) CALL DR_PRIM (T,X,Y,N,1,2)          ! (obs)
      ELSEIF (IOP.eq.10) THEN
        ICOL3= (15-ICOL)*16 +ICMAT+1               
        CALL DR_PRIM (T,X,Y,N,ICOL3,2)             !- shade material ?
      ENDIF
      END
c-----------------------------------------------------------------------
      SUBROUTINE GOURAUD_LINE (X,Y,N,ICOL,IOP, ICMAT)
C
C     This is a sub-part of 'GOURAUD' which simply draws the shapes
C     -- LINES
C     IOP is the 'style' : black/coloured/zebra etc.
c     ..  hmmm we don't actualy use IMAT do we ?
c
      REAL X(*),Y(*)
      CHARACTER T*80

      IB = 2
      ICOL2 = IB + ICOL             !- offset into the colour table
      J = ICOL/2
      K = ICOL-2*J
      IF (IOP.eq.0) THEN
         CALL DR_PRIM (T,X(N-1),Y(N-1),2,0,3)              !- black
       ELSEIF (IOP.eq.1) THEN                    
         CALL DR_PRIM (T,X(N-1),Y(N-1),2,1,3)              !- white
       ELSEIF (IOP.eq.2) THEN                        
         CALL DR_PRIM (T,X(N-1),Y(N-1),2,ICOL2,3)           !- coloured
       ELSEIF (IOP.eq.3) THEN
         CALL DR_PRIM (T,X(N-1),Y(N-1),2,MOD(ICOL,2),3)    !- zebra:1
       ELSEIF (IOP.eq.4) THEN
         CALL DR_PRIM (T,X(N-1),Y(N-1),2,MOD(ICOL+1,2),3)  !- zebra:2
       ENDIF
      END

C-----------------------------------------------------------------------
      SUBROUTINE FILL_POLY (XI,YI,N,ICOL)
C                                                            DJK 19-3-92
C     this colours-in a polygon (integer co-ords)            
C     ... with co-incident vertex elimination :-)
C    (cf the Postscript version)
C
      INTEGER*2  XI(*), YI(*) ,N,ICOL       !- all input as INT*2

      INTEGER*2  X(30), Y(30),IFAIL,NN,IPOLY,ICOLS, iymin,iymax

c    ... killing ICOL=-1 is handled outside in DR_PRIM
c     IF (ICOL.NE.-1) then  ! skip 'invisible colour = -1'
C--------- patch to remove 'zero-height' polygons 14-07-92 -------------
c.. I thought the error was just co-incident nodes ??
        ICOLS = ICOL
        NN = N
        iymin = xi(1)
        iymax = yi(1)
        DO I=1,N            !- why the copy here ??
          X(I) = XI(I)
          Y(I) = YI(I)
        ENDDO
        DO I=1,N
          iymin = min(iymin,y(i))
          iymax = max(iymax,y(i))
        ENDDO
        IF (iymin.ne.iymax) THEN
C--------------------------------------------------------------------

        CALL CREATE_POLYGON@      (X,Y,NN,IPOLY,  IFAIL)
        CALL DOSERR@ (IFAIL)
        CALL FILL  _POLYGON@            (IPOLY,ICOLS,IFAIL)
        CALL DOSERR@ (IFAIL)
        CALL DELETE_POLYGON_DEFINITION@ (IPOLY,  IFAIL)
        CALL DOSERR@ (IFAIL)

c.. debug..info
c         write(19,'(50(''-''))')
c        write(19,'(5(''  <>'',2i4))') (x(i),y(i),i=1,nn),icol

        ENDIF
c     ENDIF
      RETURN
      END

C----------------------------------------------------------------------
      SUBROUTINE DR_ARROW (pt_fr,pt_to,iangle,iscale,icol)
C
C      to draw an arrow vector (from L2.for)    c. Dan Kidger
C       eg. angle=10. (deg) , scale = 10 (%)
C
c.... nicer as a sub-function of DR_PRIM ..extends X and Y for the head! 
c.... also could have 'other' arrows eg. filled_head
c
      REAL PT_FR(*),PT_TO(*) ,  x(6),y(6)
      CHARACTER T*80

      XL=pt_to(1) - pt_fr(1)
      YL=pt_to(2) - pt_fr(2)

      IF (ISCALE.GT.0) THEN           ! skip if no arrow head
        A  =  3.14159295/180  * IANGLE 
        SA = -SIN(A) * ISCALE/100.
        CA = -COS(A) * ISCALE/100.

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

        CALL DR_PRIM (T,X,Y,5,ICOL,3)    !----- as line-draw (filled ?)

      ENDIF
      RETURN
      END
C----------------------------------------------------------------------
      SUBROUTINE GRATIC (IRESX,IRESY,NX,NY,ICOL)
C
C     draws a graticule of RESX by RESY boxes
C .. should really draw within the current window ??
C     15-6-93   ... actually this is a bit obsolete ! 

      INTEGER*2 IRESX,IRESY
      CHARACTER T*80
      REAL X(2),Y(2)
      JRESX = IRESX-1         !-- so that the far side is visible
      JRESY = IRESY-1         !--  "    "   " "   "
      DO I=0,NX                        !--- horizontal rulings
        X(1) = 0.
        X(2) = JRESX
        Y(1) = REAL(JRESY)*I / REAL(NX)
        Y(2) = Y(1)
        CALL DR_PRIM (T,X,Y,2,ICOL,3)    !----- edges
      ENDDO

      DO I=0,NY                         !--- vertical rulings       
        X(1) = REAL(JRESX)*I / REAL(NY)
        X(2) = X(1)
        Y(1) = 0.
        Y(2) = JRESY
        CALL DR_PRIM (T,X,Y,2,ICOL,3)    !----- edges
      ENDDO
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE SELECT_PS_PRINTER (FILE,IRES_X,IRES_Y
     +       ,ioff_x,ioff_y,file_dat)
C   
C     This opens a file for post-script output on unit 55
C       
      CHARACTER FILE*(*) ,FILE_DAT*(40), FILE2*(40)
     +   ,TIME@*8,EDATE@*8, FDATE@*20
      EXTERNAL TIME@,EDATE@,FDATE@
      OPEN (55,FILE=FILE)
      FILE2 = FILE_DAT
      DO I=1,LENG(FILE2)     !- reverse backslash for postscript output
        IF (FILE2(I:I).EQ.'\') FILE2(I:I)='/'
      ENDDO

C----------------------- Postscript headers ----------------------------
      WRITE(55,'(A)')
     + '%!Adobe-2.0 EPSF'
     +,'%%Creator: Danplot Finite Element Postprocessor'
     +,'%%Creationdate '//TIME@()//' '//FDATE@()
     +,'%%Title: Unknown'
     +,'%%BoundingBox 0 0 576 792'

c------------------------------ Macros ---------------------------------

      WRITE(55,'(A)') ' ' 
     +,'/fp { newpath 3 1 roll moveto {lineto} repeat'//
     +    ' closepath setrgbcolor fill } def' 
     +,'/dl { newpath 3 1 roll moveto {lineto} repeat'//
     +    ' setrgbcolor stroke } def' 
     +,'/dp { newpath 3 1 roll moveto {lineto} repeat'//
     +    ' closepath setrgbcolor stroke } def' 
     +,'/fc { newpath pop 4 2 roll pop 0 360 arc'//
     +              ' setrgbcolor fill } def'
     +,'/dt { /Helvetica findfont 30 scalefont setfont 6 1 roll'//
     +    ' moveto setrgbcolor show } def'

c    +,' ' 
c    +,'gsave clippath pathbbox grestore'  !-----use Landscape mode ?----
c    +,'4 dict begin'
c    +,'/ury exch def /urx exch def /lly exch def /llx exch def'
c    +,'90 rotate      ury urx sub llx ury add neg translate'
c    +,'end'
c    +,' ' 
c    +,'newpath llx lly moveto urx lly lineto urx ury lineto llx ury '//
c    + 'moveto closepath gsave 20 setlinewidth stroke grestore'

C----------------------- Logo + date + filename ------------------------
c.. do the logo with *style*, then filename centered then time and date
      WRITE(55,'(A)') ' ' 
     +,'gsave /Times-BoldItalic findfont 10 scalefont setfont'
     +,'/rays { 0 1.5 179 {gsave rotate 0 0 moveto 108 0 lineto stroke'
     +//  ' grestore  } for } def'
     +,'28 50 translate .25 setlinewidth  1 1 scale'
     +,'newpath 0 0 moveto (Danplot v2.0) true charpath clip'
     +,'newpath 54 -15 translate rays grestore' 

      WRITE(55,'(A)') ' ' 
     +,'/Helvetica findfont 10 scalefont setfont 575 2 div 45 moveto'
     +,' ('//FILE2(1:LENG(FILE2))//') dup '//
     +  'stringwidth pop 2 div 0 exch sub 0 rmoveto show '

      WRITE(55,'(A)') ' ' 
     +,'/Helvetica findfont  6 scalefont setfont 575 50 moveto'
     +,' ('//TIME@()//' '//EDATE@()//') dup '//
     +  'stringwidth pop 0 exch sub 0 rmoveto show'

c-------------------------- set scaling --------------------------------

      WRITE(55,'(A)') ' ' 
     +,' 72 300 div dup scale '      !--- units in 300ths inch
     +,' ' 


      IRES_X =   300. * ( 8. -.4)
      IRES_Y =   300. * (11. -.4)
      IOFF_X =   300. * (   .4  ) 
      IOFF_Y =   300. * (   .8  )

      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE CLOSE_PS_PRINTER ()
C   
C     This closes the post-script output on unit 55
C       
      WRITE(55,'(A)')
     +       ' grestore '            !- back out of rescaling
     +      ,'%%   "showpage" will be inhibited for EPS files'
c    +       '%%   finally show the page'  
c    +      ,'showpage'
      CLOSE(55)
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE DR_PS (TEXT,X,Y,N,ICOL, IOP)
C
C     This writes PostScript Output to unit 55 for lines,fill_areas etc.
C     IOP = 1 un-closed line
C           2    closed line = polygon
C           5   filled polygon
C          20   text string   (auto font & scale)
c
c..... 29-11-93 

      REAL X(*),Y(*)
      INTEGER N, ICOL,IOP
      COMMON /PALETTE/PAL
      INTEGER PAL(3,0:255)
      CHARACTER FORMAT*20, ACTION*5, TEXT*80

      IF (IOP.EQ. 1) ACTION = ' dl'        !/* draw-line */
      IF (IOP.EQ. 2) ACTION = ' dp'        !/* draw-polygon */
      IF (IOP.EQ. 5) ACTION = ' fp'        !/* fill-polygon */
      IF (IOP.EQ.10) ACTION = ' fc'        !/* fill-circle */
      IF (IOP.EQ.20) ACTION = ' dt'        !/* draw-text */

      WRITE (FORMAT,'(A,I2,A)') '(3f6.3,',2*N,'I7,I3,A)' 
      IF (IOP.EQ.20) THEN    !- handle text differently
        WRITE (FORMAT,'(A,I2,A)') '(3f6.3,',2*N,'I5,4A)' 
        WRITE (55,FORMAT) 
     +        (PAL(J,ICOL)/255.,J=1,3) 
     +       ,INT(X(1)),INT(Y(1))
     +       ,' (',TEXT(1:LENG(TEXT)),') ',ACTION

      ELSE          !- standard 'default' format
        WRITE (55,FORMAT) 
     +          (PAL(J,ICOL)/255.,J=1,3)
     +         ,(INT(X(I)),INT(Y(I)),I=N,1,-1), N-1, ACTION
      ENDIF
      RETURN
      END

C-----------------------------------------------------------------------
c     SUBROUTINE PCX_FIX_old (FILE_PCX)
C
C        *OBSOLETE --> only handles 16 colour files
C     this 'corrects' the palette in a PCX file to be the same as
C     on the screen .. but only for 16 colour files tho' :-(
C
c.... the 256 colour version is in c:\users\dan\2\ssc4.for ??
c
c     COMMON /PALETTE/PAL
c     INTEGER PAL(3,0:255)
c     INTEGER*2 IFAIL,IO
c     INTEGER*4 NBYTES_R
c     CHARACTER FILE_PCX*(*), dummy*16
c     CALL OPENRW@ (FILE_PCX,IO,IFAIL)
c     CALL READF@  (dummy,IO,16L,NBYTES_R,IFAIL)
c     DO I=0,15
c       CALL WRITEF@ (CHAR(PAL(1,I)),IO,1L,IFAIL)
c       CALL WRITEF@ (CHAR(PAL(2,I)),IO,1L,IFAIL)
c       CALL WRITEF@ (CHAR(PAL(3,I)),IO,1L,IFAIL)
c     ENDDO
c     CALL CLOSEF@ (IO,IFAIL)
c     END

C-----------------------------------------------------------------------
      SUBROUTINE PCX_FIX (FILE_PCX,PAL)
C
C     this 'corrects' the palette in a PCX file to be the same as
C     on the screen
C       --> also handles 256 colour PCX files
C        ( PAL of course does not have to come from COMMON !)
C
c     COMMON /PALETTE/PAL

      INTEGER*2  PAL(3,0:255)
      CHARACTER FILE_PCX*(*)

      CHARACTER dummy *(16)
      INTEGER*4  NBYTES_R, EOF, NEW_POS, POS, nb

      CALL OPENRW@ (FILE_PCX,IO,IERROR)

      nb = 16
      CALL READF@ (dummy,IO,nb,NBYTES_R,IERROR)   !- why not just file_pos

c     CALL FPOS@ (IO,nb,NEW_POS,IERROR)     !- ok at pos 16 (or 17?)

      DO I=0,15              !- the first copy of the palette
        nb = 1
        CALL WRITEF@ (CHAR(PAL(1,I)),IO,nb,IERROR)
        CALL WRITEF@ (CHAR(PAL(2,I)),IO,nb,IERROR)
        CALL WRITEF@ (CHAR(PAL(3,I)),IO,nb,IERROR)
      ENDDO

      IF (ICHAR(DUMMY(4:4)).eq.8) then     ! 8 bit planes .. so paletee 
c                                          ! is at the end too !
        POS = 1024
        POS = POS * POS                    !(this avoids int*2 problems)
        CALL FPOS@(IO,POS,NEW_POS,IERROR)
        EOF = NEW_POS
        POS = EOF - 769
        CALL FPOS@(IO,POS,NEW_POS,IERROR)
        nb = 1
        call READF@(DUMMY,IO,nb,NBYTES_R,IERROR)
        IF(ICHAR(DUMMY(1:1)).ne.12) THEN
          print*,' *** Palette at the end of the PCX file not found'
          return
        ENDIF
        DO I=0,255          !(can't I write 3 bytes at once?
          nb = 1
          CALL WRITEF@ (CHAR(PAL(1,I)),IO,nb,IERROR)
          CALL WRITEF@ (CHAR(PAL(2,I)),IO,nb,IERROR)
          CALL WRITEF@ (CHAR(PAL(3,I)),IO,nb,IERROR)
        ENDDO
      ENDIF
      CALL CLOSEF@ (IO,IFAIL)         !---- OOPS ! 6-11-93
      END
C----------------------------------------------------------------------

C-----------------------------------------------------------------------
      SUBROUTINE INTO_GRAPHICS ()
C
C     This goes into standard VGA mode 
C     Then sets the desired palette  ( black,white,gray, colours..)
C     and then enables the mouse
C      .. just rest VGA and then set_all_dacs ??
c
      LOGICAL LMOUSE
      INTEGER*2 IRESX,IRESY, ICOL_S

      call vga@()

C------------------ set the palette to be the video DACs ---------------
      DO I=1,15      
        IPOS_S = I
        ICOL_S = I
        CALL SET_PALETTE@ (IPOS_S,ICOL_S)
      ENDDO
      IPOS_S = 17
      ICOL_S = 2
      CALL SET_PALETTE@ (IPOS_S,ICOL_S)       !- overscan
c-----------------------------------------------------------------------
c      CALL SET_PAL ( 0,0,0,0,13)    !- set all PAL to VGA dacs !
c      CALL SET_PAL (11,0,0,0,13)    !- old 'default' set ? 
c      CALL SET_PAL (10,0,0,0,13)    !- default menus
c     CALL SET_PAL (14,0,0,0,13)    !- spectrum set ..  try this one
c     call set_all_dacs()
c     call graphics_mode_set@ (640s,480s,256L,ifail)   !- this fails!
c-----------------------------------------------------------------------

      IRESx = 640
      IRESy = 480

C----------------------- reset mouse -----------------------------------
c.. use generic_int_to int*2 ??
      CALL INITIALISE_MOUSE@ ()        !
      CALL MOUSE_SOFT_RESET@ (LMOUSE)  ! =.false. if mouse is not present
      CALL SET_MOUSE_BOUNDS@ (0,0,IRESx-1,IRESy-1)
      CALL DISPLAY_MOUSE_CURSOR@()     !
      CALL SET_MOUSE_POSITION@ ( 121 + 10, 41+ 10 )
c-----------------------------------------------------------------------

      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE RESET_SIZES (RES)
C
C     This simply resets the image size  to the 'drawing' window
C
      INTEGER RES(*)

      RES(4) = 121     !- offsets -x
      RES(5) = 41      !          -y
      RES(1) = 518     !- image width    !- is this right ?
      RES(2) = 413     !- image height   !- hmm ?
      RETURN
      END
C-----------------------------------------------------------------------
c     SUBROUTINE DEBOUNCE_BUTTONS
C
C     This subroutine will absorb all mouse button presses
C     until there is some change in the button status.
C
c     INTEGER*2 OLD_STATUS,IH,IV,STATUS
c     CALL GET_MOUSE_POS@(IH,IV,OLD_STATUS)
c     STATUS=OLD_STATUS
c     WHILE(STATUS.EQ.OLD_STATUS.AND.STATUS.NE.0)DO
c        CALL GET_MOUSE_POS@(IH,IV,STATUS)
c        ENDWHILE
c     END

C----------------------------------------------------------------------
      SUBROUTINE SET_PAL (IOP,ir,ig,ib,ncol)
C
C     set palette option IOP:  1 = simple, 2= greyscale, etc.
C
C .. SET COLORS 1-5 = black,white,d.green,l.green, brown,red
C ... -1 = user defined , -2 = 'standard' ,-3 = 'firestorm'
C      0..16 = normal colors
C----- also handles 'modifiers' eg. EVENS, REVERSE,INVERT ...
C... be VERY careful about the set of colors that this handles
C... ie. contour colours/ menu colours/ black+white  etc.

      COMMON /PALETTE/PAL
      save imenu_back
      INTEGER*2 IFAIL                   !- for read_edited_line ?
      INTEGER PAL (3,0:255)
      INTEGER PAL_DEFAULT(3,0:15),      !-  'old' Danplot colours
     +        PAL_RAINBOW(3,0:15),      !-  a preset rainbow 
     +        PAL_Menupane(3,1:6)       !- some background colours
      INTEGER*1 PAL1 (3,0:255)
      CHARACTER LINE*80 
      REAL PI
      DATA PAL_DEFAULT/  0,  0,  0,   255,255,255,   100,100,180, 
     +    0, 90, 10,   127,127,  0,   200, 20, 20,    90,  0,190,
     +  160,  0,160,   127,127,127,    63, 63, 63,
     +    0,255,  0,    40, 40,255,   255, 40,255,   255,  0,  0,
     +  255,127,  0,   255, 255, 0/

      DATA PAL_RAINBOW/  0,  0,  0,   200,255,200,    40,120, 90, 
     +    1,177,204,     1,204,177,     9,227,147,    23,243,116,
     +   45,253, 85,    71,255, 57,   100,249, 33,   131,236, 15,
     +  162,216,  4,   191,191,  0,   216,162,  4,   236,131, 15,
     +  249,100, 33/

      DATA imenu_back / 0 /
      DATA PAL_Menupane/  160, 80, 80,  80,160, 80,  80, 80,160,  
     +      160,160, 80,   80,160,160, 160, 80,160 /
 
      DATA PI/3.14159265/
      JM = 2          !- #menu colours: so COLORS  0,1,2  are reserved
C----------------------------------------------------------------------
c      DO I=1,15                        !- maybe only do once ?
c        CALL SET_PALETTE@ (INTS(I),INTS(I))
c      ENDDO
c      CALL SET_PALETTE@ (17,2)       !- overscan

C----------------------------------------------------------------------
C--------- rest all colours to 'standard' VGA ------------
      IF (IOP.EQ.0) THEN
        CALL GET_VIDEO_DAC_BLOCK@ (1,256,PAL1)
        DO I=0,255
          DO J=1,3
            PAL (J,I) = PAL1 (J,I) * 4
          ENDDO
        ENDDO
      ELSEIF (IOP.GE.1000.AND.IOP.lt.1256) THEN
C--------- single colour redefinition     (obs) ??  ------
        IVAL = IOP -1000
        PAL (1,IVAL) = IR
        PAL (2,IVAL) = IG
        PAL (3,IVAL) = IB

C-------------- reverse palette ---------------
      ELSEIF (IOP.EQ.2) THEN
        DO I=1,(ncol+1)/2 
          DO J=1,3                        !-- be careful of ISWAP !
            CALL ISWAP (PAL(J,I+jm),PAL(J,ncol+jm+1-I) )
          ENDDO
        ENDDO

c---------------------- invert all the colours ---------------------
      ELSEIF (IOP.EQ.3) THEN    
        DO I=1,ncol
          DO J=1,3
            PAL(J,jm+I) = 255 - PAL(J,jm+I)
          ENDDO
        ENDDO

C------------- reverse black and white --------------
      ELSEIF (IOP.EQ.4) THEN
        DO J=1,3
          CALL ISWAP (PAL(J,0),PAL(J,1) )
        ENDDO

C----------- reverse alternate colours --------------
      ELSE IF(IOP.EQ.5) THEN                                       
        DO I=1,ncol-1,2               !- == 'evens'
          DO J=1,3
            CALL ISWAP (PAL(J,I+1),PAL(J,I+2) )
          ENDDO
        ENDDO

c-----------------------------------------------------------------------
      ELSEIF (IOP.EQ.201) THEN
C .. shade range redefinition       <- what is this ?  'Shade' ??
        IOFF= 6
        DO I=1,ncol
          PAL(1,I+JM) = IR * (I+IOFF) /(NCOL+IOFF)
          PAL(2,I+JM) = IG * (I+IOFF) /(NCOL+IOFF)
          PAL(3,I+JM) = IB * (I+IOFF) /(NCOL+IOFF)
        ENDDO

c-----------------------------------------------------------------------
c--------------------  'preset' palettes -------------------------------

c--------- standard 'menu' colours ------
      ELSEIF (IOP.EQ.10) THEN
        DO I=0,2
          DO J=1,3
            PAL (J,I) = PAL_DEFAULT (J,I)
          ENDDO
        ENDDO

c....... the 'fudgy bit' for different menu background colours
c... ARRGH! .. so the menu colours will keep on changing !
         imenu_back = mod (imenu_back, 6) + 1 
          DO J=1,3
            PAL (J,2) = PAL_Menupane (J,imenu_back)
          ENDDO

c--------- standard 'material' colours ------
      ELSEIF (IOP.EQ.11) THEN
        DO I=3,15
          DO J=1,3
            PAL (J,I) = PAL_DEFAULT (J,I)
          ENDDO
        ENDDO

C---------- 'rainbow' colours for contouring ----------
      ELSEIF (IOP.EQ.12) THEN
        DO I=3,15               !- 'miss' menu cols
          DO J=1,3
            PAL(J,I) = PAL_RAINBOW(J,I)
          ENDDO
        ENDDO

C---------- 'Hot-Iron' colours for contouring ----------
      ELSEIF (IOP.EQ.13) THEN          != black->red->yellow->white
        DO I=1,ncol
          hue = real(I) / real(ncol)   ! hue = 0.-->1.
          PAL(1,I+jm) = min(max(0,nint(255.* 3.*(hue+.03    ))),255) ! note 
          PAL(2,I+jm) = min(max(0,nint(255.* 3.*(hue-.333333))),255) ! clipping
          PAL(3,I+jm) = min(max(0,nint(255.* 3.*(hue-.666667))),255)
        ENDDO

c---------- 'spectrum' colours for contouring -------
      ELSEIF (IOP.EQ.14) THEN                                       
        DO I=0,ncol-1
          RED = real(I) /real(ncol) * 2 * PI
          PAL (1,jm+1+I) = 128. + 127. * COS (RED)
          PAL (2,jm+1+I) = 128. + 127. * COS (RED + 2.* PI/3.)
          PAL (3,jm+1+I) = 128. + 127. * COS (RED - 2.* PI/3 )
        ENDDO

C---------- 'greyscale' colours for contouring ----------
      ELSEIF (IOP.EQ.15) THEN   
        from = .5
        to = 1.
        DO I=1,ncol
          val = from + (to-from) * real(i)/real(ncol)
          PAL(1,jm+I) = 255. * val
          PAL(2,jm+I) = 255. * val
          PAL(3,jm+I) = 255. * val
        ENDDO

c------------------ 'linear' palettes ----------------------------------
      ELSEIF (IOP.EQ.21) THEN   
        LINE = '3 15 1.     0 .5  0    1 0 0 '       ! green--> red     
        CALL SET_SHADE (LINE)

      ELSEIF (IOP.EQ.22) THEN                   
        LINE = '3 15 1.     0 0 .7     1 1 0 '       ! bure -> yellow
        CALL SET_SHADE (LINE)

      ELSEIF (IOP.EQ.23) THEN                        ! blue->black->red
        LINE = '3  9 1.     0 .5  1    0 0 0 '
        CALL SET_SHADE (line)
        LINE = '9 15 1.     0  0  0    1 0 0 '
        CALL SET_SHADE (LINE)

      ELSEIF (IOP.EQ.25) THEN                   
        LINE = '3 15 1.     .1 .1 .1     1. 1. 1. '       ! gray
        CALL SET_SHADE (LINE)

      ELSEIF (IOP.EQ.24) THEN                        ! blue->black->red
        LINE = '3 15 1.     0 0 .7    1 1 0 '
        CALL READ_EDITED_LINE@ (LINE,0,1,33,IFAIL)   !- 'esc' to quit
        IF (IFAIL.EQ.0) CALL SET_SHADE (LINE)     


C------------------- palette 'animations' etc. -------------------------

C----------- step left -------------
      ELSEIF (IOP.EQ.41) THEN
        DO J=1,3
          DO I=15,jm+1,-1 
            PAL(J,I) = PAL(J,I-1)
          ENDDO
            PAL(J,jm+1) = PAL(J,15)
        ENDDO

C----------- step right -------------
      ELSEIF (IOP.EQ.42) THEN
        DO J=1,3
          DO I=jm+1,15-1 
            PAL(J,I) = PAL(J,I+1)
          ENDDO
            PAL(J,15) = PAL(J,2)
        ENDDO

c-----------------------------------------------------------------------
      ENDIF     !-- end-of-all options

      IF (NCOL.EQ.13) THEN   !-- shade down the other colors :-)
        DO I=0,15            !-- if we do this then 256 colour
          DO K=0,15          !-- mode will show the lighting-shading
            I2 = K*16 + I
            DO J=1,3   
              PAL(J,I2) = PAL(J,I) * (16+5-K) / (16.+5.)
            ENDDO
          ENDDO
        ENDDO
      ENDIF  
      CALL SET_ALL_DACS ()
      END

C-----------------------------------------------------------------------
      SUBROUTINE SET_SHADES (IOP)
C
C     This sets one of a set of pre-defined (or user-def) linear palettes
C     .. potentialy I can have a large set of these
C     .. Is this obsolete ? (cf SET_PAL itself)
C
      CHARACTER LINE*80
      INTEGER*2 IFAIL                          ! for read_edited_line

      IF (IOP.eq.1) THEN                       ! blue --> yellow
         LINE = '2 15 1.     0 0 .7   1 1 0 '
         CALL SET_SHADE (LINE)
      ELSEIF (IOP.eq.2) THEN                   ! green--> red
         LINE = '2 15 1.     0 .5 0   1 0 0 '
         CALL SET_SHADE (LINE)
      ELSEIF (IOP.eq.3) THEN                   ! blue->black->red
         LINE = '2  8 1.     0 1 1    0 0 0 '
         CALL SET_SHADE (line)
         LINE = '8 15 1.     0 0 0    1 0 0 '
         CALL SET_SHADE (LINE)
      ELSE                                      !- user-defined
         LINE = '2 15 1.     0 0 .7    1 1 0 '
         CALL READ_EDITED_LINE@ (LINE,0,1,33,IFAIL)   !- 'esc' to quit
         IF (IFAIL.EQ.0) CALL SET_SHADE (LINE)     
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SET_SHADE (LINE)
C
C     this sets a shade range palette read-in from a CHARACTER LINE
C
      CHARACTER LINE*(*)
      COMMON /PALETTE/PAL
      INTEGER PAL(3,0:255)

      READ (LINE,*) ICF,ICT,GAMA, RF,GF,BF,  RT,GT,BT
c  ( gama correction not yet implimented )
      DO I = ICF,ICT
        FACT = REAL (ICT - I ) / REAL (ICT-ICF)
        PAL(1,I) = (RF * FACT + RT * (1.-FACT)) * 255
        PAL(2,I) = (GF * FACT + GT * (1.-FACT)) * 255
        PAL(3,I) = (BF * FACT + BT * (1.-FACT)) * 255
      ENDDO
      CALL SET_ALL_DACS ()
      END


C----------------------------------------------------------------------
      SUBROUTINE SET_ALL_DACS()
C... sets all the video dacs ...      16 or 256 colours ???
      COMMON /PALETTE/PAL
      INTEGER PAL(3,0:255)
      INTEGER*1 PAL1 (3,0:255)

        DO I=0,255
          DO J=1,3   
            PAL1(J,I) = PAL(J,I)/4 
          ENDDO
        ENDDO
        CALL SET_VIDEO_DAC_BLOCK@( 0, 256, PAL1)

      END
C----------------------------------------------------------------------
c-----------------------------------------------------------------------

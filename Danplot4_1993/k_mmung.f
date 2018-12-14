c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c             part of DANLIB ... MESH_GENERATION
c
c    Keyword droven modules ... << including the amazing *DANBLOCKS >>
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      SUBROUTINE KEY_MESH_MUNG
     +          (FOUND,KEYWORD,IO,GC,IGC,NDIM,NN,NUMS,INUMS,NEL,P)
C
C     A Keyword handler for reading/writing the basic mesh
C     DJK 24-5-93
C
      REAL GC (IGC,*)
      INTEGER NUMS(INUMS,*) , P(*)
      LOGICAL FOUND
      CHARACTER KEYWORD*(*)

      FOUND = .TRUE.   

c-----------------------------------------------------------------------

      IF (KEYWORD.EQ.'*DANBLOCKS') THEN
        CALL DANBLOCKS (IO,GC,IGC,NN,NDIM,NUMS,INUMS,NEL,P)

C---------------------- removing elements ------------------------------
      ELSEIF (KEYWORD.EQ.'*DELETE_MATERIALS') THEN
        CALL DEL_MATS (IO,GC,IGC,NN,NUMS,INUMS,NEL,P,NDIM)

      ELSEIF (KEYWORD.EQ.'*DELETE_ZERO_AREA_ELEMENTS') THEN
        CALL DEL_ZERO_AREA (IO,GC,IGC,NN,NUMS,INUMS,NEL,P,NDIM)

      ELSEIF (KEYWORD.EQ.'*DELETE_ORPHAN_NODES') THEN
        CALL DEL_O_NODES (IO,GC,IGC,NN,NUMS,INUMS,NEL,P,NDIM)
                                                
      ELSEIF (KEYWORD.EQ.'*SORT_NODES_FROM_POINT') THEN
        CALL SORT_NODES (IO,GC,IGC,NN,NUMS,INUMS,NEL,P,NDIM)
                                                
      ELSEIF (KEYWORD.EQ.'*DELETE_COINCIDENT_NODES') THEN
        CALL DEL_COIN_NODES (IO,GC,IGC,NN,NUMS,INUMS,NEL,P,NDIM)

C----------------------- modifing elements -----------------------------
c.. DANBLOCKS is so important it should appear at the top !

      ELSEIF (KEYWORD.EQ.'*4NQS_FROM_4_3NTS') THEN   !-- 'teapot.g'
      CALL MUNG_4NQ  (IO,GC,IGC,NN,NDIM,NUMS,INUMS,NEL,P)

      ELSEIF (KEYWORD.EQ.'*MIRROR_EVEN_ELEMENTS') THEN
      CALL MIRROR_EVENS (IO,GC,IGC,NN,NDIM,NUMS,INUMS,NEL,P)

      ELSEIF (KEYWORD.EQ.'*CHANGE_MATERIALS') THEN
        CALL CHANGE_MATS (IO,GC,IGC,NN,NUMS,INUMS,NEL,NDIM)
c----------------------- moving the mesh -------------------------------
      ELSEIF (KEYWORD.EQ.'*MOVE_NODES') THEN
        CALL R_NODMOV(IO,GC,IGC,NDIM,NN)

      ELSEIF (KEYWORD.EQ.'*SHIFT_MESH') THEN
        CALL SHIFT_MESH (IO,GC,IGC,NN,NDIM)

      ELSEIF (KEYWORD.EQ.'*SCALE_MESH') THEN
        CALL SCALE_MESH (IO,GC,IGC,NN,NDIM)

      ELSEIF (KEYWORD.EQ.'*X_TO_Y_TO_Z') THEN
        CALL X2Y2Z  (GC,IGC,NN,NDIM)

      ELSEIF (KEYWORD.EQ.'*X_TO_Y') THEN
        CALL X2Y    (GC,IGC,NN,NDIM)

      ELSEIF (KEYWORD.EQ.'*WRAP_AROUND_Y') THEN
        CALL WRAP_Y (GC,IGC,NN,NDIM)

      ELSEIF (KEYWORD.EQ.'*CIRCLE_A_SQUARE') THEN
        CALL WRAP_CIRCLE_SQUARE (IO,GC,IGC,NN,NDIM)

      ELSEIF (KEYWORD.EQ.'*WRAP_CIRCLE_Y') THEN
        CALL WRAP_CIRCLE_Y (GC,IGC,NN,NDIM)

      ELSEIF (KEYWORD.EQ.'*WRAP_SPHERE') THEN
        CALL WRAP_SPHERE (GC,IGC,NN,NDIM)

      ELSEIF (KEYWORD.EQ.'*WRAP_TWIST_Y') THEN
        CALL WRAP_TWIST_Y (IO,GC,IGC,NN,NDIM)

      ELSE
        FOUND = .FALSE.
      ENDIF
      RETURN
      END
C----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE DANBLOCKS (IO,GC,IGC,NN,NDIM,NUMS,INUMS,NEL,P)
C
C     This will 'condense' all nodes that share the same position
C     using FIND_NODE to get the node number of a node
C     3-4-93  now can fill 8nq's with 20nb !!  :-)
C     ... however first we need to set the z-dimensions to zero
C     ... then upadte NDIM to 3
C---> really there should be a 'daughter' subroutine where we just
C---- handle one element at a time

      PARAMETER (MAX_NOD = 27,MDIM=4)
      PARAMETER (ICOORD= MAX_NOD,IDER=MDIM)

      REAL      GC (IGC,*)        !- Nodal coord database
     +     ,WIDTHS (0:20,5)       !- input sub-element widths
     +  ,COORD_OLD (ICOORD,MDIM)  !- coords of the parent element
     +    ,COORD_L (ICOORD,3)     !-  local coords of the daughter elements
     +    ,COORD_G (MDIM)         !- global coords of a daughter node
     +       ,SAMP (1,MDIM)       !- sampling point of the 'new' nodes
     +        ,FUN (MAX_NOD)      !- shape functions of the parent
     +        ,DER (IDER,MAX_NOD) !- sf derivs of the parent (unused)

      INTEGER  NUMS (INUMS,*)     !- Element database
     +         ,NUM (MAX_NOD)     !- node numbers of the new elements
     +           ,P (*)           !- node pointers (unused)
     +           ,L (MDIM)        !- toggles for the sub-directions
     +         ,NXE (MDIM)        !- number of sub-elements in each direction

      NEL_NEW  = NEL      !- new element numbers (incremented)
      NN_NEW   = NN       !- new node numbers (updated)

      DO ILOOP = 1,9999

    1 READ (IO,*,IOSTAT=IOS) JMAT , NDIME_NEW,NOD_NEW,ITYPE  
      CALL IN_TEST (IO,IOS,*1,*999)

c------- go into '3d'  if we were only in '2d' before -------
        IF (NDIME_NEW.GT.NDIM) THEN
          DO I=1,NN
            DO J=NDIM+1,NDIME_NEW
              GC(J,I) = 0.
            ENDDO
          ENDDO
          NDIM = NDIME_NEW
        ENDIF

c-----------  get the local coords of the new element type -------------
      CALL WTHATN (NOD_NEW,NDIME_NEW,ITYPE,COORD_L,ICOORD)

      DO I=1,NDIME_NEW       !-- read the edge spacings
    2   READ (IO,*,IOSTAT=IOS)   NXE(I),(WIDTHS(J,I),J=1,NXE(I))
        CALL IN_TEST (IO,IOS,*2,*999)   

        WIDTHS(0,I) = 0.    !- null the first coord in the list

C.... careful .. I dont want to 'non-dimensionalize' if 2d-> 3d

        WIDTH_TOT = 0.
        DO J=0,NXE(I)    !-- get cumultive widths
          WIDTH_TOT   = WIDTH_TOT + WIDTHS(J,I)
          WIDTHS(J,I) = WIDTH_TOT 
        ENDDO
        DO J=0,NXE(I)   
          WIDTHS(J,I) = WIDTHS(J,I) / WIDTH_TOT !-- scale to 0. -> +1.
          WIDTHS(J,I) = WIDTHS(J,I) * 2. -1.    !-- into locals -1 -> +1
        ENDDO
      ENDDO

C-------------------- loop the 'original' elements ---------------------
      ic = 0
      DO IEL = 1,NEL
        write (*,'(i5,a,i5,a)')  iel,' /', nel,char(27)//'[A'
        CALL GET_ELEMENT (NUMS,INUMS,IEL, NUM, 
     +          NOD_OLD,NDIME_OLD,ITYPE,IMAT,IUSER1,IUSER2,IUSER3)

        IF (JMAT.EQ.IMAT) THEN               !-- found a material
          ic = ic + 1
          CALL GET_COORD (NUM,NOD_OLD,NDIM,GC,IGC,COORD_OLD,ICOORD)  

c-------------------------------------------------

          CALL PUT_ELEMENT (NUMS,INUMS,IEL, NUM, NOD_OLD,NDIME_OLD
     +      ,ITYPE  ,-1,IUSER1,IUSER2,IUSER3)    !- just kill

c------------ set up the indecs to loop ALL daughter elements ----------
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
              SAMP(1,J)      = WIDTHS(L(J),J)       !- base + offset
     +   + (WIDTHS(L(J)+1,J) - WIDTHS(L(J),J)) *(COORD_L(INODE,J)+1.)/2.
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
                COORD_G(J) = SAMP(1,J)   !- but widths ?? 
              ENDIF

C----------------------------------------------------------
            ENDDO   !- the dimensions of the new element
               
c----------------- build the nodes of the new element ------------------ 
            CALL FIND_OR_ADD_NODE (GC,IGC,NDIM,COORD_G,1,NN_NEW,JNODE)
            NUM(INODE) = JNODE
          ENDDO     !-- (loop nodes of the new element)

c--------------------- add-in the new element ---------------------------
          IF (N_SUB_EL.EQ.1) THEN
            IEL_NEW = IEL                !- so replace the old one
          ELSE
            NEL_NEW = NEL_NEW + 1        !- so build the new elements
            IEL_NEW = NEL_NEW            !- at the end of the list
          ENDIF

          CALL PUT_ELEMENT 
     +    (NUMS,INUMS,IEL_NEW, NUM, NOD_NEW,NDIME_NEW,ITYPE, IMAT,
     +     IUSER1,IUSER2,IUSER3)

C------------- end of this new element .. increment pointers -----------

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
     +  ILOOP,' :',ic,' subdivided into', N_SUB_EL,' each'
      ENDDO         !-- end of main loop

  999 CONTINUE

      NEL = NEL_NEW
      NN  = NN_NEW
      END            

C-----------------------------------------------------------------------
      SUBROUTINE MUNG_4NQ (IO,GC,IGC,NN,NDIM,NUMS,INUMS,NEL,P)
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
      SUBROUTINE R_NODMOV (IO,GC,IGC,NDIM,NN)
C
C     This moves nodes based on their 'old' and 'new' coords
C
      PARAMETER (MAX_NODES = 100)
      REAL GC(IGC,*), COLD(3), CNEW(3)
      INTEGER NUMLST(MAX_NODES) 

      WILD = 230964.
    1 DO K=1 , 9999
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
C     **** MUST CHANGE TO USE GET/PUT_ELEMENT ! ****
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*), P(*)

      DO KK = 1,999
    1   READ (IO,*,IOSTAT=IOS) IMAT_
        CALL IN_TEST (IO,IOS,*1,*999)
        NEL1 = NEL    
        IC = 0
        DO IEL=1,NEL
c         NOD  = NUMS (INUMS  ,IEL)
          IMAT = NUMS (INUMS-3,IEL)
          IF (IMAT.NE.IMAT_) THEN               !-- If a 'good' element ..
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
      SUBROUTINE DEL_ZERO_AREA (IO,GC,IGC,NN,NUMS,INUMS,NEL,P,NDIM)
C
C     this will delete any elements which have 'zero' area
C     **** MUST CHANGE TO USE GET/PUT_ELEMENT ! ****
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*), P(*)


C----------------------- Workspace arrays ------------------------------
      PARAMETER (M_NOD =  27    !-- max # nodes per element
     +          ,M_NODOF= 4     !-- max # freedoms per node
     +          ,ISMP  = 1      !-- max GP's per element
     +                      )
      PARAMETER ( IDER=M_NODOF, IJAC=4, ICOORD=M_NOD)

      REAL    FUN (M_NOD)          !-- shape funs
     +       ,DER (IDER  ,M_NOD)   !-- derivs of shape funs
     +       ,JAC (IJAC,IJAC)      !-- the 'Jacobian'
     +       ,SMP (ISMP,M_NODOF)   !-- Integration sampling points
     +       ,WTS (ISMP)           !-- Integration point 'weights'
      INTEGER NUM (M_NOD)          !-- element node numbers
      REAL  COORD(ICOORD,M_NODOF)

      NDIME = NDIM

      NEL1 = NEL    
      IC = 0
      DO IEL=1,NEL
c.. OK could use 'Danplot' area method
c... or I could use Shape funs.. hence.. get DET
        CALL GET_ELEMENT (NUMS,INUMS,IEL, NUM, NOD,NDIME,ITYPE
     +    ,IMAT,IUSER1,IUSER2,IUSER3)
        CALL GET_COORD (NUM,NOD,NDIM,GC,IGC,COORD,ICOORD)    !- coords

        CALL GET_ANY_GAUSS (SMP,ISMP,WTS,NDIME,1,NOD)
        CALL GSF  (NDIME,NOD,ITYPE, DER,IDER,FUN,SMP,ISMP,1)
        CALL MAT_MULT (DER,IDER,COORD,ICOORD,JAC,IJAC,NDIME,NOD,NDIME)
c       CALL INVERT_JAC (JAC,IJAC,DET,NDIME)
        DET=JAC(1,1)*JAC(2,2)-JAC(1,2)*JAC(2,1)

        IF (DET.GE. 0.0001) THEN               !-- If a 'good' element ..
          IC = IC + 1 
          DO J = 1,INUMS
            NUMS(J,IC)=NUMS(J,IEL)            !-- store it at IC
          ENDDO
        ENDIF
       ENDDO
       NEL = IC
       WRITE(*,'(A,I4,A,I4,A)')  '>> ',1,':',
     +               NEL1-NEL,' elements deleted'
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE CHANGE_MATS (IO,GC,IGC,NN,NUMS,INUMS,NEL,NDIM)
C
C     This will change any elements of a given material type
C     to another material type
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*), NUM(99)

      DO KK = 1,999
    1   READ (IO,*,IOSTAT=IOS) ITYPE1, ITYPE2
        CALL IN_TEST (IO,IOS,*1,*999)

        IC = 0
        DO IEL=1,NEL
          CALL GET_ELEMENT (NUMS,INUMS,IEL, NUM, 
     +          NOD,NDIME,ITYPE,IMAT,IUSER1,IUSER2,IUSER3)

          IF (IMAT.EQ.ITYPE1) THEN               !-- If a 'good' element ..
            IMAT = ITYPE2
            IC = IC + 1 
            CALL PUT_ELEMENT (NUMS,INUMS,IEL, NUM, 
     +          NOD,NDIME,ITYPE,IMAT,IUSER1,IUSER2,IUSER3)

          ENDIF
        ENDDO
        WRITE(*,'(A,I4,A,I4,A)')  '>> ',KK,':',
     +               IC,' elements changed'
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
      END
C-----------------------------------------------------------------------
      SUBROUTINE SORT_NODES (IO,GC,IGC,NN,NUMS,INUMS,NEL,P,NDIM)
C
C     This will sort the nodes wrt distance from a point
C
      REAL GC(IGC,*), CENTRE(5), DEPTHS(10000)
      INTEGER NUMS(INUMS,*), P(*)

    1 READ (IO,*,IOSTAT=IOS)  (CENTRE(J),J=1,NDIM)
      CALL IN_TEST (IO,IOS,*1,*999)
      DO I=1,NN
        D = 0.
        DO J=1,NDIM
          D = D + (GC(J,I) - CENTRE(J))**2
        ENDDO
        DEPTHS(I) = D
      ENDDO

      CALL DSORT@ (P,DEPTHS,NN)
      OPEN (97,STATUS='SCRATCH',FORM='UNFORMATTED')
      DO I=1,NN
        WRITE(97) P(I)
      ENDDO
      REWIND(97)
      DO I=1,NN
      READ(97) J
        P(J) = I
      ENDDO
      CLOSE(97)

      CALL UPDATE_GC   (GC,IGC,NDIM,NN, P)
      CALL UPDATE_NUMS (NUMS,INUMS,NEL, P)

C     WRITE(*,'(A,I4,A)')  '>> ',NN-NN_NEW,' nodes deleted'
C     NN = NN_NEW
  999 RETURN
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
C      *OBSOLETE as routines can multiply call FIIND_NODE
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
C     held in P
C     -- if P(I) = 0 then this node is removed
C       ---->  This WONT work in general as references are recursive
C
      REAL GC(IGC,*)
      INTEGER P(*)
      IC = 0
      open (97,status='scratch',form='unformatted')
c      DO INODE=1,NN
c        IF (P(INODE).GT.0) THEN
c          IC = IC + 1
c         GC(J,P(INODE)) = GC(J,INODE)
c          WRITE(97)  GC(1,P(INODE)), GC(2,P(INODE))
c          write(*,'(i8,$)') P(INODE)
c        ENDIF
c      ENDDO
c     NN = IC     !- added this 2-11-93

      DO I = 1,NN
        WRITE(97) (GC(J,I),J=1,NDIM)
      ENDDO
      REWIND(97)
      DO I=1,NN
        IF (P(I).NE.0) THEN
          IC = IC + 1
          READ(97) (GC(J,P(I)),J=1,NDIM)
        ELSE
          READ(97) (JUNK,J=1,NDIM)
        ENDIF
      ENDDO


c      REWIND(97)
c      DO I=1,IC
cc      DO I=1,NN
c        READ(97) GC(1,I),GC(2,I)
c      ENDDO
      CLOSE(97)
c     pause
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE UPDATE_NUMS (NUMS,INUMS,NEL, P)
C
C     This re-orders the node-steering in NUMS based on the NEW order
C     held in P  
C     if IOP=1 (default) then P() is the NEW number
C     if IOP=2  then P() is the position to put this number in (ie. inverse)
C
      INTEGER NUMS(INUMS,*), P(*), NUM(99), NUM2(99)
      IOP = 1
      DO IEL=1,NEL
        CALL GET_ELEMENT (NUMS,INUMS,IEL, NUM, 
     +          NOD,NDIME,ITYPE,IMAT,IUSER1,IUSER2,IUSER3)
        IF(IOP.EQ.1) THEN
          DO J=1,NOD
            NUM2(J) = P(NUM(J))       !-- if ZERO delete elem ?
          ENDDO
        ELSEIF(IOP.EQ.2) THEN
          DO J=1,NOD
            NUM2(J) = P(NUM(J))   
          ENDDO
        ELSE
          STOP '*OOPS* IN UPDATE_NUMS'
        ENDIF
        CALL PUT_ELEMENT (NUMS,INUMS,IEL, NUM2, 
     +          NOD,NDIME,ITYPE,IMAT,IUSER1,IUSER2,IUSER3)
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE SCALE_MESH (IO,GC,IGC,NN,NDIM)
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
      SUBROUTINE SHIFT_MESH (IO,GC,IGC,NN,NDIM)
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

      DO ILOOP = 1,9999

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



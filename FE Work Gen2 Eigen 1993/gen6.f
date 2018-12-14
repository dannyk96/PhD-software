      OPTIONS  (DREAL, INTL)
c     OPTIONS  (DREAL, INTL, FULLCHECK, UNDEF)

      PROGRAM DANS_GENERAL_FE_PROGRAM

      PARAMETER (
     +       IKV  =200 500    !-- Max. size of the stiffness matrix
     +      ,MNN  = 900      !-- Max. # of nodes
     +      ,MEL  = 500      !-- Max. # of elements

     +      ,MNOD = 20       !-- Max. # of nodes per element
     +      ,MDIM = 3     )  !-- Max 2d/3d 

C======================================================================
C            D A N ' S    2 D / 3 D  F E - P R O G R A M              !
C======================================================================
C                                                                     !
C                                                                     !
C                                                                     !
C======================================================================
C     REVISION HISTORY                                                !
C
C   8- 4-93   FORCES_TOT # total applied loads
C   7- 4-93  GEN6  : Built-up/ Excavation support :-)
C   1- 4-93  GEN5g : Von Mise if IOP*4 .eq. 2 (M-C is =1)             !
C  26- 3-93  GEN5f : plasticity to a subroutine : PCG support         !
C  25- 3-93  GEN5d : first 'working' copy                             !
C  11- 3-93  Major work: BUILD_KV to its own subroutine               !
C  16- 2-93  Gen5a.. back to work :-)                                 !
C  26- 1-93  GEN5 .. complete rework with keywords etc.               !
C  17-10-91  PARMETER's shuffled, FOS looping added, MUNG extended    !
C                                                                     !
C  May 1991  This is a re-hash of my code for FOS of 3d slopes        !
C----------------------------------------------------------------------
 
C     Max. element size parameters        Max. problem size params
C     ----------------------------        -------------------------
C        MN   = nodes per elem.           
C        MDOF = Dof. per mode             MEL = elements in mesh
C        MDF  = freedoms per elem         IKV = global stiffness matrix
C               
C        NS   = strains per GP            ILDS  = loaded nodes per mesh
C        IGRULE = GP's per elem           IQINC = load steps
C                                         IPRPTAB = material types
C-----------------------------------------------------------------------

       PARAMETER (
     +            MQINC  = 50          !-- Max. steps per loading
     +          ,IPRPTAB = 20          !-- Max. # of props per material
                                       
     +          ,MDF = MDIM            !-- max Dimensionality
     +       ,IGRULE = 27              !-- max # of GP's per element
     +       ,ILOADS = MNN*MDF* 7/10   !-- max #of freedoms
     +                               )
C--------------------- 'dependant' array sizes -------------------------
      PARAMETER ( 
     +                NS = MDF*2       !-- stress-terms per GP
     +             ,MSEL = IGRULE*NS)  !-- stress-terms per elem

c------------------- dummy sizes for subroutine calls ------------------
      PARAMETER (
     +               IGC = MDF         !-- width of GC (3 ?)
     +              ,INF = MDF         !-- node-freedom array
     +            ,INUMS = MNOD + 7    !-- width of NUMS (20+2 ?)
     +            ,IEVPT = MSEL        !-- size of EVPT / EVSTRS
     +           ,IKDIAG = ILOADS )    !-- length of KDIAG
C-----------------------------------------------------------------------

C ----------------- problem size depandant arrays --------------------
      REAL       KV (IKV)         !-- The Global Stiffness Matrix
     +          ,GC (IGC,MNN)     !-- the coordinates of all the nodes
     +      ,EVSTRS (IEVPT,MEL)   !-- all the stresses
     +        ,EVPT (IEVPT,MEL)   !-- all the strains
     +       ,LOADS (ILOADS)      !-- the applied loads  (at freedoms)
     +       ,DISPS (ILOADS)      !-- disps of this load step
     +   ,DISPS_TOT (MDF,MNN)     !-- the cumulative displacements
     +   ,DISPS_OLD (ILOADS)      !-- disps from the last iteration 
     +    ,BDYLDS_N (MDF,MNN)     !-- the iterated excess forces
     +    ,FORCES_N (MDF,MNN)     !-- applied nodal loads
     +  ,FORCES_TOT (MDF,MNN)     !-- total applied nodal loads 
c    +     ,REACT_N (MDF,MNN)     !-- resulting reaction forces (cf BDYLDS_N)


      INTEGER    NF (INF,MNN)     !-- the node freedom array
     +        ,NUMS (INUMS,MEL)   !-- the nodes attached to every element
     +       ,KDIAG (IKDIAG)      !-- pointers to each row of KV
     +    ,Q_NUMBER (MQINC)       !-- number of steps of this amount

      REAL   PRPTAB (IPRPTAB,10)  !-- the material properties 
     +   ,Q_PERCENT (MQINC)       !-- the amount of each load-step 

C ---------------------- standard FE arrays --------------------------


C --------------------------------------------------------------------
      INTEGER  U0                 !-- 'root' data file
     +        ,U1                 !-- 'current' data file 
     +        ,U2                 !--  Danplot file
     +        ,U3                 !--  ?
     +        ,IOPS(10)           !-- job-control options
     +        ,IDBG(10)           !-- codes for 'debug' outputs
C---------------------------------------------------------------------
      CHARACTER   FILE*80         !-- the data file name
     +         ,CMNAM@*80         !-- the command line function
     +          ,TIME@*8          !-- the time function
     +        ,KEYWORD*60         !-- the keyword token

      LOGICAL CONVERGED           !-- convergence flag

      DATA U0,U2,U3/40,12,13/   !- file numbers

C-----------------------------------------------------------------------
      DATA (IOPS(I),I=1,10) /
     +   0,  ! (1)
     +   1,  ! (2)  1 for plasticity (M_C) .. 0=elastic only
     +   0,  ! (3)  1 for Conjugate gradient method : 0 = SPARIN
     + 100,  ! (4)   Max Iterations for convergence
     +   3,  ! (5)   Convergence tolerance = 1/10**val
     +   2,  ! (6)   NGP code  NODOF**val  (0= off)
     +   1,  ! (7)  1= Mohr-Coulomb  ; 2= Von Mise 
     +   1,  ! (8)  1= loads are 'forces'; 2: loads are 'displacements'
     +   0,  ! (9) 
     +   0/  !(10) 

c.... also default # load-steps per load-case
c....      if a built-up model
c....      default Gauss-point rule (4 in 2d, 14? in 3D)     
c....           
c....           
C-----------------------------------------------------------------------
      DATA (IDBG(I),I=1,10) /
     +   1,  ! (1)     output every load step
     +   0,  ! (2)     output every iteration -
     +   0,  ! (3)                            -
     +   0,  ! (4)     node to output its disps every load step ?
     +   0,  ! (5)  
     +   0,  ! (6)
     +   0,  ! (7)
     +   0,  ! (8) 
     +   0,  ! (9) 
     +   0/  !(10) 


      DATA NN/0/, NEL/0/                   !- initially no elements
      DATA FOS_C/1./, FOS_PHI/1./          !- use the full c/phi
      DATA  Q_NUMBER (1)/1/
     +    , Q_PERCENT(1)/100./,IQ_TOT/1/   !- 100% of the load
      DATA BIG_SPR_KV / 1.E30/      !- arbitary 'Big spring'
C-----------------------------------------------------------------------
      WRITE (*,'(A)') '-------- Dan''s General FE Program -------'


c------------ get the input file name and generate output files --------
      FILE = CMNAM@()                  !-- the command-line
      IF (FILE.EQ.' ') CALL SELECT_FILE@('*.d*',FILE,*999)
      CALL UPCASE (FILE)
      OPEN (U0,FILE=FILE,STATUS='OLD') 
      OPEN (U2,FILE=FILE(1:INDEX(FILE,'.D'))//'RES')  !-- Output
      OPEN (U3,FILE=FILE(1:INDEX(FILE,'.D'))//'PL ')  !-- Danplot

      U1 = U0                         !-- start from the base file = U0
      IKEYWORD = 0
C---------------- loop back point for multiple load-steps --------------

 1000 CONTINUE
 

C----------------- Parse the data file (s) -----------------------------
        IKEYWORD = IKEYWORD + 1
        CALL GET_KEYWORD (U1,U0,KEYWORD)
        WRITE(*,'(I3,A,A)') IKEYWORD, ' ',KEYWORD

C-----------------------------------------------------------------------
      IF (KEYWORD.EQ.'*EOF') THEN        !-- end of program so exit
        GOTO 888                         
C-----------------------------------------------------------------------
      ELSEIF (KEYWORD.EQ.'*ANALYSE') THEN
        GOTO 1234              !--  run the rest of the program

C-----------------------------------------------------------------------
      ELSEIF (KEYWORD.EQ.'*DEBUG') THEN
      CALL R_OPTIONS(U1,IDBG)

      ELSEIF (KEYWORD.EQ.'*CONTROL') THEN
      CALL R_OPTIONS(U1,IOPS)

C----------------- Form and reduce the global stiffness matrix ---------
      ELSEIF (KEYWORD.EQ.'*FORM_KM') THEN

      CALL FORM_KV (GC,IGC,NN,NDIM,NUMS,INUMS,NEL,PRPTAB,IPRPTAB
     +         ,NF,INF,NODOF,N,ILOADS,KV,IKV,IR,KDIAG,IOPS(3)) 

      IF (IOPS(3).EQ.0) THEN          !-- SPARIN method
        IF (IOPS(8).EQ.2) THEN           !-- displaceent loading
          DO I=1,NN
            DO J=1,NODOF
              IFREE = NF(J,I)
              IF (IFREE.NE.0.AND.ABS(FORCES_N(J,I)).GT.1.E-8) THEN
                KV (KDIAG(IFREE)) = BIG_SPR_KV    !- big-spring
              ENDIF
            ENDDO
          ENDDO
        ENDIF
        CALL SPARIN (KV,N,KDIAG) !- else PCG
      ELSE
        IF (IOPS(8).EQ.2) 
     +  CALL MYERROR (2,'Program cant do DISP loads yet - sorry')
      ENDIF

C-----------------------------------------------------------------------
c.... this also 'triggers' the NULLing of the problem sized arrays
      ELSEIF (KEYWORD.EQ.'*DANPLOT_MESH') THEN
        CALL  R_DANPLOT (U1,GC,IGC,NDIM,NN,NUMS,INUMS,NEL)
        CALL WR_DANPLOT (U3,GC,IGC,NDIM,NN,NUMS,INUMS,NEL) !-for plotting

        NODOF = NDIM                 !-- default is the same as NDIM
        IF (NODOF.EQ.1) IH = 1
        IF (NODOF.EQ.2) IH = 4      !- include the z term
        IF (NODOF.EQ.3) IH = 6
        NGP = NDIM**IOPS(6)         !- ie. product gauss !
        CALL NULL (EVSTRS,IEVPT, IH*NGP,NEL)   !- total stresses
        CALL NULL (DISPS_TOT, MDF ,NODOF,NN)   !- total disps
        CALL NULL (FORCES_N,MDF ,NODOF,NN)     !- applied forces
        CALL NULL (FORCES_TOT,MDF ,NODOF,NN)   !- total applied forces

        CALL ON_FREEDOMS (NUMS,INUMS,NEL, NF,INF,NN,NODOF)

C-----------------------------------------------------------------------
      ELSEIF (KEYWORD.EQ.'*ADD/REMOVE_ELEMENTS') THEN
        CALL R_ADD_REM_ELEMS (U1,NUMS,INUMS,NEL)
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      ELSEIF (KEYWORD.EQ.'*MATERIAL_PROPERTIES') THEN
        CALL R_PRPTAB (U1,PRPTAB,IPRPTAB, NTYPES,NPROPS)

C-----------------------------------------------------------------------
      ELSEIF (KEYWORD.EQ.'*BC_BY_COORD') THEN      !- by coord (+wildcard)
        CALL R_BC_COORD (U1,GC,IGC,NN,NDIM,NF,INF,NODOF)

c      ELSEIF (KEYWORD.EQ.'*BC_BY_BOX') THEN        !- by box 
c        CALL R_BC_BOX (U1,GC,IGC,NN,NDIM,NF,INF,NODOF)

C---------------------- read loadings ----------------------------------

C-----------------------------------------------------------------------
      ELSEIF (KEYWORD.EQ.'*APPLY_CONSOL') THEN
        CALL APPLY_CONSOL (GC,IGC,NN,NDIM,NUMS,INUMS,NEL,
     +                    PRPTAB,IPRPTAB,EVSTRS,IEVPT,NODOF)

C-----------------------------------------------------------------------
      ELSEIF (KEYWORD.EQ.'*APPLY_GRAVITY') THEN
        CALL APPLY_GRAVITY (GC,IGC,NN,NDIM,NUMS,INUMS,NEL,
     +                    PRPTAB,IPRPTAB,FORCES_N,MDF,NODOF)
C-----------------------------------------------------------------------
      ELSEIF (KEYWORD.EQ.'*NODAL_LOADS') THEN
        IOPS(8) = 1
        CALL R_N_LOADS (U1,GC,IGC,NN,NDIM,FORCES_N,MDF,NODOF)

      ELSEIF (KEYWORD.EQ.'*NODAL_DISPLACEMENTS') THEN
        IOPS(8) = 2
        CALL R_N_LOADS (U1,GC,IGC,NN,NDIM,FORCES_N,MDF,NODOF)

C-----------------------------------------------------------------------
      ELSEIF (KEYWORD.EQ.'*FOS') THEN
        READ (U1,*) FOS_C,FOS_PHI

C-----------------------------------------------------------------------
c... kill the total disps, stresses etc.
      ELSEIF (KEYWORD.EQ.'*RESET_STATE') THEN
        CALL NULL (EVSTRS,IEVPT, IH*NGP,NEL)   !- total stresses
        CALL NULL (DISPS_TOT, MDF ,NODOF,NN)   !- total disps
C-----------------------------------------------------------------------
      ELSEIF (KEYWORD.EQ.'*LOAD_STEPS') THEN
         CALL R_LOAD_STEPS (U1,Q_NUMBER,Q_PERCENT,MQINC,IQ_TOT)

C-----------------------------------------------------------------------
      ELSE
        PRINT*,'UNKNOWN KEYWORD :'//KEYWORD
        CALL MYERROR (2,'Invalid keyword')
C-----------------------------------------------------------------------
      ENDIF

      GOTO 1000          !-- get the next keyword 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

 1234 CONTINUE                  !--- this is the 'rest of the program'


C----------------------- check array sizes -----------------------------
c.... maybe do this as we go along  ?  (and do a summary at the end ?)

      WRITE(*,'(79(''-''))')
      WRITE(*,'(A10,A10,I8,A10,I8)')   
     +   'nodes    ','NN  = ',NN      ,' MNN   =',MNN   
     +  ,'elements ','NEL = ',NEL     ,' MEL   =',MEL
     +  ,'freedoms ','N   = ',N       ,' ILOADS=',ILOADS
     +  ,'KV length','IR  = ',IR      ,' IKV   =',IKV
      WRITE(*,'(79(''-''))')

      IF   (      NN .GT. MNN
     + .OR. KDIAG(N) .GT. IKV 
     + .OR.      NEL .GT. MEL 
     + .OR.        N .GT. ILOADS ) 
     +  CALL MYERROR (2,'******* ARRAYS TOO SMALL ********')

C======================================================================
C======================= LOAD INCREMENT LOOP ==========================
C======================================================================

      CALL NULVEC (DISPS,N)    !- start at zero for first time PCG
      IY1  = 1      ! indeces into Q_NUMBER : step nuber
      IY2  = 0      !                       : so far of this step

      PTOT = 0      !-- total of the load applied so far

      N_LOAD_STEPS = 0   !-- work out the total number of load-steps
      DO I=1,IQ_TOT
        N_LOAD_STEPS = N_LOAD_STEPS + Q_NUMBER(I)
      ENDDO

      DO IY = 1,N_LOAD_STEPS

c...... first work out how much of the 'total' load to apply 
      IY2 = IY2 + 1
      IF (IY2.GT.Q_NUMBER(IY1)) THEN   !- go to the next stage
        IY1 = IY1 + 1
        IY2 = 1
      ENDIF

      STEP = Q_PERCENT (IY1) / 100.      !-- as a percentage
      PTOT = PTOT + STEP

c     PRINT*,'LOAD STEP',iy,' /',N_LOAD_STEPS,CHAR(27)//'[A'
c     PRINT*,'LOAD STEP',iy,' /',N_LOAD_STEPS,' STEP=',STEP

C======================= iteration loop ===============================

      CALL NULVEC (DISPS_OLD,N)                 !-- 'last' freedom disps
      CALL NULL   (EVPT  ,IEVPT, IH*NGP,NEL)    !- plastic strains
      CALL NULL   (BDYLDS_N,INF, NODOF,NN)      !- just for now ?

      ITER_MAX = IOPS(4)       !- eg 0.001 if val = 3
      TOL = 10.**(-IOPS(5))

      DO ITER = 1, ITER_MAX

C-------------------- form the vector of applied loads -----------------
      CALL NULVEC (LOADS,N)

      IF (IOPS(8).EQ.1) THEN
        FAC = 1.
      ELSEIF (IOPS(8).EQ.2) THEN
        FAC = BIG_SPR_KV
      ELSE
        CALL MYERROR (2,' Unknown method : Force/Displacemnt loading')
      ENDIF

      DO I=1,NN           !-- = number of 'loaded' freedoms
        DO J=1,NODOF
          IFREE = NF(J,I)
          IF (IFREE.GE.1) THEN    !--- only freedoms
            LOADS(IFREE) = 
     +      LOADS(IFREE) + FORCES_N(J,I) *STEP*FAC   !-- applied loads
          IF (ITER.GT.1)       !- no 'need' to do on the first pass :-)
     +      LOADS(IFREE) =
     +      LOADS(IFREE) + BDYLDS_N(J,I)       !-- XS bodyforces

           ENDIF
        ENDDO
      ENDDO

C======================= SOLVE THE EQUATIONS ===========================
      IF (IOPS(3).EQ.0) THEN
        CALL SPABAC (KV,LOADS,N,KDIAG)
        CALL VECCOP (LOADS,DISPS,N)
      ELSEIF (IOPS(3).EQ.1) THEN         !---- Conj. Grad. Method
        CALL SOLVE_CG (KV,LOADS,DISPS,N,KDIAG,NEL,ITERS_CG)  
      ENDIF
C=======================================================================

      IF (IOPS(2).EQ.0) THEN        !- for 'elastic' only
         CONVERGED = .TRUE.
      ELSE                          !- for 'elasto-plasticity'
        CALL CHECON3 (DISPS,DISPS_OLD,N,BIG,BIGGEST)
        CONVERGED = (BIGGEST.LT.TOL)    !-- convergence
        IF (IDBG(2).NE.0) WRITE(*,'(I5,A,F10.5,A,F10.5)') 
     +      ITER,' : Big=',BIG,' Biggest =',Biggest

        WRITE(*,'(A,I5,A,F10.5,A)')  
     +   'Iter=',ITER,' Biggest=',Biggest, CHAR(27)//'[A'
      ENDIF
      CALL VECCOP (DISPS,DISPS_OLD,N)       !- ready for the next pass

C----------------------------------------------------------------------

c... skip for now .. but REACTIONS are undefined :-(

      IOP_VM = IOPS(7)     !- for M-C / Von Mise
      CALL PLASTICITY 
     +   (GC,IGC,NN,NDIM,NUMS,INUMS,NEL,PRPTAB,IPRPTAB ,NF,INF,NODOF
     +  , DISPS,N, EVSTRS,IEVPT,EVPT,IEVPT,BDYLDS_N,INF,CONVERGED
     +  , IOP_VM,FOS_C,FOS_PHI)


C---------- end of GP and elem loops .. can do some monitoring --------

        IF (CONVERGED) THEN
          GOTO 1111
        ELSE
          IF (ITER.GE.ITER_MAX) CALL MYERROR(2, 'CONVERGENCE FAILED')
        ENDIF
      ENDDO    !-- loop the next iteration

C----------------------------------------------------------------------
C------------------- End of a load step so update arrays ---------------
C----------------------------------------------------------------------

 1111 CONTINUE   !- jump out of iteration loop

      DO I=1,NN             !-- first the total displacements
        DO J=1,NODOF
          IFREE = NF(J,I)
          IF (IFREE.NE.0)
     +    DISPS_TOT(J,I) = DISPS_TOT(J,I) + DISPS(IFREE)
        ENDDO
      ENDDO

      DO I=1,NN             !-- next the total applied forces
        DO J=1,NODOF
          FORCES_TOT(J,I) = FORCES_TOT(J,I) + FORCES_N(J,I) * STEP
        ENDDO
      ENDDO


C-------------------- Output various tables ----------------------------

C--- print the total force and the mean displacemnt of the loaded nodes
      S_D = 0.
      S_F = 0.
      IF (IY.EQ.1) S_F_OLD = 0.     !-- The 'last' reaction force
      DO J=1,NODOF
c     DO J=2,2          !- hack to vertical only for now ??
        IC = 0
        DO I=1,NN
          IF (ABS(FORCES_N(J,I)).GT.1.E-8) THEN     !-- if an applied force
            IC = IC + 1
            S_D = S_D + DISPS_TOT (J,I)
            S_F = S_F + BDYLDS_N (J,I)
          ENDIF
        ENDDO
        IF (IC.GT.0) THEN
        IF (IDBG(1).GT.0) THEN  
          IF (IY.EQ.1) WRITE(*,'(A3,3A12,A5,A9)') 
     +    'IY','Mean Disp', 'Total force','Force inc.','# Iters','TIme' 
          WRITE (*,'(I3,A,I5,A,3E12.4,A,A9)') 
     +    IY,' :', ITER,' :', S_D/IC, S_F, S_F - S_F_OLD ,' :',TIME@()
        ENDIF
       ENDIF
       ENDDO   !-- loop the freedoms
C      WRITE (U2,'( 2F20.10   )')
C     +    -S1/NL /.14,   -S2/67./(3.14159*.07**2) *4.
      S_F_OLD = S_F

C------------------- output displacements in DANPLOT style -------------

        WRITE (U3,'(A,I3,A,E12.5)')
     +    '#     Load-step No.',IY,' Total Load=',S_F
        WRITE(U3,*)NN
        DO I=1,NN
          WRITE (U3,'(I5,3E12.4)')  I,(DISPS_TOT(J,I),J=1,NODOF)
        ENDDO

C --- end of the load steeping loop.. print final state ? ------
      ENDDO    !--- loop the load-steps 

      GOTO 1000     !--- loop back for some more keywords

  888 STOP '> > > > > > > > Program completed OK < < < < < < < < <'
  999 STOP '************** PROGRAM ABORTED ****************'
      END

C-----------------------------------------------------------------------
C°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
C°°°°°°°°°°°°°°°°°° End of the Main Program °°°°°°°°°°°°°°°°°°°°°°°°°°°
C°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
C-----------------------------------------------------------------------













C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C    > > > > > > > >   Data Reading Routines   < < < < < < < 
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

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

C-----------------------------------------------------------------------
      SUBROUTINE R_N_LOADS (IO,GC,IGC,NN,NDIM,FORCES,IFORCES,NODOF)
C
C     This reads in the nodal loads based on their coordinates
C
      REAL GC(IGC,*), FORCES(IFORCES,*), POINT(5) ,data(5)
      DO I = 1,99999
    1   READ (IO,*,IOSTAT=IOS)  (POINT(J),J=1,NDIM),(data(j),j=1,nodof)
        CALL IN_TEST (IO,IOS,*1,*999)

C----------------now consult the table to find the nodes ---------------
        IFROM = 1
        DO IC = 1, 99999
          CALL FIND_NODE (GC,IGC,NDIM,POINT,IFROM,NN,INODE)
          IF (INODE.LT.1) GOTO 99
          DO J=1,NODOF
            FORCES (J,INODE) = FORCES(J,INODE) + DATA(J)
          ENDDO
          IFROM = INODE + 1
        ENDDO
  99    PRINT*, I,': Loads found=',IC-1
        ENDDO
 999  CONTINUE
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE R_BC_COORD (IO,GC,IGC,NN,NDIM,NF,INF,NODOF)
C
C     This reads in the nodal BOUNDARY CONDITIONS 
C     based on their coordinates
C
      REAL GC(IGC,*), POINT(5)
      INTEGER NF (INF,*), NF_DATA(5)
      DO I = 1,99999
    1   READ (IO,*,IOSTAT=IOS)  
     +  (POINT(J),J=1,NDIM), (NF_DATA(J),J=1,NODOF)
        CALL IN_TEST (IO,IOS,*1,*999)

C----------------now consult the table to find the nodes ---------------
        IC = 0
        IFROM = 1
        DO IC = 1, 999999
          CALL FIND_NODE (GC,IGC,NDIM,POINT,IFROM,NN,INODE)
          IF (INODE.LT.1) GOTO 99
          DO J=1,NODOF
            NF (J,INODE) = NF(J,INODE) * NF_DATA(J)  !- product-in!
          ENDDO
          IFROM = INODE + 1
        ENDDO
   99   PRINT*, I,': Boundary condition nodes found=',IC-1
        ENDDO
 999  CONTINUE 
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE R_LOAD_STEPS (IO, Q_NUMBER,Q_PERCENT,MQINC,IQ_TOT) 
C
C     This reads in the stages of appling the load
C     as the count and the percentage each time 
C        eq. 4 20. // 20 1.    == "4 of 20% then 20 of 1%"
C
      INTEGER Q_NUMBER (MQINC)
      REAL Q_PERCENT   (MQINC)
      DO I = 1,MQINC
    1   READ(IO,*,IOSTAT=IOS) Q_NUMBER(I),Q_PERCENT(I)
        CALL IN_TEST (IO,IOS,*1,*999)
      ENDDO
      CALL MYERROR (2, 'Too many load steps (IQ_TOT>MQINC)')
  999 IQ_TOT = I-1
      END
C-----------------------------------------------------------------------
      SUBROUTINE R_PRPTAB (IO,PRPTAB,IPRPTAB,NTYPES,NPROPS) 
C
C     This reads in the table of data (eg. material properties)
C
      REAL PRPTAB(IPRPTAB,*)

      NPROPS = 0
    1 READ(IO,*,IOSTAT=IOS) NTYPES     !- the # of data points per line
      CALL IN_TEST (IO,IOS,*1,*999)

      DO I =1,99999 
    2   READ(IO,*,IOSTAT=IOS) IMAT,(PRPTAB(J,IMAT),J=1,NTYPES)  !-  a line
        CALL IN_TEST (IO,IOS,*2,*999)
      ENDDO
  999 NPROPS = I-1
      PRINT*,NPROPS,' lines of material properties read'
      END


C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     > > > > > > >     Stiffness Matrix Routines   < < < < < < < < <
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      SUBROUTINE FORM_KV (GC,IGC,NN,NDIM,NUMS,INUMS,NEL,PRPTAB,IPRPTAB
     +         ,NF,INF,NODOF,N,ILOADS,KV,IKV,IR,KDIAG,IOP_PCG)
C
C     This forms the skyline stiffness matrix KV 
C     from the 'Database' of elements in NUMS and GC  
C     (Material props in PRPTAB)
C     IOP_PCG = 0 for SPARIN, =1 for conj. gradients
C

      REAL      GC (IGC,NN)     !-- the coordinates of all the nodes
     +     ,PRPTAB (IPRPTAB,*)  !-- the properties of each material type
     +         ,KV (IKV)        !-- the resulting stiffnes matrix

      INTEGER NUMS (INUMS,NEL)  !-- the topology of the elements
     +      ,KDIAG (*)          !-- pointers to the rows of KV
     +         ,NF (INF,NN)     !-- node freedom array
C----------------------- Workspace arrays ------------------------------
      PARAMETER (M_NOD =  20    !-- max # nodes per element
     +          ,M_NODOF= 2     !-- max # freedoms per node
     +          ,M_IH   = 4     !-- max # stresses per GP
     +          ,MDOF   = M_NOD*M_NODOF  !-- max # freedoms per element
     +          ,ISMP  = 27     !-- max GP's per element
     +                      )
      PARAMETER (IDER=M_NODOF, IDERIV=M_NODOF, IJAC=M_IH
     +    ,IBEE=M_IH, IKM=MDOF  ,ICOORD=M_NOD, IDEE=M_IH)

      REAL    FUN (M_NOD)          !-- shape funs
     +       ,DER (IDER  ,M_NOD)   !-- derivs of shape funs
     +     ,DERIV (IDERIV,M_NOD)   !-- derivs of shape funs - global
     +     ,COORD (ICOORD,M_NODOF) !-- nodal coords
     +       ,JAC (IJAC,IJAC)      !-- the 'Jacobian'
     +       ,DEE (IDEE,IDEE)      !-- the 'D' matrix
     +       ,BEE (IBEE,MDOF)      !-- the 'B' matrix
     +        ,KM (IKM,IKM)        !-- element stiffness matrix
     +       ,SMP (ISMP,M_NODOF)   !-- Integration sampling points
     +       ,WTS (ISMP)           !-- Integration point 'weights'

      INTEGER   G (MDOF)           !-- freedoms steering of an element
     +       ,NUM (M_NOD)          !-- node numbers of an element
C
C The concept of freedoms is ONLY relevant to the matrix solution
C so this subroutine must also return the KDIAG pointers to KV and
C resequence NF (?) to point to these freedoms  (so that applied FORCES
C may be packed into LOADS, and the resulting DISPLACEMENTS extracted)
C
C--------- an attempt at IMPLICIT NONE for tidyness -----------
      INTEGER NDIME,ITYPE,IMAT,IEL,IDOF,IGP,NGP

C----------------- first sort-out the freedoms ------------------------- 
      DO IEL=1,NEL
        CALL GET_ELEMENT (NUMS,INUMS,IEL, NUM, NOD,NDIME,ITYPE
     +      ,IMAT,IUSER1,IUSER2,IUSER3)
        IF (IMAT.GT.0) THEN                !- skip 'missing' materials
          DO I=1,NOD
            INODE = NUM(I)
            DO J=1,NODOF
              IF (NF(J,INODE).NE.0) NF(J,INODE) = 2
            ENDDO
          ENDDO
        ENDIF
      ENDDO   !- elements
C--------- switch off the '1's as they aren't here ---
      DO I=1,NN
        DO J=1,NODOF
          IF (NF(J,I).EQ.1) NF(J,I) = 0
        ENDDO
      ENDDO

      CALL SORTNF (NF,INF,NN,NODOF,N)  !-- check for a valid 'N' after ! 

      WRITE(*,'(A,I8,A,I8,A,F4.1,A)')
     +   '# of freedoms =',N,'/',ILOADS,' (',100.*N/REAL(ILOADS) ,'%)'
C---------------- form the KDIAG pointer list --------------------------
      DO I=1,N
        KDIAG(I) = 0 !-- 'null' first   (only needed for SPARIN method)
      ENDDO
      IBASE = 0
      DO IEL=1,NEL
        PRINT*,'Element',IEL,' /',NEL,CHAR(27)//'[A'

        CALL GET_ELEMENT (NUMS,INUMS,IEL, NUM, NOD,NDIME,ITYPE
     +      ,IMAT,IUSER1,IUSER2,IUSER3)
        IF (IMAT.GT.0) THEN                !- skip 'missing' materials
          IDOF = NOD * NODOF
          CALL NUM2G (NF,INF,NOD,NUM,NODOF,G) 
    
          IF (IOP_PCG.EQ.0) THEN
            CALL FORM_KDIAG (KDIAG,G,IDOF)       !- for SPARIN
          ELSEIF (IOP_PCG.EQ.1) THEN
C..... for PCG it is maybe neater to store G as we go along ? 
            CALL FORM_KDIAG_PCG (KDIAG,G,IDOF,IBASE,IKV)  !- store G's
          ENDIF
        ENDIF
      ENDDO

      IF (IOP_PCG.EQ.0) THEN
        CALL RESEQ_KDIAG (KDIAG,N  ,IKV)
        IR = KDIAG(N)
      ELSE
        IR = IBASE 
      ENDIF
C---------------- loop through the elements ---------------------------
      CALL NULVEC (KV,KDIAG(N))

      IF (NODOF.EQ.2) IH = 3    ! if axisym IH = 4 tho'
      IF (NODOF.EQ.3) IH = 6    ! be careful of FORM_BEE returning IH

      IBASE = 0     !- for PCG type KV 
      DO IEL=1,NEL
        PRINT*,'Element',IEL,' /',NEL,CHAR(27)//'[A'

        CALL GET_ELEMENT (NUMS,INUMS,IEL, NUM, NOD,NDIME,ITYPE
     +      ,IMAT,IUSER1,IUSER2,IUSER3)

        IF (IMAT.GT.0) THEN                !- skip 'missing' materials

        CALL GET_COORD (NUM,NOD,NDIM,GC,IGC,COORD,ICOORD)    !- coords
        CALL NUM2G     (NF,INF,NOD,NUM,NODOF,G)              !- freedoms

c.. OK I can now call a subroutine to just form KM
c... then do my own assembly ??
        IDOF = NOD * NODOF
        CALL NULL   (KM,IKM,IDOF,IDOF)

C----------- get material properties, hence form DEE -------------------
        E    = PRPTAB (1,IMAT)
        V    = PRPTAB (2,IMAT)
        BW   = PRPTAB (6,IMAT)               !    BW = a bit obsolete ?
        CALL FORM_DEE (DEE,IDEE,E,V,BW,NODOF,1)

c       NGP = 4**NDIME        !-- hack to full-integration for now
        NGP = 2**NDIME        !-- R.I.

        CALL GET_ANY_GAUSS (SMP,ISMP,WTS,NODOF,NGP,NOD)
        DO IGP=1,NGP
          CALL GSF  (NODOF,NOD,ITYPE, DER,IDER,FUN,SMP,ISMP,IGP)
          CALL MATMUL (DER,IDER,COORD,ICOORD,JAC,IJAC,NDIME,NOD,NDIME)
          CALL INVERT_JAC (JAC,IJAC,DET,NODOF)
c... need to modify DET if a triangle or 'axisym' 
          CALL MATMUL (JAC,IJAC,DER,IDER,DERIV,IDERIV,NDIM,NDIM,NOD)
          CALL FORM_BEE (BEE,IBEE,DERIV,IDERIV,NOD,NDIM)
c... need to extend BEE if axisymmetry
          WEIGHT=  DET * WTS(IGP)
          CALL FMBTDB (BEE,IBEE,DEE,IDEE,KM,IKM, WEIGHT, IH,IDOF)
        ENDDO   !-- GP's ---

c...... optional .. can also store as the full set of KMs for 
c...... PCG methods :-)
c... so store in KV too .. maybe use directly from KV via pointers :-)

        IF (IOP_PCG.EQ.0) THEN
          CALL FORM_SPARCE_KV (KM,IKM,IDOF,G,KDIAG, KV)
        ELSEIF (IOP_PCG.EQ.1) THEN
          CALL FORM_PCG_KV (KM,IKM,IDOF,IBASE, KV)
        ELSE
          CALL MYERROR (2,'Unknown method : SPARIN/PCG/?')
        ENDIF

        ENDIF   !-- skip of 'missing' elements
      ENDDO   !-- elements ---

      RETURN
      END

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE PLASTICITY 
     +   (GC,IGC,NN,NDIM,NUMS,INUMS,NEL,PRPTAB,IPRPTAB ,NF,INF
     +  ,NODOF, DISPS,N,EVSTRS,IEVSTRS,EVPT,IEVPT, BDYLDS_N,IBDYLDS_N
     +   ,CONVERGED, IOP_VM,FOS_C,FOS_PHI)
C
C     This finds the excess body-forces in a mesh
C      ie. loop elements & Gauss points:
C      disp inc.-> strain inc.-> stress inc. -> tot. stress
C           check stress , if too high : 
C    stress->(flow)-> flow strain-> flow stress-> flow forces->bdylds
C         if CONVERGED:
C         can update tot. stresses & do tot stress->react. forces instead
C         if IOP_VM =1 use Mohr-Coulomb; =2 use Von Mise

c------------------------ passed arrays --------------------------------
      REAL      GC (IGC,NN)     !- the coordinates of all the nodes
     +     ,PRPTAB (IPRPTAB,*)  !- the properties of each material type

     +      ,DISPS (N)           !- the freedom displacements
     +     ,EVSTRS (IEVSTRS,*)   !- total GP stresses
     +       ,EVPT (IEVPT,*)     !- total plastic strains
     +   ,BDYLDS_N (IBDYLDS_N,*) !- output nodal forces (bdylds/reactions)

      INTEGER NUMS (INUMS,NEL)  !- the topology of the elements
     +         ,NF (INF,NN)     !- node freedom array

      LOGICAL CONVERGED         !- flag if converged (so only update arrays)

C----------------------- Workspace arrays ------------------------------
      PARAMETER (M_NOD = 20     !- max # nodes per element
     +          ,M_NODOF= 3     !- max # freedoms per node
     +          ,M_IH   = 6     !- max # stresses per GP
     +          ,MDOF   = M_NOD*M_NODOF  !-- max # freedoms per element
     +          ,ISMP  = 27     !-- max GP's per element
     +                      )
      PARAMETER (IDER=M_NODOF, IDERIV=M_NODOF, IJAC=M_IH
     +    ,IBEE=M_IH, ICOORD=M_NOD, IDEE=M_IH,IFLOW=M_IH)

      REAL
     +        SMP (ISMP,M_NODOF)   !- Integration sampling points
     +       ,WTS (ISMP)           !- Integration point 'weights'
     +       ,FUN (M_NOD)          !- shape funs
     +       ,DER (IDER  ,M_NOD)   !- derivs of shape funs
     +     ,DERIV (IDERIV,M_NOD)   !- derivs of shape funs - global
     +     ,COORD (ICOORD,M_NODOF) !- nodal coords
     +       ,JAC (IJAC,IJAC)      !- the 'Jacobian'
     +       ,DEE (IDEE,IDEE)      !- the 'D' matrix
     +       ,BEE (IBEE,MDOF)      !- the 'B' matrix

     +       ,ELD (MDOF)           !- disps of a singel element
     +       ,EPS (M_IH)           !- tot strain inc.
     +    ,STRESS (M_IH)           !- STRESS increment (DEE*EPS)
     +     ,SIGMA (M_IH)           !- total STRESS   (TOTAL+STRESS)
     +      ,FLOW (IFLOW,IFLOW)    !- plastic 'flow' matrix (cf DEE)
     +       ,EVP (M_IH)           !- plastic 'flow' strain 
     + ,STRESS_FL (M_IH)           !- plastic 'flow' stress 

      INTEGER   G (MDOF)           !- freedoms steering of an element
     +       ,NUM (M_NOD)          !- node numbers of an element

C--------- an attempt at IMPLICIT NONE for tidyness -----------
      INTEGER NDIME,ITYPE,IMAT,IEL,IDOF,IGP,NGP

C---------------- loop through the elements ---------------------------
      IF (CONVERGED) CALL NULL (BDYLDS_N,IBDYLDS_N, NODOF,NN) ! hence reactions

      IF (NODOF.EQ.1) IH = 1
      IF (NODOF.EQ.2) IH = 4      !- include the z term
      IF (NODOF.EQ.3) IH = 6

      IBASE = 0     !- for PCG type KV 

      DO IEL=1,NEL
c       PRINT*,'Element',IEL,' /',NEL,CHAR(27)//'[A'

        CALL GET_ELEMENT (NUMS,INUMS,IEL, NUM, NOD,NDIME,ITYPE
     +      ,IMAT,IUSER1,IUSER2,IUSER3)

        IF (IMAT.GT.0) THEN                !- skip 'missing' materials
c.... need to check wether to skip this element or not ! (IMAT.GE.0) ?

        CALL GET_COORD (NUM,NOD,NDIM,GC,IGC,COORD,ICOORD)    !- coords
        CALL NUM2G     (NF,INF,NOD,NUM,NODOF,G)              !- freedoms

        IDOF = NOD * NODOF

C------------------ abstract the element's displacements ---------------
        IC = 0
        DO I=1,NOD
          DO J=1,NODOF
            IC=IC+1
            IFREE = NF(J,NUM(I))
            IF (IFREE.EQ.0) ELD(IC) = 0.
            IF (IFREE.NE.0) ELD(IC) = DISPS(IFREE)
          ENDDO
        ENDDO

C----------- get material properties, hence form DEE -------------------
        E    = PRPTAB (1,IMAT)
        V    = PRPTAB (2,IMAT)
        BW   = PRPTAB (6,IMAT)               !    BW = a bit obsolete ?
        CALL FORM_DEE (DEE,IDEE,E,V,BW,NODOF,1)

        C    = PRPTAB (3,IMAT)     !-- the plasticity parameters ---
        PHI  = PRPTAB (4,IMAT)
        PSI  = PRPTAB (5,IMAT)     !- FOS this too ?

        DTR  = 3.14159265 / 180.   !-- factor for 'degress to radians'
        C    = C /FOS_C  
        PHI  = ATAN (TAN(PHI*DTR) /FOS_PHI) /DTR


        NGP = 2**NDIME      !-- to match the outside world !
c.... this should really be the same as on the 'outside'!

        CALL NULL (BEE,IBEE, IH,MDOF)  !- only needed for the LAST row
        CALL GET_ANY_GAUSS (SMP,ISMP,WTS,NODOF,NGP,NOD)
        DO IGP=1,NGP
          CALL GSF    (NODOF,NOD,ITYPE, DER,IDER,FUN,SMP,ISMP,IGP)
          CALL MATMUL (DER,IDER,COORD,ICOORD,JAC,IJAC,NDIME,NOD,NDIME)
          CALL INVERT_JAC (JAC,IJAC,DET,NODOF)
c... need to modify DET if a triangle or 'axisym' 
          WEIGHT =  DET * WTS(IGP)

          CALL MATMUL (JAC,IJAC,DER,IDER,DERIV,IDERIV,NDIM,NDIM,NOD)
          CALL FORM_BEE (BEE,IBEE,DERIV,IDERIV,NOD,NDIM)
c...( need to extend BEE if axisymmetry )

          CALL MVMULT (BEE,IBEE,ELD,IH,IDOF,EPS)   !- tot. strain inc.

          DO J=1,IH                            !- remove the plastic strain
            IPT = (IGP-1)*IH + J
            EPS(J) = EPS(J) - EVPT(IPT,IEL)    !- hmm Ez is a dummy!
          ENDDO
          CALL MVMULT (DEE,IDEE,EPS,IH,IH,SIGMA)     !- stress inc.
          DO J=1,IH
            IPT = (IGP-1)*IH + J
            STRESS(J) = EVSTRS(IPT,IEL) + SIGMA(J)   !- total stress
          ENDDO

C--------- update the stresses and calculate the reactions -------------

      IF (CONVERGED) THEN
        DO K=1,IH                 
          IPT = (IGP-1)*IH + K
          EVSTRS(IPT,IEL) = STRESS(K)       !- update the total stress

          T = STRESS(K) * WEIGHT  
          IC = 0                      !-- hence the reaction forces
          DO I=1,NOD
            II = NUM(I)
            DO J=1,NODOF
              IC = IC + 1
              BDYLDS_N (J,II) = BDYLDS_N (J,II) + BEE(K,IC) * T
            ENDDO
          ENDDO
        ENDDO      !- the 4 (or 6) streses

C--------- Calc the Body-Loads if not converged and yielded ------------
      ELSE                          !--- get the total stress and check

        CALL INVARIENTS (STRESS,NODOF,SIGM,DSBAR,THETA)
        IF (IOP_VM.EQ.1) CALL MOCOUF (PHI,C, SIGM,DSBAR,THETA ,F) ! M-C 
        IF (IOP_VM.EQ.2) F=DSBAR-SQRT(3.)*C   !- Von Mise: plane-strain
        IF (IOP_VM.EQ.3) F=DSBAR-  2.    *C   !- Von Mise: triaxial
                             
        IF (F.GT.0.) THEN   

c------------------- yielded so get the XS stress ----------------------
c.. ie build the FLOW (DEE) matrix from PSI,XS-STRESS etc.
c.. hence plastic strains hence plastic stress (hence *BEE -> BDYLDS)

        CALL GET_FLOW  (DSBAR,THETA,PSI,E,V,PHI,F,STRESS,NODOF
     +                  ,FLOW,IFLOW,IOP_VM)
        CALL MVMULT  (FLOW,IFLOW,STRESS,IH,IH,EVP)       !- flow 'strain'
        CALL MVMULT  (DEE,IDEE  ,EVP  ,IH,IH,STRESS_FL)  !- flow 'stress'
        DO J=1,IH
          IPT = (IGP-1)*IH + J
          EVPT (IPT,IEL) = EVPT (IPT,IEL) + EVP(J)  !-- update pl. strains
        ENDDO

C------------------ calculate the excess Body-forces ------------------
        DO K=1,IH     !-- loop stress-components
          T = STRESS_FL(K) * WEIGHT     !- be careful about SIGMA 
          IC = 0
          DO I=1,NOD
            II = NUM(I)
            DO J=1,NODOF
              IC = IC + 1
              BDYLDS_N (J,II) = BDYLDS_N (J,II) + BEE (K,IC) * T
            ENDDO
          ENDDO
        ENDDO       !- K=1,IH
C----------------------------------------------------------------------
        ENDIF        !- only if a failed GP
      ENDIF       !- only if we were not 'converged'

        ENDDO   !-- GP's ---

        ENDIF   !-- only if element is 'present'  
      ENDDO   !-- elements ---

      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE CHECON3 (DISPS,DISPS_OLD,N,BIG,BIGGEST)
C
C     This find the largest displacement: BIG      and
C     the largest change in displacement: BIGGEST
C
      REAL DISPS(*),DISPS_OLD(*)
      BIG = 0.
      BIGGEST = 0.
      DO I=1,N
        BIG     = MAX (BIG     ,ABS(DISPS(I)))
        BIGGEST = MAX (BIGGEST, ABS(DISPS(I)-DISPS_OLD(I)))
      ENDDO
      BIGGEST = BIGGEST / BIG
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE APPLY_CONSOL (GC,IGC,NN,NDIM,NUMS,INUMS,NEL,
     +                    PRPTAB,IPRPTAB,EVSTRS,IEVSTRS,NODOF)
C
C     This forms the Initial Stresses in the Material due to both
C     self-weigth (GAMA) and overburden (CONS)
C
C     from the 'Database' of elements in NUMS and GC  
C     (Material props in PRPTAB)
C

      REAL     GC (IGC,NN)      !-- the coordinates of all the nodes
     +    ,PRPTAB (IPRPTAB,*)   !-- the properties of each material type
     +    ,EVSTRS (IEVSTRS,*)   !-- the Gauss-point stresses
      INTEGER NUMS (INUMS,NEL)  !-- the topology of the elements
C----------------------- Workspace arrays ------------------------------
      PARAMETER (M_NOD = 20     !-- max # nodes per element
     +          ,M_NODOF= 3     !-- max # freedoms per node
     +          ,M_NS   = 6     !-- max # stresses per GP
     +          ,ISMP  = 27 )   !-- max GP's per element
      PARAMETER (IDER=M_NODOF, ICOORD = M_NOD)

      REAL    FUN (M_NOD)          !-- shape funs
     +       ,DER (IDER,M_NOD)     !-- shape fun derivs
     +     ,COORD (ICOORD,M_NODOF) !-- nodal coords
     +       ,SMP (ISMP,M_NODOF)   !-- Integration sampling points
     +       ,WTS (ISMP)           !-- Integration point 'weights'
     +        ,XY (M_NODOF)        !-- coords of a Gauss point
     +    ,STRESS (M_NS)           !-- GP stress

      INTEGER NUM (M_NOD)          !-- node numbers of an element

      IF (NODOF.EQ.2) IH = 4    ! 'cos we have a 'Sz' term
      IF (NODOF.EQ.3) IH = 6    ! 

C--------------- loop through the elements ---------------------------
      DO IEL=1,NEL
        PRINT*,'Element',IEL,' /',NEL,CHAR(27)//'[A'

        CALL GET_ELEMENT (NUMS,INUMS,IEL, NUM, NOD,NDIME,ITYPE
     +      ,IMAT,IUSER1,IUSER2,IUSER3)

        IF (IMAT.GT.0) THEN                !- skip 'missing' materials
c...... hmm ..need to only do those elements that are 'new' !

        CALL GET_COORD (NUM,NOD,NDIM,GC,IGC,COORD,ICOORD)    !- coords

        EPK0 = 1
        GAMA = PRPTAB (8,IMAT)      !- hack out the Y-value
        CONS = PRPTAB (10,IMAT)      !- the consolidation (2D ?)

C------------------------ loop the G.P.'s -----------------------------
        NGP = 2**NDIME        !-- R.I.

        CALL GET_ANY_GAUSS (SMP,ISMP,WTS,NODOF,NGP,NOD)
        DO IGP=1,NGP
          CALL GSF    (NODOF,NOD,ITYPE, DER,IDER,FUN,SMP,ISMP,IGP)
          CALL MTVMULT (COORD,ICOORD,FUN,NDIM,NOD,XY)

          STRESS (1) = CONS - XY(2) * GAMA * EPK0
          STRESS (2) = CONS - XY(2) * GAMA    !- -ve 'cos downwards
          STRESS (3) = 0.
          STRESS (4) = CONS - XY(2) * GAMA * EPK0
          IF (NODOF.EQ.3) THEN
            STRESS(3) = STRESS(4)
            STRESS(4) = 0.
            STRESS(5) = 0.
            STRESS(6) = 0.
          ENDIF
          DO K=1,IH                 
            IPT = (IGP-1)*IH + K
            EVSTRS(IPT,IEL) = EVSTRS(IPT,IEL) + STRESS(K) 
          ENDDO
        ENDDO      !-- GP's
        ENDIF   !-- only if element is 'present'  
      ENDDO   !-- elements

      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE APPLY_GRAVITY (GC,IGC,NN,NDIM,NUMS,INUMS,NEL,
     +                    PRPTAB,IPRPTAB,GRAV_FORCES,IGF,NODOF)
C
C     This forms the gavity loading 
C     from the 'Database' of elements in NUMS and GC  
C     (Material props in PRPTAB)
C

      REAL     GC (IGC,NN)      !-- the coordinates of all the nodes
     +    ,PRPTAB (IPRPTAB,*)   !-- the properties of each material type
     +   ,GRAV_FORCES (IGF,*)   !-- the resulting nodal forces
      INTEGER NUMS (INUMS,NEL)  !-- the topology of the elements
C----------------------- Workspace arrays ------------------------------
      PARAMETER (M_NOD = 20     !-- max # nodes per element
     +          ,M_NODOF= 3     !-- max # freedoms per node
     +          ,M_NS   = 6     !-- max # stresses per GP
     +          ,ISMP  = 27 )   !-- max GP's per element
      PARAMETER (IDER=M_NODOF,IJAC=M_NS, ICOORD = M_NOD)

      REAL    FUN (M_NOD)          !-- shape funs
     +       ,DER (IDER,M_NOD)     !-- derivs of shape funs
     +       ,JAC (IJAC,M_NS)      !-- the 'jacobian'
     +     ,COORD (ICOORD,M_NODOF) !-- nodal coords
     +       ,SMP (ISMP,M_NODOF)   !-- Integration sampling points
     +       ,WTS (ISMP)           !-- Integration point 'weights'
     +       ,RHO (3)              !-- Gravity !
      INTEGER NUM (M_NOD)          !-- node numbers of an element
C--------- an attempt at IMPLICIT NONE for tidyness -----------
c     INTEGER NDIME,ITYPE,IMAT,IEL,IDOF,IGP,NGP

C--------------- loop through the elements ---------------------------
      DO IEL=1,NEL
        PRINT*,'Element',IEL,' /',NEL,CHAR(27)//'[A'

        CALL GET_ELEMENT (NUMS,INUMS,IEL, NUM, NOD,NDIME,ITYPE
     +      ,IMAT,IUSER1,IUSER2,IUSER3)

        IF (IMAT.GT.0) THEN                !- skip 'missing' materials
c...... hmm ..need to only do those elements that are 'new' !

        CALL GET_COORD (NUM,NOD,NDIM,GC,IGC,COORD,ICOORD)    !- coords

        DO J=1,NODOF
          RHO(J)  = PRPTAB (6+J,IMAT)      !- ie columns 7,8,9
          PRPTAB (6+J,IMAT) = 0.       !- to inhibit doing them twice!
        ENDDO

C------------------------ loop the G.P.'s -----------------------------
        NGP = 2**NDIME        !-- R.I.

        CALL GET_ANY_GAUSS (SMP,ISMP,WTS,NODOF,NGP,NOD)
        DO IGP=1,NGP
          CALL GSF    (NODOF,NOD,ITYPE, DER,IDER,FUN,SMP,ISMP,IGP)
          CALL MATMUL (DER,IDER,COORD,ICOORD,JAC,IJAC,NDIM,NOD,NDIM)
          CALL INVERT_JAC (JAC,IJAC,DET,NODOF)
c... need to modify DET if a triangle or 'axisym' 
          WEIGHT= DET * WTS(IGP)
          DO I=1,NOD
            DO J=1,NODOF
              GRAV_FORCES (J,NUM(I)) = 
     +        GRAV_FORCES (J,NUM(I)) + WEIGHT * FUN (I) * RHO(J)
            ENDDO
          ENDDO
        ENDDO      !-- GP's
        ENDIF   !-- only if element is 'present'  
      ENDDO   !-- elements

      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE ON_FREEDOMS (NUMS,INUMS,NEL, NF,INF,NN,NODOF)
C
C     This 'switches' on the freedoms at all 'present' nodes (IMAT>0)
C
      INTEGER NUMS (INUMS,NEL)    !-- eg. as max # elements
     +         ,NF (INF,NN)       !-- DITTO
     +        ,NUM (99)           !-- an elements' node numbers'
C-------------------- first switch OFF all freedoms --------------------
      DO I=1,NN
        DO J=1,NODOF
          NF(J,I) = 0
        ENDDO
      ENDDO

c------------ switch ON only the 'present' elements' freedoms ----------
      DO IEL=1,NEL
        CALL GET_ELEMENT (NUMS,INUMS,IEL, NUM, NOD,NDIME,ITYPE
     +      ,IMAT,IUSER1,IUSER2,IUSER3)
        IF (IMAT.GT.0) THEN    !-- (-ve) materials are 'absent'
          DO I=1,NOD
            DO J=1,NODOF
              NF(J,NUM(I)) = 1
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      END
c-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C      > > > > >    Numerical Integration Subroutines    < < < < <
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      SUBROUTINE GET_ANY_GAUSS (SMP,ISMP,WTS,NDIME,NGP,NOD)
C
C     a header routine to GAUSS which returns the quadrature
C     points as a single list (no implied looping) ie. in the same
C     form as NUMINT, NUMIN3 and QTURE
C     
C                                            Dan Kidger  April '91
C
C     NOD = the 'target' element to find a rule for :-)
C

C... should also return a 'Modifier' if the area of the domain is
C... less than the full integral (eq. triangles are only 1/2, tets 1/6)
C... (cf. Axisymetry where DET is further modified by the 'circumference'

      REAL SMP(ISMP,*), WTS(*)

      IF (NDIME.EQ.2.AND.                  !-- triangles
     +  (NOD.EQ.3.OR.NOD.EQ.6.OR.NOD.EQ.10.OR.NOD.EQ.15)) THEN
         CALL NUMINT (SMP,ISMP,WTS,NGP)

      ELSEIF (NDIME.EQ.3.AND.              !-- tetrahedra
     +  (NOD.EQ.4)) THEN
        CALL NUMIN3 (SMP,ISMP,WTS,NGP)

      ELSEIF (NDIME.EQ.3.AND.              !-- Irons rules 
     + (NGP.EQ.6.OR.NGP.EQ.13.OR.NGP.EQ.14.OR.NGP.EQ.15)) THEN
        CALL QTURE(SMP,ISMP,WTS,NGP)     

      ELSE                                 !-- product-Gauss
        CALL GET_PROD_GAUSS (SMP,ISMP,WTS,NDIME,NGP)
      ENDIF
      END
C-----------------------------------------------------------------------
      SUBROUTINE GET_PROD_GAUSS (SMP,ISMP,WTS,NDIME,NGP)
C
C     This recursively call GAUSS to build a product-Gauss rule
C
      REAL     SAMP (8,2)       !- one-dimensional rules
     +         ,SMP (ISMP,*)    !- the returned product rule
     +         ,WTS (*)         !- the returned product weights  
      INTEGER POINT (5)         !- 'toggles'

      NGPI = NINT (REAL(NGP)**(1./REAL(NDIME)) )

      IF (NGPI**NDIME.NE.NGP) THEN
        PRINT*,'NGP=',NGP,' is NOT a product Gauss rule !'
        CALL MYERROR (2,'in GET_PROD_GAUSS')
      ENDIF
      DO I=1,NDIME
        POINT(I)=1
      ENDDO
      CALL GAUSS (SAMP,8,NGPI)
      DO IC=1,NGP   !---- loop the points in the 'product' table

        WTS(IC)=1.
        DO I=1,NDIME      !--------- build this row -------------
          SMP(IC,I) = SAMP(POINT(I),1)
          WTS(IC)   = SAMP(POINT(I),2) * WTS(IC)
        ENDDO

        I=NDIME
    1   POINT(I) = POINT(I) + 1       !-- 'toggle' up 
        IF (POINT(I).GT.NGPI)THEN
         POINT(I)=1
         I = I - 1
         IF (I.GT.0) GOTO 1
       ENDIF
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE QTURE (SAMPL,ISAMPL,WT,NQP)
C
C     Weights and sampling points for Irons' 3D Integration rules
C
C     .... 13 point rule added 3/5/90    (so now 6,13,14,15 rules)
C
      REAL SAMPL(ISAMPL,*),WT(*)
      CALL NULL(SAMPL,ISAMPL,NQP,3)

      IF (NQP.EQ.6) THEN
        DO I=1,NQP
          WT(I)=4./3.
        ENDDO
        SAMPL(1,1)=-1.
        SAMPL(2,1)=1.
        SAMPL(3,2)=-1.
        SAMPL(4,2)=1.
        SAMPL(5,3)=-1.
        SAMPL(6,3)=1.

      ELSEIF (NQP.EQ.13) THEN
        B=-0.49584802
        C= 0.88030430
        D= 0.79562143
        E= 0.025293237
        DO I=1,6
          WT( I)= 0.54498736
        ENDDO
        DO I=7,12
          WT( I)= 0.507644216
        ENDDO
        WT(13)= 1.68421056
        DO I=1,6
          DO J=1,3
            SAMPL (I,J)=B
            SAMPL (I,MOD(I-1,3)+1) =C
          ENDDO
        ENDDO
        DO I=7,12
          DO J=1,3
            SAMPL(I,J)= D
            SAMPL(I,MOD(I-1,3)+1) =E
          ENDDO
        ENDDO
        DO J=1,3
          SAMPL(13,J)=0.
          DO I=4,9
            SAMPL(I,J)=-SAMPL(I,J)
          ENDDO
        ENDDO

      ELSEIF (NQP.EQ.14) THEN
        B=0.795822426
        C=0.758786911
        DO I=1,6
          WT(I)=0.886426593
        ENDDO
        DO I=7,NQP
          WT(I)=0.335180055
        ENDDO
        SAMPL(1,1) = -B
        SAMPL(2,1) =  B
        SAMPL(3,2) = -B
        SAMPL(4,2) =  B
        SAMPL(5,3) = -B
        SAMPL(6,3) =  B
        DO I=7,NQP
          DO J=1,3
            SAMPL(I,J) = C      !-- default value
          ENDDO
        ENDDO
        SAMPL(7,1)  = -C
        SAMPL(7,2)  = -C
        SAMPL(7,3)  = -C
        SAMPL(8,2)  = -C
        SAMPL(8,3)  = -C
        SAMPL(9,1)  = -C
        SAMPL(9,3)  = -C
        SAMPL(10,3) = -C
        SAMPL(11,1) = -C
        SAMPL(11,2) = -C
        SAMPL(12,2) = -C
        SAMPL(13,1) = -C
      ELSEIF (NQP.EQ.15) THEN
        B=1.
        C=0.674199862
        WT(1)=1.564444444
        DO I=2,7
          WT(I)=0.355555556
        ENDDO
        DO I=8,15
          WT(I)=0.537777778
        ENDDO
        SAMPL(2,1)=-B
        SAMPL(3,1)=B
        SAMPL(4,2)=-B
        SAMPL(5,2)=B
        SAMPL(6,3)=-B
        SAMPL(7,3)=B
        DO I=8,NQP
          DO J=1,3
            SAMPL(I,J)=C
          ENDDO
        ENDDO
        SAMPL(8,1) = -C
        SAMPL(8,2) = -C
        SAMPL(8,3) = -C
        SAMPL(9,2) = -C
        SAMPL(9,3) = -C
        SAMPL(10,1)= -C
        SAMPL(10,3)= -C
        SAMPL(11,3)= -C
        SAMPL(12,1)= -C
        SAMPL(12,2)= -C
        SAMPL(13,2)= -C
        SAMPL(14,1)= -C
      ELSE
        PRINT*,'Error in QTURE .. unknown integreation rule'
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE NUMINT (S,IS,WT,NGP)
C
C      This forms the sampling points and
C      weights for integration over a triangular area
C

c..  IN Fortran 90 I can use a set of 'internal' subroutines
c..  Each sets the constants as 'Parameters' 
c..  then uses simple 'array constructors' top build SMP and WTS
c..  eg.    REAL,PARAM, A / 0.3333/
c..          S= (/A,B,A,A,B,A, C,C,C,D,D,D/)      ok ?

      REAL S(IS,*),WT(*)
      IF (NGP.EQ.1) THEN
        S(1,1) = 1./3.
        S(1,2) = 1./3.
        WT(1) = 1.
      ELSEIF (NGP.EQ.3) THEN
        S(1,1) = .5
        S(1,2) = .5
        S(2,1) = .5
        S(2,2) = 0.
        S(3,1) = 0.
        S(3,2) = .5
        WT(1) = 1./3.
        WT(2) = WT(1)
        WT(3) = WT(1)
      ELSEIF (NGP.EQ.4) THEN
        S(1,1) = 1./3.
        S(1,2) = 1./3.
        S(2,1) = .6
        S(2,2) = .2
        S(3,1) = .2
        S(3,2) = .6
        S(4,1) = .2
        S(4,2) = .2
        WT(1) = -9./16.
        WT(2) = 25./48.
        WT(3) = WT(2)
        WT(4) = WT(2)
      ELSEIF (NGP.EQ.6) THEN
        S(1,1) = .816847572980459
        S(1,2) = .091576213509771
        S(2,1) = S(1,2)
        S(2,2) = S(1,1)
        S(3,1) = S(1,2)
        S(3,2) = S(1,2)
        S(4,1) = .108103018168070
        S(4,2) = .445948490915965
        S(5,1) = S(4,2)
        S(5,2) = S(4,1)
        S(6,1) = S(4,2)
        S(6,2) = S(4,2)
        WT(1) = .109951743655322
        WT(2) = WT(1)
        WT(3) = WT(1)
        WT(4) = .223381589678011
        WT(5) = WT(4)
        WT(6) = WT(4)
      ELSEIF (NGP.EQ.7) THEN
        S(1,1) = 1./3.
        S(1,2) = 1./3.
        S(2,1) = .797426985353087
        S(2,2) = .101286507323456
        S(3,1) = S(2,2)
        S(3,2) = S(2,1)
        S(4,1) = S(2,2)
        S(4,2) = S(2,2)
        S(5,1) = .470142064105115
        S(5,2) = .059715871789770
        S(6,1) = S(5,2)
        S(6,2) = S(5,1)
        S(7,1) = S(5,1)
        S(7,2) = S(5,1)
        WT(1) = .225
        WT(2) = .125939180544827
        WT(3) = WT(2)
        WT(4) = WT(2)
        WT(5) = .132394152788506
        WT(6) = WT(5)
        WT(7) = WT(5)
      ELSEIF (NGP.EQ.12) THEN
        S(1,1) = .873821971016996
        S(1,2) = .063089014491502
        S(2,1) = S(1,2)
        S(2,2) = S(1,1)
        S(3,1) = S(1,2)
        S(3,2) = S(1,2)
        S(4,1) = .501426509658179
        S(4,2) = .249286745170910
        S(5,1) = S(4,2)
        S(5,2) = S(4,1)
        S(6,1) = S(4,2)
        S(6,2) = S(4,2)
        S(7,1) = .636502499121399
        S(7,2) = .310352451033785
        S(8,1) = S(7,1)
        S(8,2) = .053145049844816
        S(9,1) = S(7,2)
        S(9,2) = S(7,1)
        S(10,1) = S(7,2)
        S(10,2) = S(8,2)
        S(11,1) = S(8,2)
        S(11,2) = S(7,1)
        S(12,1) = S(8,2)
        S(12,2) = S(7,2)
        WT(1) = .050844906370207
        WT(2) = WT(1)
        WT(3) = WT(1)
        WT(4) = .116786275726379
        WT(5) = WT(4)
        WT(6) = WT(4)
        WT(7) = .082851075618374
        WT(8) = WT(7)
        WT(9) = WT(7)
        WT(10) = WT(7)
        WT(11) = WT(7)
        WT(12) = WT(7)
      ELSEIF (NGP.EQ.16) THEN
        S(1,1) = 1./3.
        S(1,2) = 1./3.
        S(2,1) = .658861384496478
        S(2,2) = .170569307751761
        S(3,1) = S(2,2)
        S(3,2) = S(2,1)
        S(4,1) = S(2,2)
        S(4,2) = S(2,2)
        S(5,1) = .898905543365938
        S(5,2) = .050547228317031
        S(6,1) = S(5,2)
        S(6,2) = S(5,1)
        S(7,1) = S(5,2)
        S(7,2) = S(5,2)
        S(8,1) = .081414823414554
        S(8,2) = .459292588292723
        S(9,1) = S(8,2)
        S(9,2) = S(8,1)
        S(10,1) = S(8,2)
        S(10,2) = S(8,2)
        S(11,1) = .008394777409958
        S(11,2) = .263112829634638
        S(12,1) = S(11,1)
        S(12,2) = .728492392955404
        S(13,1) = S(11,2)
        S(13,2) = S(11,1)
        S(14,1) = S(11,2)
        S(14,2) = S(12,2)
        S(15,1) = S(12,2)
        S(15,2) = S(11,1)
        S(16,1) = S(12,2)
        S(16,2) = S(11,2)
        WT(1) = .144315607677787
        WT(2) = .103217370534718
        WT(3) = WT(2)
        WT(4) = WT(2)
        WT(5) = .032458497623198
        WT(6) = WT(5)
        WT(7) = WT(5)
        WT(8) = .095091634267284
        WT(9) = WT(8)
        WT(10) = WT(8)
        WT(11) = .027230314174435
        WT(12) = WT(11)
        WT(13) = WT(11)
        WT(14) = WT(11)
        WT(15) = WT(11)
        WT(16) = WT(11)
      ELSE
        PRINT*,'Error in NUMINT.. unknown integreation rule, NGP=',NGP
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE NUMIN3 (SAMP,ISAMP,WT,NGP)
C
C      This forms the sampling points and
C      weights for integration over a tetrahedron
C
      REAL SAMP(ISAMP,*),WT(*)
      IF (NGP.EQ.1) THEN
        SAMP(1,1) = .25
        SAMP(1,2) = .25
        SAMP(1,3) = .25
        WT(1) = 1.
      ELSEIF (NGP.EQ.4) THEN
        SAMP(1,1) = .58541020
        SAMP(1,2) = .13819660
        SAMP(1,3) = SAMP(1,2)
        SAMP(2,2) = SAMP(1,1)
        SAMP(2,3) = SAMP(1,2)
        SAMP(2,1) = SAMP(1,2)
        SAMP(3,3) = SAMP(1,1)
        SAMP(3,1) = SAMP(1,2)
        SAMP(3,2) = SAMP(1,2)
        SAMP(4,1) = SAMP(1,2)
        SAMP(4,2) = SAMP(1,2)
        SAMP(4,3) = SAMP(1,2)
        WT(1) = .25
        WT(2) = .25
        WT(3) = .25
        WT(4) = .25
        
      ELSEIF (NGP.EQ.5) THEN
        SAMP(1,1) = .25
        SAMP(1,2) = .25
        SAMP(1,3) = .25
        SAMP(2,1) = .5
        SAMP(2,2) = 1./6.
        SAMP(2,3) = SAMP(2,2)
        SAMP(3,2) = .5
        SAMP(3,3) = 1./6.
        SAMP(3,1) = SAMP(3,3)
        SAMP(4,3) = .5
        SAMP(4,1) = 1./6.
        SAMP(4,2) = SAMP(4,1)
        SAMP(5,1) = 1./6.
        SAMP(5,2) = SAMP(5,1)
        SAMP(5,3) = SAMP(5,1)
        WT(1) = -.8
        WT(2) = 9./20.
        WT(3) = WT(2)
        WT(4) = WT(2)
        WT(5) = WT(2)
      ELSE
        PRINT*,' Unknown integration rule in NUMIN3'
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GAUSS(SAMP,ISAMP,NGP)
C
C     This provides the weights and sampling points
C     for Gauss-Legendre quadrature   (ie. 1d hence quads & bricks)
C
      REAL SAMP(ISAMP,*)
      GO TO(1,2,3,4,5,6,7),NGP
    1 SAMP(1,1) = 0.
      SAMP(1,2) = 2.
      GOTO 100
    2 SAMP(1,1) = 1./SQRT(3.)
      SAMP(2,1) = -SAMP(1,1)
      SAMP(1,2) = 1.
      SAMP(2,2) = 1.
      GO TO 100
    3 SAMP(1,1) = .2*SQRT(15.)
      SAMP(2,1) = .0
      SAMP(3,1) = -SAMP(1,1)
      SAMP(1,2) = 5./9.
      SAMP(2,2) = 8./9.
      SAMP(3,2) = SAMP(1,2)
      GO TO 100
    4 SAMP(1,1) = .861136311594053
      SAMP(2,1) = .339981043584856
      SAMP(3,1) = -SAMP(2,1)
      SAMP(4,1) = -SAMP(1,1)
      SAMP(1,2) = .347854845137454
      SAMP(2,2) = .652145154862546
      SAMP(3,2) = SAMP(2,2)
      SAMP(4,2) = SAMP(1,2)
      GO TO 100
    5 SAMP(1,1) = .906179845938664
      SAMP(2,1) = .538469310105683
      SAMP(3,1) = .0
      SAMP(4,1) = -SAMP(2,1)
      SAMP(5,1) = -SAMP(1,1)
      SAMP(1,2) = .236926885056189
      SAMP(2,2) = .478628670499366
      SAMP(3,2) = .568888888888889
      SAMP(4,2) = SAMP(2,2)
      SAMP(5,2) = SAMP(1,2)
      GO TO 100
    6 SAMP(1,1) = .932469514203152
      SAMP(2,1) = .661209386466265
      SAMP(3,1) = .238619186083197
      SAMP(4,1) = -SAMP(3,1)
      SAMP(5,1) = -SAMP(2,1)
      SAMP(6,1) = -SAMP(1,1)
      SAMP(1,2) = .171324492379170
      SAMP(2,2) = .360761573048139
      SAMP(3,2) = .467913934572691
      SAMP(4,2) = SAMP(3,2)
      SAMP(5,2) = SAMP(2,2)
      SAMP(6,2) = SAMP(1,2)
      GO TO 100
    7 SAMP(1,1) = .949107912342759
      SAMP(2,1) = .741531185599394
      SAMP(3,1) = .405845151377397
      SAMP(4,1) = .0
      SAMP(5,1) = -SAMP(3,1)
      SAMP(6,1) = -SAMP(2,1)
      SAMP(7,1) = -SAMP(1,1)
      SAMP(1,2) = .129484966168870
      SAMP(2,2) = .279705391489277
      SAMP(3,2) = .381830050505119
      SAMP(4,2) = .417959183673469
      SAMP(5,2) = SAMP(3,2)
      SAMP(6,2) = SAMP(2,2)
      SAMP(7,2) = SAMP(1,2)
  100 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C        > > > > >      DEE and BEE  Subroutines    < < < < <
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
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

      CALL NULL (DEE,IDEE,IDEE,IDEE)

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
      SUBROUTINE FORM_BEE (BEE,IBEE,DERIV,IDERIV,NOD,NODOF)
C
C     This forms the strain-displacement matrix for 1D,2D & 3D etc.
C
      REAL BEE(IBEE,*), DERIV(IDERIV,*)

      IH = NODOF*(NODOF+1)/2                 != 1 (1D), 4 (2D), 6 (3D)
      CALL NULL (BEE,IBEE,IH ,NOD*NODOF)

      IB = 0
      DO I=1,NOD
        IV = NODOF 
        DO J=1,NODOF
          BEE (J,IB + J) = DERIV(J,I)        !- direct terms
          DO K=J+1,NODOF
            IV = IV + 1
            BEE(IV,IB + J) = DERIV(K,I)      !- shear terms
            BEE(IV,IB + K) = DERIV(J,I)
          ENDDO
        ENDDO
        IB = IB + NODOF    !- update the 'base' pointer
      ENDDO
      END            
C-----------------------------------------------------------------------
      SUBROUTINE FORM_VOL (DERIV,IDERIV,VOL,NOD,NODOF)
C
C     This forms a vector containing the derivatives of the shape funs ?
C        ie. the sum of the first NODOF rows of BEE
C     .. is it better to do this straight from DERIV or BEE ??
C
      REAL DERIV(IDERIV,*),VOL(*)
      IB = 0
      DO I=1,NOD
        DO J=1,NODOF
          VOL(IB) = DERIV (J,I)
          IB = IB + 1
        ENDDO
      ENDDO
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C      > > > > >      Basic Matrix Algebra Subroutines    < < < < <
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      SUBROUTINE NULL (A,IA,M,N)
C
C     This nulls a 2-d array
C
      REAL A(IA,*)
      DO J=1,N
        DO I=1,M              !- note: reversed loops
          A(I,J)=0.0
        ENDDO
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE NULVEC (VEC,N)
C
C     This nulls a column vector
C
      REAL VEC(*)
      DO I=1,N
        VEC(I) = 0.
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE MATCOP (A,IA,B,IB,M,N)
C
C     This copies matrix A to matrix B
C
      REAL A(IA,*), B(IB,*)
      DO J=1,N
        DO I=1,M                    ! note : reversed order
          B(I,J) = A(I,J)
        ENDDO
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE MSMULT (A,IA,C,M,N)
C
C     This multiplies a matrix by a scalar
C
      REAL A(IA,*)
      DO J=1,N
        DO I=1,M                !- note: reversed loops
          A(I,J) = A(I,J)*C
        ENDDO
      ENDDO
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
      SUBROUTINE MATMUL (A,IA,B,IB,C,IC,L,M,N)
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
C-----------------------------------------------------------------------
      SUBROUTINE VECCOP (A,B,N)
C
C     copies vector A into vector B
C
      REAL A(*), B(*)
      DO I = 1,N
        B(I) = A(I)
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE VECADD (A,B,C,N)
C
C     adds vectors    A + B = C 
C
      REAL A(*),B(*),C(*)
      DO I = 1,N
        C(I) = A(I) + B(I)
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE VECSUB (A,B,C,N)
C
C     vector subtract   :   C := A - B
C
      REAL A(*), B(*), C(*)
      DO I = 1,N
        C(I) = A(I) - B(I)
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE VSMULT (V,SCAL,N)
C
C     multiply a vector by a scalar
C
      REAL V(*)
      DO I = 1,N
        V(I) = V(I) * SCAL
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE VDOTV (V1,V2,DOTPR,N)
C
C     dot product  :  V1 * V2
C
      REAL V1(*),V2(*)
      DOTPR = 0.0
      DO I = 1,N
        DOTPR = DOTPR + V1(I)*V2(I)
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE VVMULT (V1,V2,PROD,IPROD,M,N)
C
C     This forms a vector product
C
      REAL V1(*),V2(*),PROD(IPROD,*)
      DO J=1,N
        DO I=1,M
          PROD(I,J) = V1(I)*V2(J)
        ENDDO
      ENDDO
      END

C-----------------------------------------------------------------------
C---------------------- PCG solving routines ---------------------------
C-----------------------------------------------------------------------
      SUBROUTINE FORM_PCG_KV (KM,IKM,IDOF,IC, KV)
C
C     This stores complete KM's in KV using IC as a pointer
C
      REAL KV(*), KM(IKM,IKM)

      DO J=1,IDOF
        DO I=1,IDOF
          IC = IC + 1
          KV(IC) = KM(I,J)    !- be careful about the order (tho' symetric)
        ENDDO
      ENDDO

      END
C-----------------------------------------------------------------------
      SUBROUTINE FORM_KDIAG_PCG (KDIAG,G,IDOF,IBASE,IKV) 
C
C     This stores IDOF and the G vector in KDIAG -so PCG loop is 'tight'
C
      INTEGER KDIAG (*), G(*)

      IF (IBASE+IDOF.GE.IKV) CALL MYERROR (2,'No more room in KDIAG')

      IBASE = IBASE + 1
      KDIAG (IBASE) = IDOF
      DO I=1,IDOF
        IBASE = IBASE + 1
        KDIAG(IBASE) = G(I)
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE GET_PCG_INFO (KDIAG,G,IDOF,IBASE) 
C
C     This gets IDOF and the G vector from KDIAG 
C
      INTEGER KDIAG (*), G(*)

      IBASE = IBASE + 1
      IDOF = KDIAG (IBASE)
      DO I=1,IDOF
        IBASE = IBASE + 1
        G(I) = KDIAG(IBASE) 
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE SOLVE_CG ( KV,B,X,N, KDIAG, NEL,ITERS_CG)
C
C     This solves the equations by using the Conjugate Gradient Method
C                      A x = B
C
C     .. The hard part is freedom steering .. each KM should really
C     hold with it its' own G vector so that we don't have to keep 
C     recreating it !  so we need to pass NUMS and NF to build G
C
      REAL        KV (*)        !- all the matrices
     +            ,B (*)        !- the input forces
     +            ,X (*)        !- the output disps
      INTEGER KDIAG  (*)        !- store of all the G vectors

C------------------------ Workspace arrays -----------------------------
      PARAMETER (MAX_N= 1000)   !-- for the workspace arrays
      REAL XNEW (MAX_N)         !- 'new' solution
     +       ,P (MAX_N)         !- working vector
     +       ,R (MAX_N)         !- working vector
     +       ,U (MAX_N)         !- temp work space for (KMxP)

      LOGICAL CONVERGED

c------------------------- set up the method -------------------------
c.. can do from within by "IF (ITERS=1) THEN ..."   :-)
c..... do this outside !
c      CALL NULVEC (X,N)
c      X (1) = 1.

c-----------------------------------------------------------------------
c.... this bit could also be in its own subroutine :-)
c.... given KV and X .. produces P
c     CALL MVMULT (A,IN,X,N,N,TEMP)   !- get the first guess

      CALL NULVEC (P,N)
      IBASE  = 1
      IBASE2 = 1
      DO IEL=1,NEL
        IDOF = KDIAG(IBASE)
        CALL SUB_BIT (KV(IBASE2),IDOF, IDOF, KDIAG(IBASE+1),X,P,N )
        IBASE = IBASE  + IDOF + 1        !- KDIAG pointer
        IBASE2= IBASE2 + IDOF * IDOF     !- KV pointer
      ENDDO
c-----------------------------------------------------------------------

      CALL VECSUB (B,P,R,N)           !- subtract from RHS
      CALL VECCOP (R,P,N)             !- and have 2 copies ( P and R )

c------------------------ loop the iterations --------------------------
      TOL = 10.**(-4)              !- adjust as necessary
      MAX_ITERS = N                !- Theoretical maximum needed  :-)
      MAX_ITERS = MAX_ITERS * 2    !- try this ?

      DO ITERS = 1 ,MAX_ITERS      

c-----------------------------------------------------------------------
c     CALL MVMULT(A,IN,P,N,N,U)     !- the 'hard bit' 
      CALL NULVEC (U,N)
      IBASE  = 1
      IBASE2 = 1
      DO IEL=1,NEL
        IDOF = KDIAG(IBASE)
        CALL SUB_BIT (KV(IBASE2),IDOF, IDOF, KDIAG(IBASE+1),P,U, N )
        IBASE = IBASE  + IDOF + 1        !- KDIAG pointer
        IBASE2= IBASE2 + IDOF * IDOF     !- KV pointer
      ENDDO
c-----------------------------------------------------------------------

      CALL VDOTV (R,R,UP,N)         !  UP   = R.R    R= residual?
      CALL VDOTV (P,U,DOWN,N)       !  DOWN = P.U    P= ? , U= ?
      ALPHA = UP/DOWN

      DO I=1,N
        XNEW(I) = X(I) + ALPHA * P(I)        !  X = X + à P
      ENDDO
      DO I=1,N
        R(I) = R(I) - ALPHA * U(I)           !  R = R - à U
      ENDDO

      CALL VDOTV  (R,R,UP2,N)            !  UP2 = R.R    (the new R)
      BETA = UP2/UP
      DO I=1,N
        P(I) = R(I) + BETA * P(I)
      ENDDO

c... now check for convergence beteen X and XNEW (can do 'earlier'?)
      CALL CHECON3 (XNEW ,X ,N,BIG,BIGGEST)
      CONVERGED = (BIGGEST.LT.TOL)         !-- convergence
c      WRITE(*,'(I3,A,F12.6,A,F12.6)') 
c     +       ITERS,': BIG=',BIG,' relative change =',BIGGEST

      IF (CONVERGED) GOTO 1             !- exit the loop

      CALL VECCOP (XNEW,X,N)       !- ready for the next pass
c-----------------------------------------------------------------------
      ENDDO
      CALL MYERROR (2,'Conjugant Gradient Method failed to converge')

C--------------------------- Ok converged ------------------------------
    1 CONTINUE
      WRITE(*,*) 'Iterations to convergence =', ITERS
      CALL VECCOP (XNEW,X,N)     !-- copy back into X (cf.'LOADS')
      END

C-----------------------------------------------------------------------
      SUBROUTINE SUB_BIT (KM,IKM,IDOF,G, LHS,RHS,N)
C
C     This multiplies a matrix by a vector (taking freedoms into account)
C     ie. given KM and G into multiples the 'bits' of LOADS by KM 
C     so just loop the elements, form KM and G and call this routine :-)
C
C     I could have KM being only the upper triangle to save storage ??
C
      REAL KM(IKM,IKM), LHS(*),RHS(*)
      INTEGER G(IDOF), G1, G2

      DO I=1,IDOF
        T = 0.
        G1 = G(I)
        IF (G1.NE.0) THEN
          DO J=1,IDOF
            G2 = G(J)
            IF (G2.NE.0) THEN
            T = T + KM(I,J) * LHS(G2)
            ENDIF
          ENDDO
          RHS(G1) = RHS(G1) + T           !- store the result
        ENDIF
      ENDDO
      END

C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C------------------- SPARIN matrix solving routines --------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      SUBROUTINE FORM_KDIAG (KDIAG,G,IDOF)
C
C     This finds the maximum bandwidth for each freedom (=FKDIAG)
C
      INTEGER KDIAG(*),G(*)
      DO 1 I=1,IDOF
        IWP1=1
        IF (G(I).EQ.0) GOTO 1    !- skip 'fixity'
        DO 2 J=1,IDOF
        IF (G(J).EQ.0) GOTO 2    !- skip 'fixity'
          IM=G(I)-G(J)+1           !- dummy
          IF(IM.GT.IWP1)IWP1=IM     !- ie. 'max' value
    2   CONTINUE
        K=G(I)                   !- 'dummy'
        IF (IWP1.GT.KDIAG(K)) KDIAG(K)=IWP1   !- ie 'max' value
    1 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE RESEQ_KDIAG (KDIAG,N,IKV)
C
C     This simply resequences the KDIAG pointers and checks the values
C
      INTEGER KDIAG(N)     !--- actualy may be longer than N

      IR = 0     !--------------- now reseqence the pointers ------------
      IBW_MIN  = 99999
      IBW_MAX  = 0
      DO I=1,N
        IW = KDIAG(I)
        IBW_MIN = MIN (IBW_MIN,IW)
        IBW_MAX = MAX (IBW_MAX,IW)
        IR = IR + IW
        KDIAG(I) = IR
      ENDDO
      IBW_MEAN = IR / N
      
      WRITE(*,'(3(A,I5))') 'MIN  Bandwidth=',IBW_MIN
     +                   ,' MAX  Bandwidth=',IBW_MAX 
     +                   ,' MEAN Bandwidth=',IBW_MEAN
      WRITE(*,'(A,I8,A,I8,A,F4.1,A)')
     +      'Size of KV =',IR,'/',IKV,' (',100.*IR/REAL(IKV) ,'%)'
      IF (IBW_MIN.eq.0) CALL MYERROR 
     + (1,'A zero bandwidth was found : ie. freedoms with no elements')
      IF (IR.GT.IKV) CALL MYERROR 
     +  (1,' The stiffnes matrix is too small')

      END
C-----------------------------------------------------------------------
      SUBROUTINE FORM_SPARCE_KV (KM,IKM,IDOF,G,KDIAG, KV)
C
c     This assembles the element stiffness matrix into the global matrix
c     stored as a vector accounting for a variable bandwidth (=FSPARV)
C
      REAL KV(*), KM(IKM,*)
      INTEGER KDIAG(*), G(*), G1, G2

      DO I=1,IDOF
        G1 = G(I)
        IF (G1.NE.0) THEN
          DO J=1,IDOF
            G2 = G(J)
            IF (G2.NE.0) THEN
              IW = G1-G2
              IF (IW.GE.0) THEN
                IVAL = KDIAG(G1)-IW
                KV(IVAL) = KV(IVAL)+KM(I,J)
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE SPARIN (A,N,KDIAG)
C
C     This performs Choleski reduction of the
C     variable-bandwidth stiffness matrix stored as a vector
C
      REAL A(*)
      INTEGER KDIAG(*)
      A(1) = SQRT(A(1))

      DO I=2,N                           !- loop the equations
        KI = KDIAG(I)-I
        L  = KDIAG(I-1)-KI+1
        DO 2 J=L,I
        X = A(KI+J)
        KJ = KDIAG(J)-J
        IF (J.EQ.1) GOTO 2               !- skip
        LBAR = MAX(L,KDIAG(J-1)-KJ+1)
C       LBAR = MAX0(L,LBAR)
        IF (LBAR.EQ.J) GOTO 2            !- skip
c       M=J-1
        DO K=LBAR,J-1
          X = X - A(KI+K)*A(KJ+K)
        ENDDO
    2   A(KI+J) = X/A(KJ+J)
        IF (X.LE.0) CALL MYERROR(2,'Degenerate Equations in SPARIN')
        A(KI+I) = SQRT(X)   
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE SPABAC (A,B,N,KDIAG)
C
C     This performs the Choleski back-substitution
c     on the variable bandwidth stiffness matrix
C
      REAL A(*),B(*)          !- B is the RHS and LHS !
      INTEGER KDIAG(*)
      B(1) = B(1)/A(1)        !- the first is easy :-)

      DO 1 I=2,N
        KI = KDIAG(I)-I         !- this row
        L  = KDIAG(I-1)-KI+1    !- previous row
        X = B(I)                !- get the 'force'  (factor ?)
        IF (L.EQ.I) GOTO 1      !- skip if zero-length (irrelevant!)
c       M=I-1                  
        DO J=L,(I-1)              
          X = X - A(KI+J)*B(J)  !- accumulate factor
        ENDDO
    1 B(I) = X/A(KI+I)          !- store as a 'disp' (well almost)

      DO 3 IT=2,N               !-  IT is NOT used !!
        I = N+2-IT              !- 'freedom to end at' = N-->2
        KI = KDIAG(I)-I         !- 'start' of this row ?
        X  = B(I)/A(KI+I)       !- extract, factor and hold_in_X
        B(I) = X
        L = KDIAG(I-1)-KI+1     !- start of the previous row ?
        IF (L.EQ.I) GOTO 3      !- skip if zero-length (irrelevant!)
c       M=I-1
        DO K=L,(I-1)            !- loop this row
          B(K) = B(K) - X*A(KI+K)     !- note X is B(I) saved
        ENDDO
    3 CONTINUE
      B(1)=B(1)/A(1)            !- and the first :-) (can do earlier?)
      END

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      > > > > > > >  Toolbox routines : JAC,BTDB   < < < < < < < 
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
      SUBROUTINE INVERT_JAC (JAC,IJAC,DET,NODOF)
C
C     This inverts a Jacobian (small square) matrix in 1D,2D,3D
C     this calls TWOBY2,TREEX3 as needed, putting JAC1 back into JAC
C
      REAL JAC(IJAC,*), JAC1(4,4)
      IF (NODOF.EQ.1) THEN
        DET = JAC(1,1)
        JAC(1,1) = 1./DET
      ELSEIF (NODOF.EQ.2) THEN
        CALL TWOBY2 (JAC,IJAC,JAC1,3,DET)
        CALL MATCOP (JAC1,3,JAC,IJAC,NODOF,NODOF)
      ELSEIF (NODOF.EQ.3) THEN
        CALL TREEX3 (JAC,IJAC,JAC1,3,DET)
        CALL MATCOP (JAC1,3,JAC,IJAC,NODOF,NODOF)
      ELSE
        CALL MATINV(JAC,IJAC,NODOF)
        PRINT*,'** WARNING, DET not calculated (INVERT)' 
      ENDIF
        IF(DET.LT.0.)PRINT*,'** WARNING negative Jacobian (INVERT)' 
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
      SUBROUTINE TREEX3 (JAC,IJAC,JAC1,IJAC1,DET)
C
C     This forms the inverse of a 3 by 3 matrix
C
      REAL JAC(IJAC,*),JAC1(IJAC1,*)
      DET=JAC(1,1)*(JAC(2,2)*JAC(3,3)-JAC(3,2)*JAC(2,3))
      DET=DET-JAC(1,2)*(JAC(2,1)*JAC(3,3)-JAC(3,1)*JAC(2,3))
      DET=DET+JAC(1,3)*(JAC(2,1)*JAC(3,2)-JAC(3,1)*JAC(2,2))
      JAC1(1,1)=JAC(2,2)*JAC(3,3)-JAC(3,2)*JAC(2,3)
      JAC1(2,1)=-JAC(2,1)*JAC(3,3)+JAC(3,1)*JAC(2,3)
      JAC1(3,1)=JAC(2,1)*JAC(3,2)-JAC(3,1)*JAC(2,2)
      JAC1(1,2)=-JAC(1,2)*JAC(3,3)+JAC(3,2)*JAC(1,3)
      JAC1(2,2)=JAC(1,1)*JAC(3,3)-JAC(3,1)*JAC(1,3)
      JAC1(3,2)=-JAC(1,1)*JAC(3,2)+JAC(3,1)*JAC(1,2)
      JAC1(1,3)=JAC(1,2)*JAC(2,3)-JAC(2,2)*JAC(1,3)
      JAC1(2,3)=-JAC(1,1)*JAC(2,3)+JAC(2,1)*JAC(1,3)
      JAC1(3,3)=JAC(1,1)*JAC(2,2)-JAC(2,1)*JAC(1,2)
      DO 1 K=1,3
      DO 1 L=1,3
      JAC1(K,L)=JAC1(K,L)/DET
    1 CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE MATINV (A,IA,N)
C
C     This forms the inverse of a matrix by Gauss-Jordan transformation
C
      REAL A(IA,*)
      DO 1 K=1,N
      CON=A(K,K)
      A(K,K)=1.
      DO 2 J=1,N
    2 A(K,J)=A(K,J)/CON
      DO 1 I=1,N
      IF(I.EQ.K)GOTO 1
      CON=A(I,K)
      A(I,K)=0.
      DO 3 J=1,N
    3 A(I,J)=A(I,J)-A(K,J)*CON
    1 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE FMBTDB (BEE,IBEE,DEE,IDEE,KM,IKM,FACT,IH,IDOF)
C
C     ** this is a mega-simplification of the long BT * DEE * BEE
C        and MSMULT, MATADD proceedures. Simply pass it the BEE and
C        DEE matrices and this willl do the rest.
C     It has condensed looping, it avoids unnecessary work when
C     zero terms are found in DEE ( hence not much slow-down when doing
C     SRI). It also exploits the symmetry of KM thus halving the effort.
C
C                                         Dan KIdger    May 1991
C
C     21-12-92 modified to use its own workspace
C
      REAL  BEE(IBEE,*), DEE(IDEE,*), KM(IKM,*), WORK(60)

      DO K=1,IH

C-------- first form one column of BT * DEE in the workspace -----------
        CALL NULVEC (WORK,IDOF)
        DO I=1,IH
          F = DEE(I,K) * FACT
          IF (ABS(F).GT. 1.E-12) THEN     !- skip blanks for speed
            DO J=1,IDOF
              WORK(J) =  WORK(J) + F * BEE(I,J)
            ENDDO
          ENDIF
        ENDDO

C--------  now product this column with BEE and sum into KM ------------
        DO I=1,IDOF
          DO J=I,IDOF
            KM(I,J) = KM(I,J) + WORK(I) * BEE(K,J)
            KM(J,I) = KM(I,J)                          !-- symmetry
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C           > > > > > >   Plasticity  Subroutines < < < < < <
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

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
      SUBROUTINE VMPL (E,V,STRESS,PL)
C
C      This forms the plastic matrix for a Von-Mises material in 2D
C
      REAL STRESS(*),TERM(4),PL(4,4)
      SX = STRESS(1)
      SY = STRESS(2)
      TXY= STRESS(3)
      SZ = STRESS(4)
      DSBAR = SQRT((SX-SY)**2+(SY-SZ)**2+(SZ-SX)**2+6.*TXY**2)/SQRT(2.)
      EE = 1.5*E/((1.+V)*DSBAR*DSBAR)
      TERM(1) = (2.*SX-SY-SZ)/3.
      TERM(2) = (2.*SY-SZ-SX)/3.
      TERM(3) = TXY
      TERM(4) = (2.*SZ-SX-SY)/3.
      DO I = 1,4
        DO J = I,4
          PL(I,J) = TERM(I)*TERM(J)*EE
          PL(J,I) = PL(I,J)   !-- symmetry
        ENDDO
      ENDDO
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE MOCOPL (PHI,PSI,E,V,STRESS,PL,IPL)
C
C     This forms the plastic stress/strain matrix        -- 2D
C     for a Mohr-Coulomb material  (phi,psi in degrees)
C
C     What on earth is this !!   Is it covered by my MC-flow ?
C     .. or is it the 'simple' matrix for initial stress ?
C
C     hmm 3D is fairly easy as 'uncoupled' :
C     just abstract more of the stresses and loop more of ROW and COL
C
      REAL STRESS(4), ROW(6), COL(6), PL(IPL,IPL)

      SX = STRESS(1)
      SY = STRESS(2)
      TXY = STRESS(3)
      SZ = STRESS(4)

      PI = 4.*ATAN(1.)
      PHIR = PHI*PI/180.
      PSIR = PSI*PI/180.
      SNPH = SIN(PHIR)
      SNPS = SIN(PSIR)
      SQ3 = SQRT(3.)
      CC = 1.-2.*V
      DX = (2.*SX-SY-SZ)/3.
      DY = (2.*SY-SZ-SX)/3.
      DZ = (2.*SZ-SX-SY)/3.
      D2 = SQRT(-DX*DY-DY*DZ-DZ*DX+TXY*TXY)
      D3 = DX*DY*DZ-DZ*TXY*TXY
      TH = MIN (MAX (-1.,-3.*SQ3*D3/(2.*D2**3)),1. )

C      IF (TH.GT. 1.) TH = 1.
C      IF (TH.LT.-1.) TH = -1.
C      TH = ASIN(TH)/3. 

      SNTH = SIN(ASIN(TH)/3.)
      IF (ABS(SNTH).GT..49) THEN     !-- close to the 'Kink'
        SIG = -1.
        IF (SNTH.LT.0.) SIG = 1.
        RPH = SNPH*(1.+V)/3.
        RPS = SNPS*(1.+V)/3.
        CPS = .25*SQ3/D2*(1.+SIG*SNPS/3.)
        CPH = .25*SQ3/D2*(1.+SIG*SNPH/3.)
        COL(1) = RPH+CPH*((1.-V)*DX+V*(DY+DZ))
        COL(2) = RPH+CPH*((1.-V)*DY+V*(DZ+DX))
        COL(3) = CPH*CC*TXY
        COL(4) = RPH+CPH*((1.-V)*DZ+V*(DX+DY))
        ROW(1) = RPS+CPS*((1.-V)*DX+V*(DY+DZ))
        ROW(2) = RPS+CPS*((1.-V)*DY+V*(DZ+DX))
        ROW(3) = CPS*CC*TXY
        ROW(4) = RPS+CPS*((1.-V)*DZ+V*(DX+DY))
        EE = E/((1.+V)*CC*(RPH*SNPS+2.*CPH*CPS*D2*D2*CC))
      ELSE                       !------ 'normal' state of affairs
        ALP = ATAN(ABS((SX-SY)/(2.*TXY)))
        CA = COS(ALP)
        SA = SIN(ALP)
        DD = CC*SA
        S1 = 1.
        S2 = 1.
        IF ((SX-SY).LT.0.) S1 = -1.
        IF ( TXY.   LT.0.) S2 = -1.
        COL(1) = SNPH+S1*DD
        COL(2) = SNPH-S1*DD
        COL(3) = S2*CC*CA
        COL(4) = 2.*V*SNPH
        ROW(1) = SNPS+S1*DD
        ROW(2) = SNPS-S1*DD
        ROW(3) = S2*CC*CA
        ROW(4) = 2.*V*SNPS
        EE = E/(2.*(1.+V)*CC*(SNPH*SNPS+CC))
      END IF
      DO I = 1,4
        DO J = 1,4
          PL(I,J) = EE*ROW(I)*COL(J)
        ENDDO
      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE GET_FLOW (DSBAR,THETA,PSI,E,V,PHI,F,STRESS,NODOF,
     +  FLOW,IFLOW,IOP_VM)
C
C    This calculates the 'flow' matrix for a viscoplastic material
C     if IOP_VM = 1  : for Mohr-Coulomb
C     if IOP_VM = 2  : for Von Mises' model
C
      PARAMETER (IM = 6)    !- max # of stress terms

      REAL M1(IM,IM), M2(IM,IM), M3(IM,IM)  !-  componenent matrices
     +     ,FLOW (IFLOW,IFLOW)              !-  --> into flow matrix
     +    ,STRESS (*)                        !-   the given stress

      DTR = 3.14159265 / 180.   !-- factor for 'degress to radians'

      IF (NODOF.EQ.2) THEN
        IH = 4                               !- be careful !
        CALL FORMM  (STRESS,M1,M2,M3,IM)
      ELSEIF (NODOF.EQ.3) THEN
        IH = 6
        CALL FORMM3 (STRESS,M1,M2,M3,IM)
      ELSE
        STOP 'MC_FLOW - NOT 2D OR 3D'
      ENDIF

      IF (IOP_VM.EQ.1) THEN                         !- Mohr-Coulomb
        CALL MOCOUQ (PSI,DSBAR,THETA,DQ1,DQ2,DQ3)  
        DT  =  4.*(1.+V)*(1.-2.*V)  / (E*(1.-2.*V+SIN(PHI*DTR)**2 ))
      ELSEIF (IOP_VM.EQ.1) THEN
        DQ1 = 0.
        DQ2 = 1.5/DSBAR      !-- for Von Mise  only
        DQ3 = 0.
        DT  =  4.*(1.+V) / (3.*E)  !- Von Mise
      ELSE
        CALL MYERROR (2,'Unknown Von Mise/Mohr-Coulomb rule')
      ENDIF
      DO J=1,IH
        DO I=1,IH   
          FLOW(I,J) = F*DT*  (M1(I,J)*DQ1 +M2(I,J)*DQ2 +M3(I,J)*DQ3)
        ENDDO
      ENDDO

      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE MOCOUF (PHI,C,SIGM,DSBAR,THETA,F)
C
C     This calculates the value of the yield function
C     for a mohr-coulomb material (phi in degrees)
C       --- ie. < 0. is not failed; > 0. is ??
C
      PHIR = PHI*4.*ATAN(1.)/180.
      SNPH = SIN (PHIR)
      CSPH = COS (PHIR)
      SNTH = SIN (THETA)
      CSTH = COS (THETA)
      F = SNPH*SIGM + DSBAR*(CSTH/SQRT(3.)-SNTH*SNPH/3.) -C*CSPH
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE MOCOUQ (PSI,DSBAR,THETA,DQ1,DQ2,DQ3)
C
c     This forms the derivatives of a Mohr-Coulomb potential function 
c     with respect to the three invariants     (psi in degrees) 
C
      PSIR = PSI*4.*ATAN(1.)/180.
      SNTH = SIN(THETA)
      SNPS = SIN(PSIR)
      SQ3 = SQRT(3.)
      DQ1 = SNPS
      IF (ABS(SNTH).GT..49) THEN       !- close to the 'cusp' ?
        C1 = 1.
        IF(SNTH.LT.0.)C1 = -1.
        DQ2 = (SQ3*.5-C1*SNPS*.5/SQ3)*SQ3*.5/DSBAR
        DQ3 = 0.
      ELSE
        CSTH = COS(THETA)
        CS3TH = COS(3.*THETA)
        TN3TH = TAN(3.*THETA)
        TNTH = SNTH/CSTH
        DQ2 = SQ3*CSTH/DSBAR*((1.+TNTH*TN3TH)+SNPS*(TN3TH-TNTH)/SQ3)*.5
        DQ3 = 1.5*(SQ3*SNTH+SNPS*CSTH)/(CS3TH*DSBAR*DSBAR)
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FORMM (STRESS,M1,M2,M3,IM)
C
C     This forms the derivatives of the invariants  
C     with respect to the stresses  - 2D
C
      REAL STRESS(*),M1(IM,IM), M2(IM,IM), M3(IM,IM)

      CALL NULL (M1,IM,IM,IM)
      CALL NULL (M2,IM,IM,IM)
      CALL NULL (M3,IM,IM,IM)

      SX = STRESS(1)
      SY = STRESS(2)
      TXY= STRESS(3)
      SZ = STRESS(4)
      DX = (2.*SX-SY-SZ)/3.
      DY = (2.*SY-SZ-SX)/3.
      DZ = (2.*SZ-SX-SY)/3.
      SIGM = ( SX+SY+SZ)/3.

      M1(1,1) = 1.       !-- hmm could really set all to 1.
      M1(1,2) = 1.       !-- then null the 3rd row & col
      M1(2,1) = 1.
      M1(1,4) = 1.
      M1(4,1) = 1.
      M1(2,2) = 1.
      M1(2,4) = 1.
      M1(4,2) = 1.       !- note: only 7/16 zeros 
      M1(4,4) = 1.
      DO I = 1,4
        DO J = 1,4
          M1(I,J) = M1(I,J)/(9.*SIGM)
        ENDDO
      ENDDO
      M2(1,1) = .6666666666666666     !- diagonal    (FORMM3 does /3.
      M2(2,2) = .6666666666666666                 !   at the end :-)
      M2(4,4) = .6666666666666666
      M2(3,3) = 2.                    !- shear term
      M2(2,4) = -.3333333333333333
      M2(4,2) = -.3333333333333333    !- off-diagonals
      M2(1,2) = -.3333333333333333
      M2(2,1) = -.3333333333333333
      M2(1,4) = -.3333333333333333    !- note: only 6/16 zeros !
      M2(4,1) = -.3333333333333333

      M3(1,1) = DX/3.
      M3(2,4) = DX/3.
      M3(4,2) = DX/3.
      M3(2,2) = DY/3.
      M3(1,4) = DY/3.
      M3(4,1) = DY/3.
      M3(4,4) = DZ/3.
      M3(1,2) = DZ/3.
      M3(2,1) = DZ/3.
      M3(3,3) =-DZ
      M3(3,4) = TXY*(-2./3.)            !- bracketsd added 18-3-93
      M3(4,3) = TXY*(-2./3.)
      M3(1,3) = TXY/3.
      M3(3,1) = TXY/3.
      M3(2,3) = TXY/3.                 !- note: only 6/16 zeros
      M3(3,2) = TXY/3.
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE FORMM3 (STRESS,M1,M2,M3,IM)
C
C     This forms the derivatives of the invariants
C     with respect to the stresses  - 3D
C
C  ---> OK for 2d/3d... replace '3' with NODOF (careful with axisym)
C
      REAL STRESS(*), M1(IM,IM), M2(IM,IM), M3(IM,IM)

      CALL NULL (M1,IM,IM,IM)
      CALL NULL (M2,IM,IM,IM)
      CALL NULL (M3,IM,IM,IM)

      SX  = STRESS(1)
      SY  = STRESS(2)
      SZ  = STRESS(3)
      TXY = STRESS(4)
      TYZ = STRESS(5)
      TZX = STRESS(6)
      SIGM= (SX+SY+SZ)/3.
      DX = SX-SIGM         !--- hmm only the deviatorics are relevant
      DY = SY-SIGM
      DZ = SZ-SIGM

      DO I = 1,3
        DO J = 1,3         !-- only J=I:3 is needed (symmetry)
          M1(I,J) = 1./(3.*SIGM)   !- direct-stres terms
        ENDDO
      ENDDO
      DO I = 1,3
        M2(I  ,I  ) = 2.           !- direct-stress terms
        M2(I+3,I+3) = 6.           !- shear -stress terms
      ENDDO
      M2(1,2) = -1.                !- off-diagonal terms
      M2(1,3) = -1.
      M2(2,3) = -1.

      M3(1,1) = DX                 !- loop ?
      M3(1,2) = DZ
      M3(1,3) = DY
      M3(1,4) = TXY
      M3(1,5) = -2.*TYZ
      M3(1,6) = TZX                !- or use IF round the 3d bits
      M3(2,2) = DY
      M3(2,3) = DX
      M3(2,4) = TXY
      M3(2,5) = TYZ
      M3(2,6) = -2.*TZX
      M3(3,3) = DZ
      M3(3,4) = -2.*TXY
      M3(3,5) = TYZ
      M3(3,6) = TZX
      M3(4,4) = -3.*DZ
      M3(4,5) =  3.*TZX
      M3(4,6) =  3.*TYZ
      M3(5,5) = -3.*DX
      M3(5,6) =  3.*TXY
      M3(6,6) = -3.*DY
      DO I = 1,6
        DO J = I,6
          M1(I,J) = M1(I,J)/3.
          M1(J,I) = M1(I,J)       !-- symmetry
          M2(I,J) = M2(I,J)/3.    !-- and scaling 
          M2(J,I) = M2(I,J)       
          M3(I,J) = M3(I,J)/3.
          M3(J,I) = M3(I,J)       
        ENDDO
      ENDDO
      RETURN
      END

C-------------------------------------------------------------------
C-------------------------------------------------------------------
C    > > > > >     Some  mesh and freedom utilities     < < < < < 
C-------------------------------------------------------------------
C-------------------------------------------------------------------
      SUBROUTINE SORTNF (NF,INF,NN,NODOF,N)
C
C     this resequences all the non-zero terms in NF to an ascending 
C     sequence and returns the total number of freedoms in N
C
      INTEGER NF(INF,*)
      N = 0
      DO I=1,NN
        DO J=1,NODOF
          IF (NF(J,I).NE.0) THEN
            N = N + 1
            NF(J,I) = N
          ENDIF
        ENDDO
      ENDDO
      END
C-----------------------------------------------------------------------
      SUBROUTINE NUM2G (NF,INF,NOD,NUM,NODOF,G)
C
C     This copies the NF data into the steering G (via NUM)
C
      INTEGER NF(INF,*),NUM(*),G(*)
      IC=0
      DO I=1,NOD
        DO J=1,NODOF
          IC=IC+1
          G(IC) = NF(J,NUM(I))
        ENDDO
      ENDDO
      END

C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C   > > > > > > >  Basic File Handling routines < < < < < < < <  
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      SUBROUTINE GET_KEYWORD (IO,IO_BASE,KEYWORD)
C
C     This gets the next keyword from the data file
C     '*EOF' is returned when all data is exhausted
C
      CHARACTER KEYWORD*(*), LINE*80

    1 READ (IO,'(A)',IOSTAT=IOS) LINE
      CALL IN_TEST (IO,IOS,*1,*999)  
      IF (IO.LT.IO_BASE)    GOTO 999
      IF (LINE(1:1).NE.'*') GOTO 1
      KEYWORD  = LINE(1:INDEX(LINE,' '))    !-- strip off any comments
      CALL UPCASE (KEYWORD)                 !-- Uppercase is easier to handle
      RETURN
 999  KEYWORD ='*EOF'         !- end of the all the data
      END
c-----------------------------------------------------------------------
      SUBROUTINE IN_TEST (IO,IOSTAT,*,*)
C
C     this tests the nature of a data read error
C     return to *1 if 'soft', and *2 if 'fatal'
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
c        WRITE (*,'(A)') LINE
c        RETURN 1
      ELSEIF (LINE(1:1).EQ.'%') THEN !----An 'echoed' comment line -----
        WRITE (90,'(A)') LINE        !       to the '.PLT' file
        RETURN 1
      ELSEIF (LINE(1:1).eq.'*') THEN !------- A new Keyword -----------
        BACKSPACE (IO)  ! point at the '*'
        RETURN 2
      ELSEIF (LINE(1:1).eq.'<') THEN !------- read_from_another_file ---
        IO = IO + 1
        OPEN (IO,FILE = LINE(2:INDEX(LINE,'  ')),IOSTAT=IOSTAT)
        IF (IOSTAT.NE.0) THEN        !- error= 'new data file not found'
          CALL MYERROR (2,'missing data file-'//LINE)
          RETURN 2    !-- crash out 
        ENDIF
        RETURN 1      !-- this line was missing !!!  14-9-92
c- - - - - - - - - - - - - - - - - - - 
      ELSE                    !--------- 'rubbish' found in file -------
        CALL MYERROR(2,'error in data file format ')
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
C-----------------------------------------------------------------------
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
      SUBROUTINE R_ADD_REM_ELEMS (IO,NUMS,INUMS,NEL)
C
C     This simply switches on or off the given elements
C
      INTEGER NUMS (INUMS,NEL)    !- ie. >=NEL :-)

      DO ILOOP=1,9999
    1   READ (IO,IOSTAT=IOS) JMAT
        CALL IN_TEST (IO,IOS,*1,*999)

        IC = 0
        DO IEL=1,NEL
          CALL GET_ELEMENT 
     +    (NUMS,INUMS,NEL, NUM, NEN,NDIME,ITYPE, IMAT,IU1,IU2,IU3)

          IF (IMAT.EQ.ABS(JMAT)) THEN
            IC = IC + 1
            CALL ADD_ELEMENT 
     +      (NUMS,INUMS,NEL, NUM, NEN,NDIME,ITYPE, JMAT,IU1,IU2,IU3)
          ENDIF
        ENDDO
        PRINT*, ILOOP,' :',JMAT,' : elements found =',IC
      ENDDO

 999  RETURN
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
C
      REAL GC(IGC,*),COORD(ICOORD,*)
      INTEGER NUM(*)
      DO I=1,NOD
        DO J=1,NDIM
          COORD(I,J) = GC(J,NUM(I))
        ENDDO
      ENDDO
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
      END

C---------------- end of Subroutine 'families' -------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------


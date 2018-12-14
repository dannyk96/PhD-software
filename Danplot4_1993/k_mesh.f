C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C    >>>> 'Self-Contained' routines for reading various mesh file types
C          .. Dan Kidger  Sept'93
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE KEY_MESH_READ
     +          (FOUND,KEYWORD,IO,GC,IGC,NDIM,NN,NUMS,INUMS,NEL)
C
C     A Keyword handler for reading/writing the basic mesh
C     DJK 24-5-93
C
      REAL GC (IGC,*)
      INTEGER NUMS(INUMS,*), IO
      LOGICAL FOUND
      CHARACTER KEYWORD*(*), FILE*70
      DATA IO2/80/

      FOUND = .TRUE.

c------------------------- Mesh Reading -------------------------------
C...... 'Old' style mesh data ....
      IF (KEYWORD.EQ.'*BASIC_SLOPE_MESH') THEN    !- 'old' style
        CALL R_READGM (IO,GC,IGC,NN,NDIM,NUMS,INUMS,NEL)

      ELSEIF (KEYWORD.EQ.'*ELEMENT_TYPES') THEN   !- 'old' style
        CALL R_GETELT (IO,NUMS,INUMS,NEL,NDIM)

c....... and complete mesh descriptions ...
      ELSEIF (KEYWORD.EQ.'*IMPORT_OFF') THEN
        CALL R_OFF (IO,GC,IGC,NN,NDIM,NUMS,INUMS,NEL)

      ELSEIF (KEYWORD.EQ.'*IMPORT_PL') THEN               !- synonym
        CALL R_DANPLOT (IO,GC,IGC,NDIM,NN,NUMS,INUMS,NEL)
      ELSEIF (KEYWORD.EQ.'*READ_DANPLOT') THEN            !- synonym
        CALL R_DANPLOT (IO,GC,IGC,NDIM,NN,NUMS,INUMS,NEL)

      ELSEIF (KEYWORD.EQ.'*IMPORT_SRF') THEN      !- 'SRF' format
        CALL R_SRF (IO,GC,IGC,NN,NDIM,NUMS,INUMS,NEL)

      ELSEIF (KEYWORD.EQ.'*IMPORT_3DEDIT') THEN   !- '3D EDITF' format
        CALL R_3DEDIT (IO,GC,IGC,NN,NDIM,NUMS,INUMS,NEL)


C-------------------------- Mesh Writing -------------------------------
C... output file is automaticaly generated from the current input file
C... by changing the file extension 
C-- What if the file has already been opened by the main program ?
C-- so test for the current file being open ? and only open it
C-- if it has NOT already been opened ?

      ELSEIF (KEYWORD.EQ.'*WRITE_DANPLOT') THEN
        INQUIRE (UNIT=IO,NAME=FILE)
        OPEN(IO2,FILE=FILE(1:INDEX(FILE,'.'))//'PL')
        CALL WR_DANPLOT (IO2,GC,IGC,NDIM,NN,NUMS,INUMS,NEL)
        CLOSE(IO2)

      ELSEIF (KEYWORD.EQ.'*WRITE_OFF') THEN
        INQUIRE (UNIT=IO,NAME=FILE)
        OPEN (IO2,FILE=FILE(1:INDEX(FILE,'.'))//'GEO')
        CALL WR_OFF (IO2,GC,IGC,NUMS,INUMS,NDIM,NN,NEL)
        CLOSE(IO2)

      ELSEIF (KEYWORD.EQ.'*WRITE_RAYSHADE') THEN  
        INQUIRE (UNIT=IO,NAME=FILE)
        print*,'file=',file
        OPEN (IO2,FILE=FILE(1:INDEX(FILE,'.'))//'RAY')
        CALL WR_RAYSHADE (IO2,GC,IGC,NUMS,INUMS,NDIM,NN,NEL)
        CLOSE(IO2)

C------------------------- else unknown --------------------------------
      ELSE
        FOUND = .FALSE.
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
c.. ( R_DANPLOT/WR_DANPLOT were here 6 July 93)
c... now returned here 10-9-93 and 'other' file formats added too
C-----------------------------------------------------------------------
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
      SUBROUTINE R_OFF (IO,GC,IGC,NN,NDIM,NUMS,INUMS,NEL)
C
C     This imports a mesh in 'Object File Format'
C     ( These are 2d polygons (mostly 3 node triangles) in a 3d mesh
C     .. no need to test for comments as they never have any!
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*), NUM(30)

      NDIM = 3         !- force into 3D
      NN_OLD  = NN
      NEL_OLD = NEL
    1 READ (IO,*,IOSTAT=IOS)  NN_NEW,NEL_NEW,NEDGES
      CALL IN_TEST(IO,IOS,*1,*999)
      READ (IO,*)  ((GC(J,I+NN_OLD),J=1,ndim),I=1,NN)
      WRITE(*,'(A,I4,A)')  '>> ',NN_NEW,' nodes read'
      NN = NN_OLD + NN_NEW


      NDIME = 2   !-- always 2D facets
      ITYPE = 9   !-- use type 9 for OFF format (ie. do disps, and 'flat')

      DO IEL=1,NEL
        READ (IO,*) NOD,(NUM(J),J=1,NOD)

        DO J=1,NOD
          NUM(J) = NUM(J) + NN_OLD   !-- offsets
        ENDDO
        CALL PUT_ELEMENT 
     +  (NUMS,INUMS,NEL_OLD+IEL, NUM, NOD,NDIME,ITYPE, IMAT,IEL,1,1)

      ENDDO
      CLOSE (IO)
      IO = IO - 1     !-- exit from this file  :-)
      WRITE(*,'(A,I4,A)')  '>> ',NEL,' elements read'
  999 CONTINUE
      NEL= NEL_NEW
      END

C----------------------------------------------------------------------
      SUBROUTINE WR_OFF (IO,GC,IGC,NUMS,INUMS,NDIM,NN,NEL)
C
C     This writes out the mesh in 'Object_File_Format' format
C   .. need to add a zero 'z' component if only 2d
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*), IO

      WRITE(IO,'(I5)') NN,NEL,12345
      IF (NDIM.EQ.2) THEN
        DO I=1,NN
          WRITE(IO,'(3F12.5)') (GC(J,I),J=1,NDIM),0.
        ENDDO
      ELSEIF (NDIM.EQ.3) THEN
        DO I=1,NN
          WRITE(IO,'(3F12.5)') (GC(J,I),J=1,NDIM)
        ENDDO
      ENDIF
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

        CALL PUT_ELEMENT 
     +  (NUMS,INUMS,NEL_OLD+IEL, NUM, NEN,NDIME,ITYPE, IMAT,IEL,1,1)

      ENDDO
c     CLOSE (IO)
      WRITE(*,'(A,I4,A)')  '>> ',NEL_NEW,' elements read'
      NEL= NEL_OLD + NEL_NEW
  999 CONTINUE
      END

C-----------------------------------------------------------------------
      SUBROUTINE WR_DANPLOT (IO,GC,IGC,NDIM,NN,NUMS,INUMS,NEL)
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
      SUBROUTINE R_SRF (IO,GC,IGC,NN,NDIM,NUMS,INUMS,NEL)
C
C     This imports a mesh in 'SRF Format'
C     ( These are 2d polygons (mostly 3 node triangles) in a 3d mesh
C     commesnt lines begin with '#'
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*), NUM(30), IO

      NN_OLD  = NN
      NEL_OLD = NEL

        NDIM  = 3
        NDIME = 2   !-- always 2D facets
        ITYPE = 9   !-- use type 9 for OFF format (ie. do disps, and 'flat')

      READ (IO,*,iostat=ios) id               !- dummy for the title line 
      CALL IN_TEST(IO,IOS,*2,*999)
   2  continue ! -fudged continue !

      READ (IO,*)    !- dummy for the title line
      READ (IO,*) id !- dummy for the version number
      READ (IO,*)  nmats,NN_new,NEL_new, id,id

      do i=1,nmats
        READ (IO,*)    !- dummy for the material properties
      enddo 

    1   READ (IO,*)   ((GC(J,I+NN_OLD),J=1,3),I=1,NN_NEW)
        CALL IN_TEST(IO,IOS,*1,*999)
 
        WRITE(*,'(A,I4,A)')  '>> ',NN_NEW,' nodes read'

        NN = NN_OLD + NN_NEW

      DO IEL=1,NEL_NEW
        READ (IO,*) NOD,IMAT,(NUM(J),J=NOD,1,-1) !-reverse order

        DO J=1,NOD
          NUM(J) = NUM(J) + NN_OLD   !-- offsets
        ENDDO
        CALL PUT_ELEMENT 
     +  (NUMS,INUMS,NEL_OLD+IEL, NUM, NOD,NDIME,ITYPE, IMAT,IEL,1,1)

      ENDDO
      CLOSE (IO)
      IO = IO - 1     !-- exit from this file  :-)
      WRITE(*,'(A,I4,A)')  '>> ',NEL_NEW,' elements read'
  999 CONTINUE
      NEL= NEL_OLD + NEL_NEW
      END

C-----------------------------------------------------------------------
      SUBROUTINE R_3DEDIT (IO,GC,IGC,NN,NDIM,NUMS,INUMS,NEL)
C
C     This imports a mesh in '3D_EDIT Format'
C     ( These are 2d polygons (3nt or 4nq's) in a 3d domain
C     order is NN, N_TRI,N_QUAD, then nodes x,y,z
C     then quads then triangles NUM, (then color, and '2','0')
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*), NUM(30) , IO

      NN_OLD  = NN
      NEL_OLD = NEL

        NDIM  = 3
        NDIME = 2   !-- always 2D facets
        ITYPE = 9   !-- use type 9 for OFF format (ie. do disps, and 'flat')

    1 READ (IO,*,IOSTAT=IOS)  NN_NEW, NEL_TRI,NEL_QUAD
      CALL IN_TEST(IO,IOS,*1,*999)

c...................... nodal coordinates ...........................
      READ (IO,*)   ((GC(J,I+NN_OLD),J=1,3),I=1,NN_NEW)
      WRITE(*,'(A,I4,A)')  '>> ',NN_NEW,' nodes read'
      NN = NN_OLD + NN_NEW

c...................... element steering ...........................
      NEL_NEW = NEL_TRI+NEL_QUAD
      NOD = 4
      DO IEL=1,NEL_NEW
        IF (IEL.gt.NEL_QUAD) NOD = 3
        READ (IO,*) (NUM(J),J=NOD,1,-1) , IMAT, id,id       

        DO J=1,NOD
          NUM(J) = NUM(J) + NN_OLD   !-- offsets
        ENDDO
        CALL PUT_ELEMENT 
     +  (NUMS,INUMS,NEL_OLD+IEL, NUM, NOD,NDIME,ITYPE, IMAT,IEL,1,1)

      ENDDO
      CLOSE (IO)
      IO = IO - 1     !-- exit from this file  :-)
      WRITE(*,'(A,I4,A)')  '>> ',NEL_NEW,' elements read'
  999 CONTINUE
      NEL= NEL_OLD + NEL_NEW
      END

C-----------------------------------------------------------------------
      SUBROUTINE WR_RAYSHADE (IO,GC,IGC,NUMS,INUMS,NDIM,NN,NEL)
C  
C     This writes out the mesh in the RAYSHADE raytracing package format
C
C    note: this format is essiantialy just a set of polygon descriptions
C    ... only really usefuly for 'curved' objects.. eg OFF.
C    ... therefore need a converter to write out a DANPLOT 3D mesh
C    ... as the OFF file of its Facets
C
C      Dan Kidger   21 July 1993
C
      REAL GC(IGC,*)
      INTEGER NUMS(INUMS,*), IO,  NUM (27)

      DO IEL=1,NEL
        CALL GET_ELEMENT
     +  (NUMS,INUMS,iEL, NUM, NOD,NDIME,ITYPE, IMAT,I1,i2,i3)

        WRITE(IO,'(A, 99f11.3)')   'poly', 
     +  ((GC(jj,NUM(J)),jj=1,ndim),J=NOD,1,-1)  !- reverse order

      ENDDO

      WRITE(*,'(A,I4,A)')  '>> ',NEL,' elements written'
      RETURN
      END  
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE R_READGM(IO,GC,IGC,NN,NDIM,NUMS,INUMS,NEL)
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
            CALL SLOGE2(IP,IQ,NXS,NYS,NYE,TOP,BOT,DEPTH,COORD,20,NUM)
          ELSEIF (NDIM.EQ.3) THEN           !-- 8nq only
            CALL SLOGET(IP,IQ,IS,NXS,NYS,NYE,NZE,TOP,BOT,DEPTH,BRDTH,
     +      COORD,20,NUM,NOD)      !-- 8nb/14nb/20nb
          ELSE 
            CALL MYERROR(2,'NDIM not valid (not 2 or 3)')  
          ENDIF
        NDIME = NDIM  !------------- update the element table ----------
        IMAT  = 1
        ITYPE = 1
        CALL PUT_ELEMENT 
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


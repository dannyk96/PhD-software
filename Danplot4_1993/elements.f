c... this is the DANPLOT version .. 
c     with ADD_ELEMENT changed to PUT_ELEMENT
c     and GET/PUT _EL_IMAT appended
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C    > > > > > > > >     Element Database routines < < < < < < < < <
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

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
            CALL PUT_ELEMENT 
     +      (NUMS,INUMS,NEL, NUM, NEN,NDIME,ITYPE, JMAT,IU1,IU2,IU3)
          ENDIF
        ENDDO
        PRINT*, ILOOP,' :',JMAT,' : elements found =',IC
      ENDDO

 999  RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE PUT_ELEMENT (NUMS,INUMS,IEL, NUM, NOD,NDIME,ITYPE
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
      SUBROUTINE GET_EL_IMAT (NUMS,INUMS,IEL, IMAT)
C
C     simply returns the material number of a given element (cf GET_ELEMENT)
C
      INTEGER NUMS(INUMS,*)

      IMAT    = NUMS (INUMS-3,IEL) 
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE PUT_EL_IMAT (NUMS,INUMS,IEL, IMAT)
C
C     simply sets the material number of a given element (cf PUT_ELEMENT)
C
      INTEGER NUMS(INUMS,*)

      NUMS (INUMS-3,IEL) = IMAT
      RETURN
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


C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C   > > > > > > >  Basic File Handling routines < < < < < < < <  
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
c... also provide routines for getting the command-line arguements
c... and using SELECT_FILE() if desired

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
C     This tests the nature of a data read error
C     return to *1 if 'soft', and *2 if 'fatal'
C     21-9-93 fixed '< file', io status !
C
      CHARACTER LINE*80
      LOGICAL FILE_OPEN

      IF (IOSTAT.EQ.0) THEN
        RETURN                          !-------- No error! -----------
      ELSEIF (IOSTAT.LT.0) THEN
        IO = IO - 1           !------ try the previous unit --------
        INQUIRE (IO,OPENED=FILE_OPEN)       !*NOT* 'EXIST' (cos always true)
        IF (FILE_OPEN) THEN   !------------ 'end-of-file'----------------
          RETURN 1              !('soft' cos previous unit is open)
        ELSE                  !---------- 'file was not open'------------
          RETURN 2              !(hence all files must have been done)
        ENDIF
      ELSE     !- check other possible read-errors

        BACKSPACE (IO)           !(should never get an error here)
        READ (IO,'(A)') LINE     !(or here either)

        IF (LINE(1:1).EQ.'='      !--- allowable comment lines 
     +    .OR.  LINE(1:1).EQ.'#'
     +    .OR.  LINE(1:1).EQ.'c'
     +    .OR.  LINE(1:1).EQ.'C'
     +    .OR.  LINE(1:1).EQ.' '
     +    .OR.  LINE(1:1).EQ.'/') THEN
          RETURN 1
c        ELSEIF (LINE(1:1).EQ.'!') THEN !----An 'echoed' comment line -----
c          WRITE (*,'(A)') LINE
c          RETURN 1
        ELSEIF (LINE(1:1).EQ.'%') THEN !----An 'echoed' comment line -----
          WRITE (90,'(A)') LINE        !       to the '.PLT' file
          RETURN 1
        ELSEIF (LINE(1:1).eq.'*') THEN !------- A new Keyword -----------
          BACKSPACE (IO)               !- point back to the '*'
          RETURN 2
        ELSEIF (LINE(1:1).eq.'<') THEN !------- read_from_another_file ---
          IO = IO + 1
          OPEN (IO,FILE = LINE(3:INDEX(LINE,'   '))
     +     ,IOSTAT=IOSTAT,STATUS='old')
          IF (IOSTAT.NE.0) THEN        !--- 'new data file not found'
            CALL MYERROR (2,'missing data file-'//LINE)
            RETURN 2    !-- crash out 
          ENDIF
          RETURN 1      !-- this line was missing !!!  14-9-92
c-------------------- rubbish found in the file -----------------------
        ELSE    
          print*,'>>LINE was :',LINE(1:60)
          CALL MYERROR(2,'>>Error in data file format ')
          RETURN 2
        ENDIF
      ENDIF      !- test differnt values of IOSTAT
c  99 CONTINUE  
        STOP 'shouldn''t have got to this linbe in IN_TEST'
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
        STOP
      ELSEIF (ILEVEL.EQ.3) THEN
        WRITE (*,'(A,2X,A)') '**INTERNAL ERROR** :',STRING
        STOP
      ELSE
        WRITE (*,'(A,I2)') 'hmm.. unknown error level =',ILEVEL
      ENDIF
      END


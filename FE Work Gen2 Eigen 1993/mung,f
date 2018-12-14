*$noextensions
      PROGRAM MUNG
c
c     this will convert a Salford FTN77 program into
c     one more suited to for WATFRO to handle
c    ie. strip off '!' comments
c    and put '*$noextensions' at the top to disable non-ansi warnings
c
      CHARACTER LINE*80

      OPEN (10,FILE = 'IN.FOR',STATUS='OLD')
      OPEN (11,FILE = 'OUT.FOR')

c---- remove the FTN77 'options' line
      READ(10,*)

      WRITE(11,'(A)') '*$noextensions'
      DO I=1,99999
        READ(10,'(A)',END=999) LINE

c------------- remove the '@' symbol -----------------
        IPOS = INDEX(LINE,'@')       
        IF (IPOS.GE.1) LINE(IPOS:IPOS) = '_'

c------------- strip off trailing comments -----------------
        IPOS = INDEX (LINE,'!') -1   
        IF (IPOS.LE.0) IPOS = 80
        IF (IPOS.EQ.INDEX(LINE,'!''')) IPOS = 80

        WRITE (11,'(A)') LINE(1:IPOS)
      ENDDO
  999 CONTINUE
      END


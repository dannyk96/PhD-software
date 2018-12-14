
#        BOOK 6.1 (slope stability) test problem
#        -------------------------------
#
#        ... try MUNGing into 20nb's   !


*READ_DANPLOT
  2
 13
  1     0.  105.
  2    51.  105.
  3    51.    0.
  4     0.    0. 
  5     0.   50.
  6     0.   70.
  7    20.   70.
  8    20.   30.
  9     0.   30.
 10     0.   60. 
 11    10.   60.
 12    10.   40.
 13     0.   40.

  7
  1    4   6  1  2  7        1
  2    4   7  2  3  8        1
  3    4   8  3  4  9        1
  4    4  10  6  7 11        2
  5    4  11  7  8 12        2
  6    4  12  8  9 13        2
  7    4  10 11 12 13        3

*DANBLOCKS
 1    2  8 1    !- 8nq's  - outer ring
 4   4*1. 
 1   1*1.

 2    2  8 1    !- 8nq's  - inner ring
 4   4*1. 
 1   1*1.

 3    2  8 1    !- 8nq's  - core
 4   4*1. 
 4   4*1.

*DELETE_MATERIALS
 -1                !- the orphan ones

*CIRCLE_A_SQUARE
  0.  50.    20.    22.      !- = Xc,Yc, R-old, R-new

*DANBLOCKS
 1    2  8 1    !- 8nq's  - outer ring
 1   1*1. 
 3   1. 2. 4.

 2    2  8 1    !- 8nq's  - inner ring
 1   1*1. 
 5   1. 3. 5. 3. 1.

 3    2  8 1    !- 8nq's  - core
 1   1*1. 
 1   1*1.

*DELETE_MATERIALS
 -1                !- the orphan ones

*DELETE_ORPHAN_NODES
c*X_TO_Y
*WRITE_DANPLOT



      library 'c:\users\cbs\library\cgraphx.obj'
      library 'c:\users\cbs\library\carlib2.obj'
      library 'c:\users\cbs\library\carlib3.obj'
      options(intl,dreal)
      PROGRAM TIRGEO2

C     This program will produce the geometry data for the Thameside 
C     Industrial Route embankment (Dartford Road Embankment)
C     Revision of the previous TIRGEO.FOR

      PARAMETER(IT=2,NOD=8,INEL=500,IX=20,IY=10,INF=5000,nodoft=3)
      PARAMETER(IMAT=500,ISEQ=500)
      INTEGER NUMS(NOD,INEL),NUMS1(NOD,INEL),NF(NODOFT,INF)
      INTEGER NMAT(IMAT),NSEQ(ISEQ)
      REAL XX(IX),YY(IY),CORNOD(IT,INF),CORNO1(IT,INF),COORD(NOD,IT)
      real coordt(it,nod)
      integer*2 k2

      open(unit=7,file='tirgeo2.dat')

C     (Ground Layers)
      NELX=13
      NELY=8
      NN=0
      NEL=0

      XX(1)  =  0.0
      XX(2)  =  1.1
      XX(3)  =  3.4
      XX(4)  =  5.7
      XX(5)  =  8.0
      XX(6)  = 10.3
      XX(7)  = 12.6
      XX(8)  = 14.9
      XX(9)  = 17.0
      XX(10) = 19.1
      XX(11) = 21.5
      XX(12) = 24.5
      XX(13) = 27.5
      XX(14) = 31.0


      YY(1)  = -7.6
      YY(2)  = -7.1
      YY(3)  = -5.6
      YY(4)  = -4.0
      YY(5)  = -2.5
      YY(6)  = -1.0
      YY(7)  =  0.4
      YY(8)  =  1.0
      YY(9)  =  1.4

      CALL BLOCK1(NELX,NELY,XX,YY,NUMS1,CORNO1,NN1,NEL1)
      CALL ADDBLK(NUMS,NUMS1,CORNOD,CORNO1,NN,NEL,NN1,NEL1,IT,NOD)


C     (40 cm Drainage Layer)
      NELX=9
      NELY=1

      XX(1)  =  0.0
      XX(2)  =  1.1
      XX(3)  =  3.4
      XX(4)  =  5.7
      XX(5)  =  8.0
      XX(6)  = 10.3
      XX(7)  = 12.6
      XX(8)  = 14.9
      XX(9)  = 17.0
      XX(10) = 19.1

      YY(1)  =  1.4
      YY(2)  =  1.8

      CALL BLOCK1(NELX,NELY,XX,YY,NUMS1,CORNO1,NN1,NEL1)
      CORNO1(1,29)=18.7
      CORNO1(1,48)=18.3
      CORNO1(1,47)=(CORNO1(1,48)+CORNO1(1,46))*0.5
      CALL ADDBLK(NUMS,NUMS1,CORNOD,CORNO1,NN,NEL,NN1,NEL1,IT,NOD)



C     (Geocell Block)
      NELX=7
      NELY=1
      
      XX(1)  =  0.0
      XX(2)  =  1.1
      XX(3)  =  3.4
      XX(4)  =  5.7
      XX(5)  =  8.0
      XX(6)  = 10.3
      XX(7)  = 12.6
      XX(8)  = 14.9

      YY(1)  =  1.8
      YY(2)  =  2.8

      CALL BLOCK1(NELX,NELY,XX,YY,NUMS1,CORNO1,NN1,NEL1)
c      DX=0.2
c      DY=0.2
c      CALL ENVELO(NELX,NELY,DX,DY,CORNO1,NUMS1,NEL1,NN1)
      CALL ADDBLK(NUMS,NUMS1,CORNOD,CORNO1,NN,NEL,NN1,NEL1,IT,NOD)


C     (Toe end of the geocell layer)
      NELX=2
      NELY=1

      XX(1)  = 14.9
      XX(2)  = 17.0
      XX(3)  = 18.3

      YY(1)  = 1.8
      YY(2)  = 2.8

      CALL BLOCK1(NELX,NELY,XX,YY,NUMS1,CORNO1,NN1,NEL1)
      CORNO1(1,11)=16.3
      CORNO1(1,13)=(CORNO1(1,5)+CORNO1(1,11))*0.5
      CORNO1(2,13)=(CORNO1(2,5)+CORNO1(2,11))*0.5
      CORNO1(1,12)=(CORNO1(1,13)+CORNO1(1,11))*0.5
      CORNO1(2,12)=(CORNO1(2,13)+CORNO1(2,11))*0.5
      CORNO1(1,8)=(CORNO1(1,5)+CORNO1(1,13))*0.5
      CORNO1(2,8)=(CORNO1(2,5)+CORNO1(2,13))*0.5
      CORNO1(1,7)=(CORNO1(1,3)+CORNO1(1,11))*0.5
      CORNO1(1,10)=(CORNO1(1,11)+CORNO1(1,9))*0.5
      CALL ADDBLK(NUMS,NUMS1,CORNOD,CORNO1,NN,NEL,NN1,NEL1,IT,NOD)


C     (20 cm Sand layer on top of Geocell)

      NELX=8
      NELY=1

      XX(1)  =  0.0
      XX(2)  =  1.1
      XX(3)  =  3.4
      XX(4)  =  5.7
      XX(5)  =  8.0
      XX(6)  = 10.3
      XX(7)  = 12.6
      XX(8)  = 14.9
      XX(9)  = 16.3

      YY(1)  =  2.8
      YY(2)  =  3.0

      CALL BLOCK1(NELX,NELY,XX,YY,NUMS1,CORNO1,NN1,NEL1)
      CORNO1(1,26)=16.1
      CORNO1(1,43)=15.9
      CORNO1(1,42)=(CORNO1(1,43)+CORNO1(1,41))*0.5
      CALL ADDBLK(NUMS,NUMS1,CORNOD,CORNO1,NN,NEL,NN1,NEL1,IT,NOD)



C     (1.2 m Embankment Fill - Right Carriage Way, 1 layer)
      NELX=7
      NELY=1

      XX(1)  =  1.1
      XX(2)  =  3.4
      XX(3)  =  5.7
      XX(4)  =  8.0
      XX(5)  = 10.3
      XX(6)  = 12.6
      XX(7)  = 14.9
      XX(8)  = 15.9

      YY(1)  =  3.0
      YY(2)  =  4.2

      CALL BLOCK1(NELX,NELY,XX,YY,NUMS1,CORNO1,NN1,NEL1)
      CORNO1(1,26)=3.5                               ! LEFT SLOPE
      CORNO1(1,24)=CORNO1(1,26)                      !     .
      CORNO1(1,27)=(CORNO1(1,26)+CORNO1(1,28))*0.5   !     .
      CORNO1(1,25)=(CORNO1(1,24)+CORNO1(1,26))*0.5   !     .
      CORNO1(1,16)=(CORNO1(1,1)+CORNO1(1,24))*0.5    !     .
      CORNO1(1,17)=(CORNO1(1,3)+CORNO1(1,26))*0.5    !     .

      CORNO1(1,36)=13.5                              ! RIGHT SLOPE
      CORNO1(1,38)=(CORNO1(1,36)+CORNO1(1,15))*0.5   !     .
      CORNO1(1,37)=(CORNO1(1,38)+CORNO1(1,36))*0.5   !     .
      CORNO1(1,23)=(CORNO1(1,38)+CORNO1(1,15))*0.5   !     .
      CORNO1(1,22)=(CORNO1(1,36)+CORNO1(1,13))*0.5   !     .
      CORNO1(2,38)=(CORNO1(2,36)+CORNO1(2,15))*0.5   !     .
      CORNO1(2,37)=(CORNO1(2,38)+CORNO1(2,36))*0.5   !     .
      CORNO1(2,23)=(CORNO1(2,38)+CORNO1(2,15))*0.5   !     .
      CORNO1(1,35)=(CORNO1(1,36)+CORNO1(1,34))*0.5   !     .

      CALL ADDBLK(NUMS,NUMS1,CORNOD,CORNO1,NN,NEL,NN1,NEL1,IT,NOD)



C     (1.2 m Embankment Fill between the Left and Right Carriage Way)
      NELX=1
      NELY=1

      XX(1)  =  0.0
      XX(2)  =  1.1

      YY(1)=3.0
      YY(2)=4.2

      CALL BLOCK1(NELX,NELY,XX,YY,NUMS1,CORNO1,NN1,NEL1)
      CORNO1(1,8)=3.5
      CORNO1(1,7)=(CORNO1(1,8)+CORNO1(1,6))*0.5
      CORNO1(1,5)=(CORNO1(1,3)+CORNO1(1,8))*0.5

      CALL ADDBLK(NUMS,NUMS1,CORNOD,CORNO1,NN,NEL,NN1,NEL1,IT,NOD)



C     (1.5 m Surcharge Load -first half :0.8 m)
      NELX=6
      NELY=1

      XX(1)  =  0.0
      XX(2)  =  3.5
      XX(3)  =  5.7
      XX(4)  =  8.0
      XX(5)  = 10.3
      XX(6)  = 12.6
      XX(7)  = 13.5

      YY(1)=4.2
      YY(2)=5.0

      CALL BLOCK1(NELX,NELY,XX,YY,NUMS1,CORNO1,NN1,NEL1)
      CORNO1(1,31)=11.9                              ! RIGHT SLOPE
      CORNO1(1,33)=(CORNO1(1,31)+CORNO1(1,13))*0.5   !     .
      CORNO1(1,32)=(CORNO1(1,33)+CORNO1(1,31))*0.5   !     .
      CORNO1(1,20)=(CORNO1(1,33)+CORNO1(1,13))*0.5   !     .
      CORNO1(1,19)=(CORNO1(1,31)+CORNO1(1,11))*0.5   !     .
      CORNO1(2,33)=(CORNO1(2,31)+CORNO1(2,13))*0.5   !     .
      CORNO1(2,32)=(CORNO1(2,33)+CORNO1(2,31))*0.5   !     .
      CORNO1(2,20)=(CORNO1(2,33)+CORNO1(2,13))*0.5   !     .
      CORNO1(1,30)=(CORNO1(1,31)+CORNO1(1,29))*0.5   !     .
      CALL ADDBLK(NUMS,NUMS1,CORNOD,CORNO1,NN,NEL,NN1,NEL1,IT,NOD)


C     (1.5 m Surcharge Load -final half :0.7 m)
      NELX=5
      NELY=1

      XX(1)  =  0.0
      XX(2)  =  3.5
      XX(3)  =  5.7
      XX(4)  =  8.0
      XX(5)  = 10.3
      XX(6)  = 11.9

      YY(1)=5.0
      YY(2)=5.7

      CALL BLOCK1(NELX,NELY,XX,YY,NUMS1,CORNO1,NN1,NEL1)
      CORNO1(1,28)=10.5                              ! RIGHT SLOPE
      CORNO1(1,27)=(CORNO1(1,28)+CORNO1(1,26))*0.5   !     .
      CORNO1(1,17)=(CORNO1(1,28)+CORNO1(1,11))*0.5   !     .
      CALL ADDBLK(NUMS,NUMS1,CORNOD,CORNO1,NN,NEL,NN1,NEL1,IT,NOD)



C     (Introduce Restrains)
      do i=1,nn      ! All free except Pore Pressure
        nf(1,i)=1    !
        nf(2,i)=1    !
        nf(3,i)=0    !
      end do         !

      do i=1,nel
        nf(3,nums(1,i))=1 ! Free Pore Pressure
        nf(3,nums(3,i))=1 !
        nf(3,nums(5,i))=1 !
        nf(3,nums(7,i))=1 !
      end do

      eps=0.001 ! Tolerance
      wte=1.4   ! Water Table Elevation

      do i=1,nn
        if (cornod(1,i).lt.0.0 +eps) nf(1,i)=0 ! Left Rollers
        if (cornod(1,i).gt.31.0-eps) nf(1,i)=0 ! Right Rollers
        if (cornod(2,i).lt.-7.6+eps) nf(1,i)=0 ! Bottom Fixed
        if (cornod(2,i).lt.-7.6+eps) nf(2,i)=0 !     ..
        if (cornod(2,i).lt.-7.6+eps) nf(3,i)=0 ! Bottom Drained
        if (cornod(2,i).gt. wte-eps) nf(3,i)=0 ! Above WT Drained
      end do
      
C     (Forming NF)
      k=0
      do i=1,nn
        do j=1,nodoft
          if (nf(j,i).gt.0) then
            k=k+1
            nf(j,i)=k
          end if
        end do
      end do

C     (Assign Material Type)
      DO I=1,NEL
        IF (CORNOD(2,(NUMS(3,I))).LT.-5.6+EPS) THEN
          IF (CORNOD(1,(NUMS(5,I))).LT. 19.1+EPS) THEN
            NMAT(I)=1
          ELSE
            NMAT(I)=2
          END IF
        ELSE IF (CORNOD(2,(NUMS(3,I))).LT.-1.0+EPS) THEN
          IF (CORNOD(1,(NUMS(5,I))).LT. 19.1+EPS) THEN
            NMAT(I)=3
          ELSE
            NMAT(I)=4
          END IF
        ELSE IF (CORNOD(2,(NUMS(3,I))).LT. 1.0+EPS) THEN
          IF (CORNOD(1,(NUMS(5,I))).LT. 19.1+EPS) THEN
            NMAT(I)=5
          ELSE
            NMAT(I)=6
          END IF
        ELSE IF (CORNOD(2,(NUMS(3,I))).LT. 1.4+EPS) THEN
          IF (CORNOD(1,(NUMS(5,I))).LT. 19.1+EPS) THEN
            NMAT(I)=7
          ELSE
            NMAT(I)=8
          END IF
        ELSE IF (CORNOD(2,(NUMS(1,I))).GT. 1.8 -EPS .AND.
     +           CORNOD(2,(NUMS(3,I))).LT. 2.8 +EPS .AND.
     +           CORNOD(1,(NUMS(5,I))).LT. 14.9+EPS) THEN
          NMAT(I)=9
        ELSE
          NMAT(I)=10
        END IF
      END DO

C     (Assign Construction Sequence)
      DO I=1,NEL
        IF (CORNOD(2,(NUMS(3,I))).LT. 1.4+EPS) THEN
          NSEQ(I)=0
        ELSE IF (CORNOD(2,(NUMS(3,I))).LT. 1.8+EPS) THEN
          NSEQ(I)=1
        ELSE IF (CORNOD(2,(NUMS(3,I))).LT. 2.8+EPS) THEN
          NSEQ(I)=2
        ELSE IF (CORNOD(2,(NUMS(3,I))).LT. 3.0+EPS) THEN
          NSEQ(I)=3
        ELSE IF (CORNOD(2,(NUMS(3,I))).LT. 4.2+EPS) THEN
          NSEQ(I)=4
          IF (CORNOD(1,(NUMS(3,I))).LT. 0.0+EPS) NSEQ(I)=5
        ELSE IF (CORNOD(2,(NUMS(3,I))).LT. 5.0+EPS) THEN
          NSEQ(I)=6
        ELSE
          NSEQ(I)=7
        END IF
      END DO


C     (Write Title)
      WRITE(7,*)'***************************************************'
      WRITE(7,*)'*         GEOMETRY OF DARTFORD EMBANKMENT         *'
      WRITE(7,*)'***************************************************'
      WRITE(7,*)'~~~~~~~~~~~~~~~~~~ ELEMENT DATA ~~~~~~~~~~~~~~~~~~~'

      do i=1,nel
        write(7,'(11i6)')i,(nums(j,i),j=1,nod),NMAT(I),NSEQ(I)
      end do
      WRITE(7,*)'~~~~~~~~~~~~~ NODAL AND FREEDOM DATA ~~~~~~~~~~~~~~~~~'
      do i=1,nn
        write(7,'(4i5,2f15.5)')i,nf(1,i),nf(2,i),nf(3,i),
     +                         cornod(1,i),cornod(2,i)
      end do

      CALL DEVICE_VGA
      call pspace(0.,1.,0.,1.)
      call map(-5.,35.,-15.,25.)

      DO J=1,NEL
        CALL ASCORD(J,NUMS,NOD,CORNOD,IT,COORD,NOD)
        CALL MATRAN(COORDT,IT,COORD,NOD,NOD,IT)
        CALL DREL2D(COORDT,NOD,NMAT(J),4)
      END DO

      END

c      INCLUDE 'CGRAPHX.FOR'

      SUBROUTINE BLOCK1(NELX,NELY,XX,YY,NUMS,CORNOD,NN,NEL)
C     This routine creates a block of mesh NELX by NELY 8-noded elements
C     XX is the X-coords and YY is the Y-coords of the grids. The routi-
C     ne returns the element conectivity NUMS and the nodal coord CORNOD
C     NN is the total number of the nodes in the block.

      PARAMETER (NOD=8,IT=2)
      REAL    XX(*),YY(*),CORNOD(IT,*)
      INTEGER NUMS(NOD,*)

      NEL=0
      DO I=1,NELY
        DO J=1,NELX
          NEL=NEL+1
          NUMS(1,NEL)=(3*NELX+2)*(I-1)+(J-1)*2+1  
          NUMS(2,NEL)=(3*NELX+2)*(I-1)+NELX*2+1+J 
          NUMS(3,NEL)=(3*NELX+2)*I+(J-1)*2+1      
          NUMS(4,NEL)=NUMS(3,NEL)+1
          NUMS(5,NEL)=NUMS(3,NEL)+2
          NUMS(6,NEL)=NUMS(2,NEL)+1
          NUMS(7,NEL)=NUMS(1,NEL)+2
          NUMS(8,NEL)=NUMS(1,NEL)+1
      
          CORNOD(1,NUMS(1,NEL))=XX(J)
          CORNOD(1,NUMS(2,NEL))=XX(J)
          CORNOD(1,NUMS(3,NEL))=XX(J)
          CORNOD(1,NUMS(5,NEL))=XX(J+1)
          CORNOD(1,NUMS(6,NEL))=XX(J+1)
          CORNOD(1,NUMS(7,NEL))=XX(J+1)
          CORNOD(1,NUMS(4,NEL))=(XX(J)+XX(J+1))/2.
          CORNOD(1,NUMS(8,NEL))=(XX(J)+XX(J+1))/2.

          CORNOD(2,NUMS(1,NEL))=YY(I)
          CORNOD(2,NUMS(7,NEL))=YY(I)
          CORNOD(2,NUMS(8,NEL))=YY(I)
          CORNOD(2,NUMS(3,NEL))=YY(I+1)
          CORNOD(2,NUMS(4,NEL))=YY(I+1)
          CORNOD(2,NUMS(5,NEL))=YY(I+1)
          CORNOD(2,NUMS(2,NEL))=(YY(I)+YY(I+1))/2.
          CORNOD(2,NUMS(6,NEL))=(YY(I)+YY(I+1))/2.
        END DO
      END DO

      NN=(3*NELX+2)*NELY+2*NELX+1
            
      RETURN
      END

      SUBROUTINE BLOCK2(NELX,NELY,QMAP,COORD,NUMS,CORNOD,NN,NEL)
C     This routine creates a block of 8-noded elements NELX by NELY
C     and mapped the mesh using isoparametric shape function linearly
C     (QMAP=.FALSE.) or quadratically (QMAP=.TRUE.). COORD contains the
C     coordinates of the four corners of the block. The routine returns
C     the element connectivity NUMS and nodal coordinates CORNOD. NN is 
C     the total node number in the block.

      PARAMETER(NOD=8,IT=2,NODLIM=10000)
      REAL COORD(NOD,IT),CORNOD(IT,*),FUN(NOD),XIETA(IT,NODLIM)
      INTEGER NUMS(NOD,*)
      LOGICAL QMAP

      IF (QMAP) THEN
        NPT=8
      ELSE
        NPT=4
      END IF
      DXI=2./NELX
      DETA=2./NELY
      DO I=1,NELY
        DO J=1,NELX
          NEL=NEL+1
          NUMS(1,NEL)=(3*NELX+2)*(I-1)+(J-1)*2+1
          NUMS(2,NEL)=(3*NELX+2)*(I-1)+NELX*2+1+J
          NUMS(3,NEL)=(3*NELX+2)*I+(J-1)*2+1
          NUMS(4,NEL)=NUMS(3,NEL)+1
          NUMS(5,NEL)=NUMS(3,NEL)+2
          NUMS(6,NEL)=NUMS(2,NEL)+1
          NUMS(7,NEL)=NUMS(1,NEL)+2
          NUMS(8,NEL)=NUMS(1,NEL)+1

          XIETA(1,NUMS(1,NEL))=-1.+DXI*(J-1)
          XIETA(1,NUMS(2,NEL))=-1.+DXI*(J-1)
          XIETA(1,NUMS(3,NEL))=-1.+DXI*(J-1)
          XIETA(1,NUMS(5,NEL))=-1.+DXI*J
          XIETA(1,NUMS(6,NEL))=-1.+DXI*J
          XIETA(1,NUMS(7,NEL))=-1.+DXI*J
          XIETA(1,NUMS(4,NEL))=-1.+DXI*(J-1)+DXI*0.5
          XIETA(1,NUMS(8,NEL))=-1.+DXI*(J-1)+DXI*0.5

          XIETA(2,NUMS(1,NEL))=-1.+DETA*(I-1)
          XIETA(2,NUMS(7,NEL))=-1.+DETA*(I-1)
          XIETA(2,NUMS(8,NEL))=-1.+DETA*(I-1)
          XIETA(2,NUMS(3,NEL))=-1.+DETA*I
          XIETA(2,NUMS(4,NEL))=-1.+DETA*I
          XIETA(2,NUMS(5,NEL))=-1.+DETA*I
          XIETA(2,NUMS(2,NEL))=-1.+DETA*(I-1)+DETA*0.5
          XIETA(2,NUMS(6,NEL))=-1.+DETA*(I-1)+DETA*0.5

        END DO
      END DO

      NN=NN+(3*NELX+2)*NELY+2*NELX+1

      DO I=1,NN
        IF (QMAP) THEN
          CALL FUN8Q(XIETA(1,I),XIETA(2,I),FUN)
        ELSE
          CALL FUN4Q(XIETA(1,I),XIETA(2,I),FUN)
        END IF
        CORNOD(1,I)=0.0
        CORNOD(2,I)=0.0
        DO J=1,NPT
          CORNOD(1,I)=CORNOD(1,I)+FUN(J)*COORD(J,1)
          CORNOD(2,I)=CORNOD(2,I)+FUN(J)*COORD(J,2)
        END DO
      END DO

      RETURN
      END


      SUBROUTINE ADDBLK(NUMS,NUMS1,CORNOD,CORNO1,NN,NEL,NN1,NEL1,IT,NOD)
C     Adds a block of mesh (NUMS1;CORNO1;NN1;NEL1) into 
C     global block (NUMS;CORNOD;NN;NEL)

      PARAMETER(IMATCH=100,IPAIR=2,IDIM=3,tol=0.0001)
      REAL CORNOD(IT,*),CORNO1(IT,*),DIF(IDIM)
      INTEGER NUMS(NOD,*),NUMS1(NOD,*),MATCH(IPAIR,IMATCH)

C     (ADDING BLOCK)
      DO I=1,NEL1
        NEL=NEL+1
        DO J=1,NOD
          NUMS(J,NEL)=NUMS1(J,I)+NN
        END DO
      END DO

      DO I=1,NN1
        NN=NN+1
        DO J=1,IT
          CORNOD(J,NN)=CORNO1(J,I)
        END DO
      END DO

C     (FIND COINCIDENCE POINTS)
      II=0
      NNLBLK=NN-NN1
      DO I=NNLBLK+1,NN
        DO J=1,NNLBLK
          DO K=1,IT
            DIF(K)=ABS(CORNOD(K,I)-CORNOD(K,J))
          END DO
          IF (DIF(1).LT.TOL .AND. DIF(2).LT.TOL) THEN
            II=II+1
            MATCH(1,II)=J
            MATCH(2,II)=I
c            write(6,*)ii,match(1,ii),match(2,ii)
            GO TO 100
          END IF
        END DO
  100   CONTINUE
      END DO

C     (RENUMBERING)
      NELLOW=NEL-NEL1
      DO IM=1,II

        DO I=NELLOW+1,NEL
          DO J=1,NOD
            IF (NUMS(J,I).EQ.MATCH(2,IM)) NUMS(J,I)=MATCH(1,IM)
            IF (NUMS(J,I).GT.MATCH(2,IM)) NUMS(J,I)=NUMS(J,I)-1
          END DO
        END DO
        
        DO I=IM,II
          IF (MATCH(2,I).GT.MATCH(2,IM)) MATCH(2,I)=MATCH(2,I)-1
        END DO

        DO I=MATCH(2,IM),NN-1
          DO J=1,IT
            CORNOD(J,I)=CORNOD(J,I+1)
          END DO
        END DO
        NN=NN-1

      END DO          

      RETURN
      END


      SUBROUTINE ENVELO(NELX,NELY,DX,DY,CORNOD,NUMS,NEL1,NN1)
C     Creates elements enveloping a block of mesh 

      PARAMETER (NOD=8,IT=2,INN=10000,IBL=1000)
      REAL CORNOD(IT,*),CORNO0(IT,INN),CORNO1(IT,INN),CORNO2(IT,INN),
     +     CORNO3(IT,INN),CORNO4(IT,INN)

      INTEGER NUMS(NOD,*),NUMS0(NOD,IBL),NUMS1(NOD,IBL),NUMS2(NOD,IBL),
     +        NUMS3(NOD,IBL),NUMS4(NOD,IBL)

      NN=(NELX*3+2)*NELY+NELX*2+1
      DO I=1,NN
        DO J=1,IT
          CORNO0(J,I)=CORNOD(J,I)
        END DO
      END DO

C     (ADJUST COORDINATES OF THE MAIN BLOCK)
      NEL=0
      DO NY=1,NELY
        DO NX=1,NELX
          NEL=NEL+1
          DO I=1,NOD
            NUMS0(I,NEL)=NUMS(I,NEL)
          END DO

          IF (NY.EQ.1) THEN
            CORNO0(2,NUMS0(1,NEL))=CORNOD(2,NUMS(1,NEL))+DY  
            CORNO0(2,NUMS0(7,NEL))=CORNOD(2,NUMS(7,NEL))+DY  
            CORNO0(2,NUMS0(8,NEL))=CORNOD(2,NUMS(8,NEL))+DY  
            CORNO0(2,NUMS0(2,NEL))= (CORNO0(2,NUMS0(1,NEL))+
     +                               CORNO0(2,NUMS0(3,NEL)))*0.5
            CORNO0(2,NUMS0(6,NEL))= (CORNO0(2,NUMS0(5,NEL))+
     +                               CORNO0(2,NUMS0(7,NEL)))*0.5
          END IF
          IF (NY.EQ.NELY) THEN
            CORNO0(2,NUMS0(3,NEL))=CORNOD(2,NUMS(3,NEL))-DY  
            CORNO0(2,NUMS0(4,NEL))=CORNOD(2,NUMS(4,NEL))-DY  
            CORNO0(2,NUMS0(5,NEL))=CORNOD(2,NUMS(5,NEL))-DY  
            CORNO0(2,NUMS0(2,NEL))= (CORNO0(2,NUMS0(1,NEL))+
     +                               CORNO0(2,NUMS0(3,NEL)))*0.5
            CORNO0(2,NUMS0(6,NEL))= (CORNO0(2,NUMS0(5,NEL))+
     +                               CORNO0(2,NUMS0(7,NEL)))*0.5
          END IF
          IF (NX.EQ.1) THEN
            CORNO0(1,NUMS0(1,NEL))=CORNOD(1,NUMS(1,NEL))+DX  
            CORNO0(1,NUMS0(2,NEL))=CORNOD(1,NUMS(2,NEL))+DX  
            CORNO0(1,NUMS0(3,NEL))=CORNOD(1,NUMS(3,NEL))+DX  
            CORNO0(1,NUMS0(8,NEL))= (CORNO0(1,NUMS0(1,NEL))+
     +                               CORNO0(1,NUMS0(7,NEL)))*0.5
            CORNO0(1,NUMS0(4,NEL))= (CORNO0(1,NUMS0(3,NEL))+
     +                               CORNO0(1,NUMS0(5,NEL)))*0.5
          END IF
          IF (NX.EQ.NELX) THEN
            CORNO0(1,NUMS0(5,NEL))=CORNOD(1,NUMS(5,NEL))-DX  
            CORNO0(1,NUMS0(6,NEL))=CORNOD(1,NUMS(6,NEL))-DX  
            CORNO0(1,NUMS0(7,NEL))=CORNOD(1,NUMS(7,NEL))-DX  
            CORNO0(1,NUMS0(8,NEL))= (CORNO0(1,NUMS0(1,NEL))+
     +                               CORNO0(1,NUMS0(7,NEL)))*0.5
            CORNO0(1,NUMS0(4,NEL))= (CORNO0(1,NUMS0(3,NEL))+
     +                               CORNO0(1,NUMS0(5,NEL)))*0.5
          END IF
        END DO
      END DO

C     (CREATE THE ENVELOPE AT THE BOTTOM)
      DO I=1,NELX
        NUMS1(1,I)=(I-1)*2+1
        NUMS1(2,I)=NELX*2+1+I
        NUMS1(3,I)=NELX*3+2+(I-1)*2+1
        NUMS1(4,I)=NUMS1(3,I)+1
        NUMS1(5,I)=NUMS1(3,I)+2
        NUMS1(6,I)=NUMS1(2,I)+1
        NUMS1(7,I)=NUMS1(1,I)+2
        NUMS1(8,I)=NUMS1(1,I)+1

        DO J=1,IT
          CORNO1(J,NUMS1(1,I))=CORNOD(J,NUMS(1,I))
          CORNO1(J,NUMS1(8,I))=CORNOD(J,NUMS(8,I))
          CORNO1(J,NUMS1(7,I))=CORNOD(J,NUMS(7,I))
          CORNO1(J,NUMS1(3,I))=CORNO0(J,NUMS0(1,I))
          CORNO1(J,NUMS1(4,I))=CORNO0(J,NUMS0(8,I))
          CORNO1(J,NUMS1(5,I))=CORNO0(J,NUMS0(7,I))
          CORNO1(J,NUMS1(2,I))= (CORNO1(J,NUMS1(1,I))+
     +                           CORNO1(J,NUMS1(3,I)))*0.5
          CORNO1(J,NUMS1(6,I))= (CORNO1(J,NUMS1(5,I))+
     +                           CORNO1(J,NUMS1(7,I)))*0.5
        END DO
      END DO

C     (CREATE THE ENVELOPE AT THE LEFT EDGE)
      DO I=1,NELY
        NUMS2(1,I)=5*(I-1)+1
        NUMS2(2,I)=NUMS2(1,I)+3
        NUMS2(3,I)=NUMS2(1,I)+5
        NUMS2(4,I)=NUMS2(3,I)+1
        NUMS2(5,I)=NUMS2(3,I)+2
        NUMS2(6,I)=NUMS2(2,I)+1
        NUMS2(7,I)=NUMS2(1,I)+2
        NUMS2(8,I)=NUMS2(1,I)+1

        NXY=(I-1)*NELX+1
        DO J=1,IT
          CORNO2(J,NUMS2(1,I))=CORNOD(J,NUMS(1,NXY))
          CORNO2(J,NUMS2(2,I))=CORNOD(J,NUMS(2,NXY))
          CORNO2(J,NUMS2(3,I))=CORNOD(J,NUMS(3,NXY))
          CORNO2(J,NUMS2(7,I))=CORNO0(J,NUMS0(1,NXY))
          CORNO2(J,NUMS2(6,I))=CORNO0(J,NUMS0(2,NXY))
          CORNO2(J,NUMS2(5,I))=CORNO0(J,NUMS0(3,NXY))
          CORNO2(J,NUMS2(4,I))= (CORNO2(J,NUMS2(3,I))+
     +                           CORNO2(J,NUMS2(5,I)))*0.5
          CORNO2(J,NUMS2(8,I))= (CORNO2(J,NUMS2(1,I))+
     +                           CORNO2(J,NUMS2(7,I)))*0.5
        END DO
      END DO

C     (CREATE THE ENVELOPE AT THE RIGHT EDGE)
      DO I=1,NELY
        NUMS3(1,I)=5*(I-1)+1
        NUMS3(2,I)=NUMS3(1,I)+3
        NUMS3(3,I)=NUMS3(1,I)+5
        NUMS3(4,I)=NUMS3(3,I)+1
        NUMS3(5,I)=NUMS3(3,I)+2
        NUMS3(6,I)=NUMS3(2,I)+1
        NUMS3(7,I)=NUMS3(1,I)+2
        NUMS3(8,I)=NUMS3(1,I)+1

        NXY=NELX*I
        DO J=1,IT
          CORNO3(J,NUMS3(3,I))=CORNO0(J,NUMS0(5,NXY))
          CORNO3(J,NUMS3(2,I))=CORNO0(J,NUMS0(6,NXY))
          CORNO3(J,NUMS3(1,I))=CORNO0(J,NUMS0(7,NXY))
          CORNO3(J,NUMS3(5,I))=CORNOD(J,NUMS(5,NXY))
          CORNO3(J,NUMS3(6,I))=CORNOD(J,NUMS(6,NXY))
          CORNO3(J,NUMS3(7,I))=CORNOD(J,NUMS(7,NXY))
          CORNO3(J,NUMS3(4,I))= (CORNO3(J,NUMS3(3,I))+
     +                           CORNO3(J,NUMS3(5,I)))*0.5
          CORNO3(J,NUMS3(8,I))= (CORNO3(J,NUMS3(1,I))+
     +                           CORNO3(J,NUMS3(7,I)))*0.5
        END DO
      END DO

C     (CREATE THE ENVELOPE AT THE TOP)
      NELXY=NELX*(NELY-1)
      DO I=1,NELX
        NUMS4(1,I)=(I-1)*2+1
        NUMS4(2,I)=NELX*2+1+I
        NUMS4(3,I)=NELX*3+2+(I-1)*2+1
        NUMS4(4,I)=NUMS4(3,I)+1
        NUMS4(5,I)=NUMS4(3,I)+2
        NUMS4(6,I)=NUMS4(2,I)+1
        NUMS4(7,I)=NUMS4(1,I)+2
        NUMS4(8,I)=NUMS4(1,I)+1

        DO J=1,IT
          CORNO4(J,NUMS4(1,I))=CORNO0(J,NUMS0(3,I+NELXY))
          CORNO4(J,NUMS4(8,I))=CORNO0(J,NUMS0(4,I+NELXY))
          CORNO4(J,NUMS4(7,I))=CORNO0(J,NUMS0(5,I+NELXY))
          CORNO4(J,NUMS4(3,I))=CORNOD(J,NUMS(3,I+NELXY))
          CORNO4(J,NUMS4(4,I))=CORNOD(J,NUMS(4,I+NELXY))
          CORNO4(J,NUMS4(5,I))=CORNOD(J,NUMS(5,I+NELXY))
          CORNO4(J,NUMS4(2,I))= (CORNO4(J,NUMS4(1,I))+
     +                           CORNO4(J,NUMS4(3,I)))*0.5
          CORNO4(J,NUMS4(6,I))= (CORNO4(J,NUMS4(5,I))+
     +                           CORNO4(J,NUMS4(7,I)))*0.5
        END DO
      END DO

      NN=0
      NEL=0
      NN1=NELX*3+2+NELX*2+1
      NEL1=NELX
      CALL ADDBLK(NUMS,NUMS1,CORNOD,CORNO1,NN,NEL,NN1,NEL1,IT,NOD)
      NN1=(NELX*3+2)*NELY+NELX*2+1
      NEL1=NELX*NELY
      CALL ADDBLK(NUMS,NUMS0,CORNOD,CORNO0,NN,NEL,NN1,NEL1,IT,NOD)
      NN1=NELY*5+3
      NEL1=NELY
      CALL ADDBLK(NUMS,NUMS2,CORNOD,CORNO2,NN,NEL,NN1,NEL1,IT,NOD)
      NN1=NELY*5+3
      NEL1=NELY
      CALL ADDBLK(NUMS,NUMS3,CORNOD,CORNO3,NN,NEL,NN1,NEL1,IT,NOD)
      NN1=NELX*5+3
      NEL1=NELX
      CALL ADDBLK(NUMS,NUMS4,CORNOD,CORNO4,NN,NEL,NN1,NEL1,IT,NOD)
      NEL1=NEL
      NN1=NN

      RETURN
      END



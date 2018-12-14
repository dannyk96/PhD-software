      library 'c:\users\cbs\library\cgraphx.obj'
      library 'c:\users\cbs\library\carlib2.obj'
      options (dreal,intl,check,undef)
      PROGRAM TBLOCK

C     Program for testing random generation routines

      PARAMETER (IT=2,NOD=8,NODOFT=3,INF=500,INEL=300)
      CHARACTER FILGEO*60
      INTEGER   NUMS(NOD,INEL),NF(NODOFT,INF),NOMAT(INF),NOSTG(INF)
      REAL      CORNOD(IT,INF),COORD(NOD,IT),COORDT(IT,NOD)
      INTEGER*2 K2
      
      FILGEO='c:\users\cbs\prog\cubzac\geometry.dat'
      CALL RDGEO_(FILGEO,NNEL,NN,NUMS,NOD,NF,NODOFT,CORNOD,IT,NOMAT,
     +            NOSTG)
      
      XMAX=-1.E20
      XMIN= 1.E20
      YMAX=-1.E20
      YMIN= 1.E20

      DO I=1,NN
        XMAX=MAX(XMAX,CORNOD(1,I))
        XMIN=MIN(XMIN,CORNOD(1,I))
        YMAX=MAX(YMAX,CORNOD(2,I))
        YMIN=MIN(YMIN,CORNOD(2,I))
      END DO
      

      call device_vga
      call pspace(0.0,1.0,0.0,1.0)
      range=max((xmax-xmin),(ymax-ymin))
      xmid=(xmax+xmin)/2.
      ymid=(ymax+ymin)/2.
      dpic=range*0.55
      call map(xmid-dpic,xmid+dpic,ymid-dpic,ymid+dpic)

      do nel=1,nnel
        call ascord(nel,nums,nod,cornod,it,coord,nod)
        do i=1,nod
          do j=1,it
            coordt(j,i)=coord(i,j)
          end do
        end do
        call drel2d(coordt,nod,3s,4s)
      end do

      CALL MATFRM(2,NOmat,NNEL,NUMS,NOD,cornod,10)
      
c      if (npt.gt.0) then
c        do i=1,npt
c          i2x(i)=n2ofx(cornod(1,nodfrm(i)))
c          i2y(i)=n2ofy(cornod(2,nodfrm(i)))
c        end do
c        call polyline@(i2x,i2y,ints(npt),10s)
c      end if
      
      call get_key@(k2)
      CALL TEXT_MODE@
      END



      SUBROUTINE MATFRM(ICRIT,MTCRIT,NNEL,NUMS,NOD,CORNOD,ICOL)

      PARAMETER  (NPAIR=2000,ifrpt=500,it=2)
      real       cornod(it,*)
      INTEGER    NUMS(NOD,*),MTCRIT(*),NODSEG(2,NPAIR),NODFRM(ifrpt)
      integer*2  n2ofx,n2ofy,i2x(ifrpt),i2y(ifrpt)


      L=0
      DO I=1,NNEL
        IF (MTCRIT(I).EQ.ICRIT) THEN
          DO J=1,NOD
            K=MOD(J,NOD)+1
            L=L+1
            NODSEG(1,L)=MIN(NUMS(J,I),NUMS(K,I))
            NODSEG(2,L)=MAX(NUMS(J,I),NUMS(K,I))
          END DO
        END IF
      END DO

      IF (L.NE.0) THEN
        LPAIR=L
        DO I=1,LPAIR-1
          IF (NODSEG(1,I).NE.0) THEN
            DO J=I+1,LPAIR 
              IF (NODSEG(1,I).EQ.NODSEG(1,J) .AND.
     +            NODSEG(2,I).EQ.NODSEG(2,J) )THEN
                NODSEG(1,I)=0
                NODSEG(1,J)=0
                NODSEG(2,I)=0
                NODSEG(2,J)=0
                GO TO 50
              END IF
            END DO
          END IF
   50     CONTINUE
        END DO
        
        NODFRM(1)=NODSEG(1,1)
        NODFRM(2)=NODSEG(2,1)
        NODSEG(2,1)=0
        NODEND=NODFRM(2)
        K=2

   60   CONTINUE
        DO J=2,LPAIR
          IF (NODSEG(1,J).NE.0) THEN
            IF (NODSEG(1,J).EQ.NODEND .OR. NODSEG(2,J).EQ.NODEND) THEN
              K=K+1
              IF (NODSEG(1,J).EQ.NODEND) THEN
                NODFRM(K)=NODSEG(2,J)
              ELSE
                NODFRM(K)=NODSEG(1,J)
              END IF
              NODSEG(1,J)=0
              NODSEG(2,J)=0
              NODEND=NODFRM(K)
              GO TO 60
            END IF
          END IF
        END DO
      
        DO I=1,K
          I2X(I)=N2OFX(CORNOD(1,NODFRM(I)))
          I2Y(I)=N2OFY(CORNOD(2,NODFRM(I)))
        END DO
        CALL POLYLINE@(I2X,I2Y,INTS(K),INTS(ICOL))
      END IF

      RETURN
      END




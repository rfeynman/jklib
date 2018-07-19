	program ts

      real*8 a(5,5), b(5,5), c(5,5), det
      integer*4 LTMP(5),MTMP(5)

      do 1 i=1,5
      LTMP(i)=1
      MTMP(i)=1
      do 1 j=1,5
      a(i,j)= i/3.+j/5. + i*j/2.
   1  continue
c     do 2 i=1,5
c     a(i,i)=i
c  2  continue

      do 5 i=1,5
      write(*, 100) (a(i,j), j=1,5)
   5  continue   
      write(*, 101) 

      do 6 i=1,5
      do 6 j=1,5
      b(i,j)=a(i,j)
   6  continue

      do 7 i=1,5
      write(*, 100) (b(i,j), j=1,5)
   7  continue   
      write(*, 101) 

      CALL DMINV(A,5,DET,LTMP,MTMP)


      do 8 i=1,5
      write(*, 100) (a(i,j), j=1,5)
   8  continue   
      write(*, 101) 

	
      do 3 i=1,5
      do 3 j=1,5
      c(i,j)=0.d0
      do 3 k=1,5
      c(i,j)=c(i,j)+a(i,k)*b(k,j)
   3  continue

      do 4 i=1,5
      write(*, 100) (c(i,j), j=1,5)
   4  continue   
 100  format(1x, 1p5d17.8)
 101  format(1x)
      stop
      end


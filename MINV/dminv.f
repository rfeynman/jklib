C   14/11/86 611282107  MEMBER NAME  DMINV    (PETLIB)      FORTRAN77
C
C     ..................................................................
C
C        SUBROUTINE DMINV
C
C        PURPOSE
C           INVERT A MATRIX
C
C        USAGE
C           CALL DMINV(A,N,D,L,M)
C
C        DESCRIPTION OF PARAMETERS
C           A - INPUT MATRIX, DESTROYED IN COMPUTATION AND REPLACED BY
C               RESULTANT INVERSE.
C           N - ORDER OF MATRIX A
C           D - RESULTANT DETERMINANT
C           L - WORK VECTOR OF LENGTH N
C           M - WORK VECTOR OF LENGTH N
C
C        REMARKS
C           MATRIX A MUST BE A GENERAL MATRIX
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           NONE
C
C        METHOD
C           THE STANDARD GAUSS-JORDAN METHOD IS USED. THE DETERMINANT
C           IS ALSO CALCULATED. A DETERMINANT OF ZERO INDICATES THAT
C           THE MATRIX IS SINGULAR.
C
C        TRANSLATED INTO FORTRAN77 11.86 JKL
C     ..................................................................
C
      SUBROUTINE DMINV(A,N,D,L,M)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),L(1),M(1)
      write(*,fmt='(1x, 25d10.3)') (a(i), i=1,N*N)
      D=1.D0
      NK=-N
      DO 80 K=1,N
c write(*, fmt='(4hk   ,i8)') k 
         NK=NK+N
         L(K)=K
         M(K)=K
         KK=NK+K
c write(*, fmt='(4hkk  ,i8)') kk
         BIGA=A(KK)
c      write(*,101) biga, kk
c 101	format(1x,d15.8, i8)
         DO 21 J=K,N
c write(*, fmt='(4hj   ,i8)') j 
            IZ=N*(J-1)
c write(*, fmt='(4hiz  ,i8)') iz
            DO 20 I=K,N
               IJ=IZ+I
c write(*, fmt='(4hij  ,i8)') ij
               IF(DABS(BIGA).LT.DABS(A(IJ))) THEN
                  BIGA=A(IJ)
                  L(K)=I
                  M(K)=J
               ENDIF
   20       CONTINUE
   21    CONTINUE
c       write(*, fmt='(1x, 5i8)') (l(is), is=1,5)
c       write(*, fmt='(1x, 5i8)') (m(is), is=1,5)
         J=L(K)
         IF(J.GT.K) THEN
            KI=K-N
            DO 30 I=1,N
               KI=KI+N
c write(*, fmt='(4hki  ,i8)') ki
               HOLD=-A(KI)
c		       write(*,102) hold
c 102	format(1x,'h1 ', d15.8)
               JI=KI-K+J
c write(*, fmt='(4hji  ,i8)') ji
               A(KI)=A(JI)
               A(JI) =HOLD
   30       CONTINUE
         ENDIF
         I=M(K)
         IF(I.GT.K) THEN
            JP=N*(I-1)
            DO 40 J=1,N
               JK=NK+J
c write(*, fmt='(4hjk  ,i8)') jk
               JI=JP+J
c write(*, fmt='(4hji  ,i8)') ji
               HOLD=-A(JK)
               A(JK)=A(JI)
               A(JI) =HOLD
   40       CONTINUE
         ENDIF
         IF(BIGA.EQ. 0.D0) THEN
            D=0D0
            RETURN
         ENDIF
         DO 55 I=1,N
            IF(I.NE.K) THEN
               IK=NK+I
c write(*, fmt='(4hik  ,i8)') ik
               A(IK)=A(IK)/(-BIGA)
            ENDIF
   55    CONTINUE
         DO 66 I=1,N
            IK=NK+I
c write(*, fmt='(4hik  ,i8)') ik
            HOLD=A(IK)
            IJ=I-N
c write(*, fmt='(4hij  ,i8)') ij
            DO 65 J=1,N
               IJ=IJ+N
c write(*, fmt='(4hij  ,i8)') ij
               IF(I.NE.K.AND.J.NE.K) THEN
                  KJ=IJ-I+K
c write(*, fmt='(4hkj  ,i8)') kj
                  A(IJ)=HOLD*A(KJ)+A(IJ)
               ENDIF
   65       CONTINUE
   66    CONTINUE
         KJ=K-N
c write(*, fmt='(4hkj  ,i8)') kj
         DO 75 J=1,N
            KJ=KJ+N
c write(*, fmt='(4hkj  ,i8)') kj
            IF(J.NE.K) A(KJ)=A(KJ)/BIGA
   75    CONTINUE
         D=D*BIGA
c write(*, fmt='(4hkk  ,i8)') kk
         A(KK)=1D0/BIGA
   80 CONTINUE
      DO 90 K=N-1,1,-1
         I=L(K)
         IF(I.GT.K) THEN
            JQ=N*(K-1)
            JR=N*(I-1)
            DO 110 J=1,N
               JK=JQ+J
c write(*, fmt='(4hjk  ,i8)') jk
               HOLD=A(JK)
               JI=JR+J
c write(*, fmt='(4hji  ,i8)') ji
               A(JK)=-A(JI)
               A(JI) =HOLD
  110       CONTINUE
         ENDIF
         J=M(K)
         IF(J.GT.K) THEN
            KI=K-N
            DO 130 I=1,N
               KI=KI+N
c write(*, fmt='(4hki  ,i8)') ki
               HOLD=A(KI)
               JI=KI-K+J
c write(*, fmt='(4hji  ,i8)') ji
               A(KI)=-A(JI)
               A(JI) =HOLD
  130       CONTINUE
         ENDIF
   90 CONTINUE
      RETURN
      END

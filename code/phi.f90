!!!! This program creates a file phi.dat containing the matrix phi 
!!!! for excluded volume potential
       PROGRAM CALC
       IMPLICIT NONE
       Integer, Parameter :: k4b = Selected_int_kind(9)
       Integer, Parameter :: SNGL = SELECTED_REAL_KIND(4)
       INTEGER, PARAMETER :: DOBL = SELECTED_REAL_KIND(8)
       INTEGER, PARAMETER :: FPREC = DOBL
       integer :: i, j, k, l
       REAL (FPREC) :: X, A, B, NUM, DEN, ANS1, c, num1, num2
       !real(8),  parameter :: PI  = 4 * atan (1.0_8)
       INTEGER, PARAMETER :: N = 250
       real(FPREC) :: xx(N,N), YY(N,N)
       Do i =1, N
         DO j = 1,N
            XX(i,j) = 0.8d0  ! write value of phi, which u want
         END DO
       END DO
       OPEN (UNIT=11,file="phi.dat",status="unknown")
       Do k =1, N
        ! DO j = 1,4
            write(11,*) XX(k,:)
         !END DO
       END DO
      close (11)
      !OPEN (UNIT=12,FILE="phi.dat", status="old")
     ! Do l = 1,N
       ! read(12,*) YY (l,:) 
      !End DO 
      !close (12)
      !write (*,*) YY(1,2), YY(20,21), YY(99,100)
      END PROGRAM CALC

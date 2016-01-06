program LU_Factorization_MatMult_TRIPLE_DO
!*****************************************************************************80
!Main program for the LU factorization and matrix multiplication tests  program.
!TRIPLE DO LOOPS are used for both matrix multiplication and LU decomposition subroutines
!Scope: to test the performance given by each change in order of the DO loops. 
!
! The end goal is to give a hint regarding the best practice order 
! for a user tailored algorithm that cannot be reduced at the standard libraries 
!  Licensing:This code is distributed under the GNU GPL license. 
!Built & Revisited: RV- 2015 
! This version is testing only the sequential algorithms. 
! Further development implies blocks version and OpenMP implementations


  implicit none
  integer i_LU,i_mm
  integer  n
  real*8 wtime_LU(1:7),wtime_mm(0:7),err
  character(5) order,seqen(0:7)
  common /legend/seqen
  seqen(0)='MMULT'; ! Intrinsic matrix multiplication option
  seqen(1)='I-J-K'; seqen(2 )='I-K-J';
  seqen(3)='J-K-I'; seqen(4 )='J-I-K';
  seqen(5)='K-I-J'; seqen(6)='K-J-I';
  seqen(7)='LAPCK'; ! if LAPACK package libraries are available; use MKL library

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Test TRIPLE DO LOOPS'
  write ( *, '(a)' ) ' Time Performance for LU Decomposition'
  write(*,'(a, 3x,6(2x,1x,1a,1x))' ) 'Nsize',seqen(1:6)
  do n=500,2000,500
    i_mm=0
    do  i_LU=1,6
    call test_case(n,i_mm,i_LU,wtime_mm(i_mm),wtime_LU(i_LU),err)
    enddo 
  write(*,'(i5,2x,6(2x,f7.3)    )' )  n,wtime_lu(1:6)
  enddo
  write ( *, '(a)' )' Time Performance for LU Matrix Multiplication'
  write(*,'(a, 3x,7(2x,1x,1a,1x))' ) 'Nsize',seqen(0:6)
  do n=500,2000,500
    i_LU=3
    do  i_mm=0,6
      call test_case(n,i_mm,i_LU,wtime_mm(i_mm),wtime_LU(i_LU),err)
    enddo 
  write(*,'(i5,2x,7(2x,f7.3)    )' )  n,wtime_mm(0:6)
  enddo

end

!=========================================
subroutine test_case ( n,i_mm,i_LU,wtime_mm,wtime_LU,err)
!*****************************************************************************80
! Test:  LU factorization option is given by i_LU; 
! Matrices multiplication option is given by i_mm
! Output: wall time wtime_xx for each subroutine   
  use omp_lib
  real :: a(1:n,1:n)
  real :: a0(1:n,1:n),al(1:n,1:n),au(1:n,1:n)
  integer i,j,i1,i_LU,i_mm
  integer  n
  real*8 wtime_LU,wtime_mm,err
  character(5) order,seqen(0:7)
  common /legend/seqen

!  Generate the matrix A - diagonal dominant
  call generate_matrix(n,a)
  a0=a ! Store the original matrix
  ORDER=seqen(i_LU)
  wtime_LU = omp_get_wtime ( )
  CALL ordered_LU_Factorization(order,n,a) 
  wtime_LU = omp_get_wtime ( ) - wtime_LU
  
!  Check the LU decomposition LU= original A
  wtime1 = omp_get_wtime ( )
 ! Zero initialization
  al(1:n,1:n)=0.;   au(1:n,1:n)=0.
! Split the L and U matrices  
  al(1,1)=1. ;  au(1,1:n)= a(1,1:n) 
  do i=2,n
     al(i,i)=1.; al(i,1:i-1)= a(i,1:i-1)
     au(i,i:n)= a(i,i:n)
  enddo 
  a(1:n,1:n)=0.;
    ORDER=seqen(i_mm)        
   wtime_mm = omp_get_wtime ( )
   call  ordered_MATMULT(order,n,n,n,al,au,a)   
   wtime_mm = omp_get_wtime ( ) - wtime_mm 
   err=0; 
   do i=1,n
   err = sum (  ( a(i,1:n)-a0(i,1:n))**2)+err
   enddo
   err =sqrt(err)/n
   
  
  return
end
!================================================

subroutine generate_matrix (  n, a)

!*****************************************************************************80
! USER CHOICE
! generates fractional values for the non-diagonal terms of the matrix 
  implicit none
  integer  n
  real :: a(n,n)
  integer  i, j
 !  Set the matrix A.
!
  do j = 1, n
    do i = 1, n
      a(i,j) = float(2*j-1)/float(j+3*i)
    end do
    a(j,j)=10
  end do
 
  return
end

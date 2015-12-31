!       Matrix Multiplication using TRIPLE DO LOOPS
!       It is testing the performance given by each change in order of the DO loops.
! The intent is not to optimize the matrix multiplication algorithm,
! when -opt-matmul is used all loops are indentified as a matmul intrinsic subroutine.
! The goal is to give a hint regarding the best practice order 
! for a user tailored algorithm that cannot be reduced at the standard libraries 
! Matrices: A(m,l), B(l,n), C(m,n) 
        SUBROUTINE ordered_MATMULT(order,m,l,n,A,B,C)
        CHARACTER(5) order
        INTEGER  m,l, n 
        REAL*4 A(m,*), B(l,*), C(m,*) !assumed shape memory allocation
        integer I, J, K
        
        IF(order.eq.'I-J-K')then
        DO I=1, m; DO J=1, n; DO  K=1, l
            C(I,J) = C(I,J)+  A(I,K) * B(K,J)            
        ENDDO; ENDDO;   ENDDO 
        ENDIF
        
        IF(order.eq.'I-K-J')then
        DO I=1, m;   DO  K=1, l ; DO J=1, n
            C(I,J) = C(I,J) + A(I,K) * B(K,J)
        ENDDO;ENDDO;ENDDO
        ENDIF    
        
        IF(order.eq.'K-I-J')then
        DO  K=1, l ;  DO I=1, m; DO J=1, n;  
            C(I,J) = C(I,J) + A(I,K) * B(K,J)
        ENDDO;ENDDO; ENDDO
        ENDIF  
        
        IF(order.eq.'K-J-I')then
        DO  K=1, l ; DO J=1, n; DO I=1, m;   
            C(I,J) = C(I,J) + A(I,K) * B(K,J)
        ENDDO;ENDDO;ENDDO
        ENDIF
        
        IF(order.eq.'J-I-K')then
        DO J=1, n; DO I=1, m;   DO  K=1, l ;
            C(I,J) = C(I,J) + A(I,K) * B(K,J)              
        ENDDO;ENDDO; ENDDO
        
        ENDIF
        IF(order.eq.'J-K-I')then
        DO J=1, N;  DO  K=1, l ; DO I=1, m;   
            C(I,J) = C(I,J) + A(I,K) * B(K,J)
        ENDDO;ENDDO; ENDDO
        ENDIF

!       Reference case       
        IF(order.eq.'MMULT')then
          C(1:m,1:n) = C(1:m,1:n) + MATMUL(A(1:m,1:l), B(1:l,1:n))
        ENDIF      
!       If MKL -library is available use BLAS3 -subroutine        
!        IF(order.eq.'LAPCK')then
!        call SGEMM('n','n',m,n,l,1.,a,m,b,l,0.,c,m)
!        ENDIF

      RETURN
      END

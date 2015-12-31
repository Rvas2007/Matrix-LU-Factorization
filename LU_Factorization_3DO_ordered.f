!*     LU -Factorization with  no pivoting for a square matrix; 
!*     TRIPLE DO LOOPS are used. Tests the performance given by each change in order of the DO loops. 
!*     Output: L - lower triangular matrix with unit diagonal and U -upper triangular matrix stored output A
!*     LAPACK's subroutine SGETRF - is used as reference and benchmark if available for performance analysis       
! The end goal is to give a hint regarding the best practice order 
! for a user tailored algorithm that cannot be reduced at the standard libraries 
! Note: Best performance should be given by the option "J-K-I" -optimized memory access 

        SUBROUTINE ordered_LU_Factorization(order,m,a)
        INTEGER  m
        REAL*4 A(m,*)!assumed shape memory allocation  
!        INTEGER ipv(m) ! pivot vector; -since no rows interchanges are performed => ipv(i)=i
        character(5) order
        INTEGER   I, J, K,INFO

        IF(order.eq.'I-J-K')then
        DO I=2, m; 
            DO J=1, I-1; 
                DO  K=1,J-1;
                a(I,J) = a(I,J)-  a(I,K) * a(K,J)   
                ENDDO; 
            a(I,J) = a(I,J)/a(J,J)     
            ENDDO;  
            DO J=I,m; ! J-K independent (interchangeable) loops, limits depend only of I 
                DO  K=1,I-1;
                a(I,J) = a(I,J)-  a(I,K) * a(K,J)   
                ENDDO;    
            ENDDO;
        ENDDO 
        ENDIF
        
        IF(order.eq.'I-K-J')then
        DO I=2, m; 
            DO K=1, I-1; 
            a(I,K) = a(I,K)/a(K,K)
                DO  J=K+1,m;
                a(I,J) = a(I,J)-  a(I,K) * a(K,J)            
                ENDDO; 
            ENDDO;  
        ENDDO 
        ENDIF    
        
        IF(order.eq.'J-I-K')then
        DO J=1, m; 
            DO I=2, J; 
                DO  K=1,I-1;
                a(I,J) = a(I,J)-  a(I,K) * a(K,J)            
                ENDDO;
            ENDDO;    
            DO  I=J+1,m;  ! I-K -independent loops, limits dependent only of J
                DO  K=1,J-1;
                a(I,J) = a(I,J)-  a(I,K) * a(K,J)            
                ENDDO;
            a(I,J) = a(I,J)/a(J,J)      
            ENDDO; 
        ENDDO;   
        ENDIF
        
        IF(order.eq.'J-K-I')then
        DO J=1, m;
            DO  K=1,J-1;
                DO  I=K+1,m;
                a(I,J) = a(I,J)-  a(I,K) * a(K,J)            
                ENDDO;
            ENDDO;  
            DO  I=J+1,m;
            a(I,J) = a(I,J)/a(J,J)      
            ENDDO; 
        ENDDO;   
        ENDIF

        IF(order.eq.'K-I-J')then
        DO  K=1, m-1 ;  
            DO I=K+1, m; 
            a(I,K) = a(I,K)/a(K,K)
                DO J=K+1, m;   ! I-J independent loops, limits dependent of K 
                a(I,J) = a(I,J)-  a(I,K) * a(K,J)
                ENDDO;
            ENDDO;
        ENDDO
        ENDIF    
        
        IF(order.eq.'K-J-I')then
         DO  K=1, m-1 ;  
            DO I=K+1, m; 
            a(I,K) = a(I,K)/a(K,K)
            ENDDO
            DO J=K+1, m;   ! J-I -independent loops - good candidate for medium grains parallelization  
                DO I=K+1, m; 
                a(I,J) = a(I,J)-  a(I,K) * a(K,J)
                ENDDO;
            ENDDO;
        ENDDO
        ENDIF
! If available use for reference LAPACK LU decomposition subroutine -unblocked code -        
!        IF(order.eq.'LAPCK')then
!        CALL SGETF2( m, m, a, m, ipv, INFO )
!        ENDIF
      RETURN
      END
!++++++++++++++++++++++++++++++++++++

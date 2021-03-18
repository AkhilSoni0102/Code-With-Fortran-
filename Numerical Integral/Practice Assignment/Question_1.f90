program  Akhil_Soni_180122004
    implicit none
    real:: a, b, Ans_T, Ans_S 
    Integer:: N, i 
    write(*,*) "Enter the Value of N" 
    read(*,*) N 
    read(*,*) a, b 
    
    call Trapezoidal_Rule(a,b,Ans_T,N)
    write(*,*) Ans_T 
    call Simpsons_Rule(a,b,Ans_S,N)
    write(*,*) Ans_S 
    contains 
    subroutine Trapezoidal_Rule(a,b,Ans_T,N)
        implicit none 
        real:: a, b, Ans_T, h
        integer:: i, N
        if(N < 1) then 
            Ans_T = 0.0
            return 
        end if 
        h = (b-a)/N 
        Ans_T = (f(a) + f(b))/2 
        do i = 1, N-1
            Ans_T = Ans_T + f(a+i*h)
        end do 
        Ans_T = h * Ans_T 
    end subroutine Trapezoidal_Rule

    SUBROUTINE Simpsons_Rule(a, b, Ans_S, N)
        IMPLICIT NONE 
        REAL:: a, b, Ans_S, h
        INTEGER:: N, i
        IF(N < 2) THEN
            Ans_S = 0.0
            RETURN  
        ENDIF
        Ans_S = 0.0
        h = (b-a)/N
        Ans_S = f(a) + f(b)
        DO i = 1, N/2
            Ans_S = Ans_S + 4 * f(a + (2*i - 1) * h) 
        ENDDO
        DO i = 1, (N+1)/2 - 1
            Ans_S = Ans_S + 2 * f(a + (2 * i * h))
        ENDDO
        Ans_S = (h/3) * Ans_S
    END SUBROUTINE Simpsons_Rule

    REAL FUNCTION f(x)
        IMPLICIT NONE 
        REAL:: x 
        f = 1./(1 + x**2)
        return
    END FUNCTION f
end program  Akhil_Soni_180122004
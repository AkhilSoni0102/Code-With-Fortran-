PROGRAM Akhil_Soni_180122004
    IMPLICIT NONE 
    REAL:: a, b, Ans_S, Coeff
    Real, parameter::  pi=3.141592
    INTEGER::  N   
    READ(*, *) N 
    READ(*, *) a, b

    do i = 1,20
        CALL Simpsons_Rule(a, b, Ans_S, N, i)
        write(*,*) Ans_S
    end do 
    SUBROUTINE Simpsons_Rule(a, b, Ans_S, N, x)
        IMPLICIT NONE 
        REAL:: a, b, Ans_S, h, x
        INTEGER:: N, i
        IF(N < 2) THEN
            Ans_S = 0.0
            RETURN  
        ENDIF
        h = (b-a)/N
        Ans_S = f(a) + f(b)
        DO i = 1, N/2
            Ans_S = Ans_S + 4 * f(a + (2*i - 1) * h) 
        ENDDO
        DO i = 1, (N+1)/2 - 1
            Ans_S = Ans_S + 2 * f(a + 2 * i * h)
        ENDDO
        Ans_S = (h/3) * Ans_S
    END SUBROUTINE Simpsons_Rule

    REAL FUNCTION f(x)
        IMPLICIT NONE 
        REAL:: x 
        f = (1./pi) * (cos(Theta - (x*sin(Theta)))
        return
    END FUNCTION f

END PROGRAM Akhil_Soni_180122004
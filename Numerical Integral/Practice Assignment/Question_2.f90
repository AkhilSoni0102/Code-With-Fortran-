PROGRAM Akhil_Soni_180122004
    IMPLICIT NONE 
    REAL:: a, b, Ans_T, Ans_S, Coeff
    INTEGER::  N   
    READ(*, *) N 
    READ(*, *) a, b

    CALL Trapezoidal_Rule(a, b, Ans_T, N)
    Coeff = sqrt(1./Ans_T)
    WRITE(*, *) Coeff

    CALL Simpsons_Rule(a, b, Ans_S, N)
    Coeff = sqrt(1./Ans_S)
    WRITE(*, *)  Coeff

    CONTAINS
    SUBROUTINE Trapezoidal_Rule(a, b, Ans_T, N)
        IMPLICIT NONE 
        REAL:: a, b, Ans_T, Xi, h
        INTEGER:: N, i
        IF(N < 1) THEN
            Ans_T = 0.0
            RETURN  
        ENDIF
        h = (b-a)/N
        Ans_T = (f(a) + f(b)) / 2.0
        DO i = 1, N-1
            Xi = a + i * h
            Ans_T = Ans_T + f(Xi)
        ENDDO
        Ans_T = h * Ans_T
    END SUBROUTINE Trapezoidal_Rule

    SUBROUTINE Simpsons_Rule(a, b, Ans_S, N)
        IMPLICIT NONE 
        REAL:: a, b, Ans_S, h
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
        f = EXP(-x**2)**2
        return
    END FUNCTION f

END PROGRAM Akhil_Soni_180122004
MODULE precision_Use
    ! dp = double precision
    INTEGER, PARAMETER:: dp = SELECTED_REAL_KIND(12)
END MODULE precision_Use

PROGRAM False_Position
    USE precision_Use
    IMPLICIT NONE 

    REAL(KIND = dp):: x0, x1, x2, f0, f1, f2, a, b, error, root
    WRITE(*, *) "Enter value of interval (a and b)"
    READ *, a, b

    x1 = a 
    x2 = b
    f1 = Func(x1)
    f2 = Func(x2)

    IF (f1*f2 .GT. 0.0) THEN 
        WRITE(*, *) "x1 and x2 do not bracket any root"
        RETURN 
    ENDIF
    error = 1e-6
    DO
        x0 = x1 - (x2 - x1)*f1 / (f2 - f1)
        f0 = Func(x0)
        IF (f0 .EQ. 0.0) THEN
            root = x0
            WRITE(*, *) "Root is: ", root
            RETURN
        ENDIF

        IF(f1 * f0 .LT. 0.0_dp) THEN 
            x2 = x0 
        ELSE 
            x1 = x0 
        ENDIF

        IF (abs((x2 - x1)/x2) .LT. error) THEN 
            root = (x1 + x2) / 2.0_dp
            WRITE(*, *) "Root is: ", root
            RETURN 
        ELSE 
            CONTINUE 
        ENDIF
    ENDDO 

    CONTAINS
    REAL(KIND = dp) FUNCTION Func(x)
        USE precision
        IMPLICIT NONE 
        REAL(KIND = dp):: x
        Func = x + 10 - x * Cosh(50/x)
        RETURN 
    END FUNCTION Func

END PROGRAM Bisection_Method
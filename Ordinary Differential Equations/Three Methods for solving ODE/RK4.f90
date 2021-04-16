MODULE precision_Use
    ! dp = double precision
    INTEGER, PARAMETER:: dp = SELECTED_REAL_KIND(12)
END MODULE precision_Use

PROGRAM RK4
    USE precision_Use
    IMPLICIT NONE 
    REAL(KIND = dp):: x0, y0, y1, h, b, a, s1, s2, s3, s4
    WRITE(*, *) "Enter value of interval a and b"
    READ *, a, b 
    write(*,*)"Initial value of the function: x0 and y0"
    read*, x0, y0 
    write(*,*)"Enter the value of h: "
    read*, h
    do 
        s1 = h*f_dash(x0,y0)
        s2 = h*f_dash(x0 + h/(2.), y0 + (s1/2.))
        s3 = h*f_dash(x0 + h/(2.), y0 + (s2/2.))
        s4 = h*f_dash(x0 + h, y0 + s3)
        y1 = y0 + ((s1 + 2*s2 + 2*s3 + s4)/6.)
        x0 = x0 + h 
        y0 = y1 
        if(x0 == b) then 
            write(*,*) "value: ", y1 
            return 
        endif
    end do 
    CONTAINS
    REAL(KIND = dp) FUNCTION f_dash(x, y)
        USE precision_Use
        IMPLICIT NONE 
        REAL(KIND = dp):: x, y 
        f_dash = x*x + y*y0
        RETURN 
    END FUNCTION f_dash
End program RK4
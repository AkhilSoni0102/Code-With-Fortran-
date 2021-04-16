MODULE precision_Use
    ! dp = double precision
    INTEGER, PARAMETER:: dp = SELECTED_REAL_KIND(12)
END MODULE precision_Use

PROGRAM RK2
    USE precision_Use
    IMPLICIT NONE 
    REAL(KIND = dp):: x0, y0, y1, h, b, a, sl, sr, x1 
    WRITE(*, *) "Enter value of interval a and b"
    READ *, a, b 
    write(*,*)"Initial value of the function: x0 and y0"
    read*, x0, y0 
    write(*,*)"Enter the value of h: "
    read*, h
    do 
        x1 = x0 + h 
        sl = f_dash(x0,y0)
        sr = f_dash(x1,y0 + h*sl)
        y1 = y0 + h * ((sl + sr)/2.)
        x0 = x1
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
        f_dash = x*x + y*y
        RETURN 
    END FUNCTION f_dash
End program RK2
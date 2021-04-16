MODULE precision_Use
    ! dp = double precision
    INTEGER, PARAMETER:: dp = SELECTED_REAL_KIND(12)
END MODULE precision_Use

PROGRAM Eulers_Method
    USE precision_Use
    IMPLICIT NONE 
    REAL(KIND = dp):: y_new,y_old,v_old,v_new,t, h, b 
    b = 1.
    t = 0.
    h = 0.25
    y_old = 1.
    v_old = 2.
    do 
        t = t + h 
        y_new = y_old + h * v_old
        v_new = v_old + h * ((-1)*v_old + sin(t*y_old))
        y_old = y_new 
        v_old = v_new 
        write(*,*) "value: ", v_new
        if(t == b) then
            return 
        end if 
    end do 
    CONTAINS
    REAL(KIND = dp) FUNCTION f_dash(x, y)
        USE precision_Use
        IMPLICIT NONE 
        REAL(KIND = dp):: x, y 
        f_dash = x*x + y*y
        RETURN 
    END FUNCTION f_dash
End program Eulers_Method
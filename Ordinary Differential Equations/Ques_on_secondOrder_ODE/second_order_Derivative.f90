MODULE precision_Use
    ! dp = double precision
    INTEGER, PARAMETER:: dp = SELECTED_REAL_KIND(12)
END MODULE precision_Use

MODULE Constants
    use precision_use 
    real, parameter:: pi = acos(-1.0)
end module Constants

PROGRAM Eulers_Method
    USE precision_Use
    use constants
    IMPLICIT NONE 
    integer:: i
    REAL(KIND = dp):: x_new,x_old,p_old,p_new, t, h
    t = 0.0
    h = 0.02*2*pi 
    x_old = 1.0
    p_old = 0.0
    OPEN(unit = 10, file = "output_Ques2.txt")
    write(10,*) 0,1,0,1
    do i = 1, 600
        t = t + h 
        x_new = x_old + h * ((-1)*p_old)
        p_new = p_old + h * ((-1)*x_old)
        x_old = x_new 
        p_old = p_new 
        WRITE(10, *) t/2.0*pi, x_old, p_old, (x_old*x_old + p_old*p_old)
    end do 
    CLOSE(10)
End program Eulers_Method
Program Akhil 
    IMPLICIT NONE
    real:: delta, X_not, E_trans, Nu, K_not, x, G_x 
    real, parameter:: pi = 3.141592
    CHARACTER*10:: filename
    INTEGER:: i
    delta = 0.25
    X_not = 10.0 
    Nu = 4056.5 
    E_trans = 0.008 
    x = 0.5 
    K_not = SQRT((2*Nu*E_trans) - (1./(2*delta*delta)))
    WRITE(filename, fmt = "(a, i0)") "output.", 1
    OPEN(UNIT = 10, FILE = filename, FORM = "formatted")
    do i = 1,128  
        G_x = (((1/pi*delta*delta)*(1./2.))*EXP((-1*((X-X_not)**2))/(delta*delta)))
        X = X + 0.15 
        WRITE(10, *)X, G_x
    end do 
    close(10)
end
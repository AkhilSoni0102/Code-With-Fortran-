!Compute the inegral of 1/t using the trapezoidal rule for limit 1 to 1.2 for N = 2
!Using Trapezioidal Rule
Program Akhil_Soni_180122004
    implicit none 
    real:: a, b
    integer:: i, N 
    real:: ans 
    a = 1
    b = 1.2
    N = 2
    ans = (1. + (1./b))/2.
    do i = 1,N-1
        ans = ans + 1./(1. + (i*(0.1)))
    end do 
    ans = ans * 0.1 
    print*, ans 
end program Akhil_Soni_180122004
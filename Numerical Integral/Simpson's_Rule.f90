!using Simpsons rule
program Akhil_Soni_180122004
    implicit none 
    real:: a, b 
    integer:: N, i 
    real:: ans 
    a = 1.
    b = 1.2 
    ans =  (0.2) / 6. 
    ans = 1. + (1./b) + (4. * (2./(a+b)))
    ans = ans * 0.1 
    ans = ans / 3.
    print*, ans 
end program Akhil_Soni_180122004
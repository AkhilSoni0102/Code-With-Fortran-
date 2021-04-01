!Consider the case of a finite well. If you remember, there are three regions to be
!considered and piecewise solutions are obtained. And taking boundary conditions into
!account, we obtain two transcedental equations, one for ”even” solutions and the other
!for ”odd” solutions. The equation for even solution is

!Solve this equation to find a root, using Newton-Raphson method. Start with a
!value of 3.5 for . Use a value of 4.0 for ρ.

Program Newton_Raphson_method
    Implicit None 
    Real:: x0,x1,f0,f,deriv_f
    x0 = 3.5 
    f = func(x0)
    deriv_f = deriv_func(x0)
    do 
        f = func(x0)
        deriv_f = deriv_func(x0)
        x1 = x0 - (f/deriv_f)
        f0 = func(x1)
        if(f0 == 0.) then 
            print*, "Ans", x1 
            return 
        endif 
        if(abs((x1-x0)/x1) < 1e-6) then 
            print*, "Ans", x1 
            return 
        else 
            x0 = x1 
        endif 
    enddo 

    Contains     
    real function func(x)
    real::x
    func = (x*tan(x)) - (sqrt(16.- x*x))
    end
    real function deriv_func(x)
    real:: x 
    deriv_func = ((x)/(cos(x)*cos(x))) + (tan(x)) + ((x)/(sqrt(16.0 - x*x)))
    end 
end program Newton_Raphson_method

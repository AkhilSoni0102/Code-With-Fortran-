Program Secant_Method
    Implicit None 
    Real:: x0,x1, x2, f0, f1, f2
    x1 = 3.2 
    x2 = 3.4
    f1 = func(x0)
    f2 = func(x1)
    do
        x0=x2-((x2-x1)/(f2-f1))*f2
        f0=func(x0)
        if(abs((x0-x2)/x0)<1e-6) then 
        print *,"Ans: ", x0
        return
        else 
        x1 = x2
        f1 = f2
        x2 = x0
        f2 = f0
        endif
    enddo

    Contains     
    real function func(x)
    real::x
    func = (x*tan(x)) - (sqrt(16.- x*x))
    end
end program Secant_Method

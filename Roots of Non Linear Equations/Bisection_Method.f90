Module Precision_Use 
    Integer, Parameter:: Dp = selected_real_kind(12)
end Module Precision_Use 

Program Akhil_Bisection_Method
    Use Precision_Use 
    Implicit None 
    Real(Kind = Dp):: x_0, x_1, x_2, f0,f1,f2, a, b, error, Root
    Integer:: c 
    c = 1
    error = 1e-6
    a = 0.5 
    b = 2.0 

    x_1 = a 
    x_2 = b 
    f1 = Funct(x_1)
    f2 = Func(x_2)
    if (f1*f2 > 0.0) then
        write(*, *) "x_1 and x_2 do not bracket the root"
        return 
    do
        x_0 = (x_1 + x_2)/2.0_DP
        f0 = Func(x_0)
        if (f0 == 0.0) then
            root = x_0
            write(*, *) "Root is: ", root
            return 
        endif

        if(f1 * f0 < 0.0_Dp) then
            x_2 = x_0 
        else 
            x_1 = x_0 
        endif

        if(abs((x_2 - x_1)/x_2) .LT. error) then
            root = (x_1 + x_2) / 2.0_Dp
            write(*, *) "Root: ", root
            return  
        else  
            c = c + 1 
        endif
    endo

    Contains 
    Real(Kind = Dp) Function Func(x)
        USE Precision_Use
        Implicit None 
        Real(Kind = Dp):: x
        Func = x**3 - 2 * Sin(x)
        Return 
    End Function Func

End Program Akhil_Bisection_Method
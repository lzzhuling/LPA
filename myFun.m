function s=myFun(dFt,dF,directF,u,d)   
    function F=func(d)
        F=dFt*max((dF*d+directF),0)+2*u*d;
    end
s=fsolve(@func,d);
end
function sfun=solvefun(A,Q,Qvar)

    function [f,Jac]=sfunf(x)
        f   =   exp(x).*(A*[1;exp(x)])-Q-Qvar.*(exp(2*x));
        Jac =   diag(exp(x))*A(:,2:3)*diag(exp(x))+diag(exp(x))*diag(A*[1;exp(x)])-diag(Qvar.*exp(2*x))*2;
    end
sfun=@sfunf;
end
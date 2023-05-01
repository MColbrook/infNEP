function [Z] = Newtons_Simple(x0,f,df,tol,maxit)

Z=x0;

ind=1:length(x0);
for j=1:maxit
    if isempty(ind)
        break
    else
        Z(ind)=Z(ind)-f(Z(ind))./df(Z(ind));
        ind=find(abs(f(Z))>tol);
    end
end
end
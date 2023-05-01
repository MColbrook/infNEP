function S = leg_normalize(n, lambda)

C = ones(1,n);
for ii = 1:round(2*lambda-1)
    C = C.*((1:n)+ii-1);
end

C = sqrt(((2^(1-2*lambda)*pi/(gamma(lambda)^2))./((0:n-1)+lambda)).*C);

S = spdiags(C',0,n,n);


end
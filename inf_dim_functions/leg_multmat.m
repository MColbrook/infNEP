function M = leg_multmat(n, f, lambda)
% MULTMAT  multiplication matrices for ultraS
%
%  M = MULTMAT(N, F, LAMBDA) forms the nxn multiplication matrix
%  representing the multiplication of F in the C^{(LAMBDA)} basis.
% 
%  M = MULTMAT(N, F, LAMBDA) also works when F is a vector of Chebyshev
%  coefficients.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

a = legcoeffs(f);

% Multiplying by a scalar is easy.
if ( numel(a) == 1 )
    M = a*speye(n);
    return
end

% Prolong or truncate coefficients
if ( numel(a) < n )
    a = [a ; zeros(n - numel(a), 1)];   % Prolong
else
    a = a(1:n);                         % Truncate.
end


% Convert to C^{lam}
a = ultraS.convertmat(n, 0.5, lambda - 1) * a;

m = 2*n; 
M0 = speye(m);

d1 = [1 (2*lambda : 2*lambda + m - 2)]./ ...
    [1 (2*((lambda+1) : lambda + m - 1))];
d2 = (1:m)./(2*(lambda:lambda + m - 1));
B = [d2' zeros(m, 1) d1'];
Mx = spdiags(B,[-1 0 1], m, m);
M1 = 2*lambda*Mx;

% Construct the multiplication operator by a three-term recurrence: 
M = a(1)*M0;
M = M + a(2)*M1;
for nn = 1:length(a) - 2
    M2 = 2*(nn + lambda)/(nn + 1)*Mx*M1 - (nn + 2*lambda - 1)/(nn + 1)*M0;
    M = M + a(nn + 2)*M2;
    M0 = M1;
    M1 = M2;
    if ( abs(a(nn + 3:end)) < eps ), break, end
end
M = M(1:n, 1:n); 


end
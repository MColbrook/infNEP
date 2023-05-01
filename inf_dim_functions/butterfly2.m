function [coeffs,fun,F] = butterfly2(n,c)
%BUTTERFLY  Quartic matrix polynomial with T-even structure.
%  [COEFFS,FUN,F] = nlevp('butterfly',N,C) generates the coefficient matrices
%  of an M^2-by-M^2 quartic matrix polynomial
%  P(lambda) = lambda^4*A_4 + lambda^2*A_3 + lambda^2*A_2 + lambda*A_1 + A_0
%  for which A_4 and A_2 are real and symmetric and
%  A_3 and A_1 are real and skew-symmetric (assuming C is real).
%  M is chosen such that M^2 approximates N.
%  The spectrum has a butterfly shape.
%  The default is N = 64.  C is an 10-by-1 vector of parameters.
%  The matrices are returned in a cell array: COEFFS = {A0, A1, A2, A3, A4}.
%  FUN is a function handle to evaluate the monomials 1,lambda,...,lambda^4
%  and their derivatives.
%  F is the function handle: lambda^4*A_4 + lambda^2*A_3 +
%  lambda^2*A_2 + lambda*A_1 + A_0.
%  This problem has the properties pep, real, parameter-dependent, T-even,
%  scalable, sparse, banded (8 bands).

%  Reference:
%  V. Mehrmann and D. Watkins. Polynomial eigenvalue problems with
%  Hamiltonian structure. Electron. Trans. Numer. Anal., 13:106-118, 2002.

if nargin < 1 || isempty(n)
    n = 64;
end  
if nargin < 2, c = [0.6 1.3 1.3 0.1 0.1 1.2 1.0 1.0 1.2 1.0]; end

% m = floor(sqrt(n));
% if abs(n-m^2)> abs(n-(m+1)^2),
%     m = m+1;
% end
m=n;

%N = gallery('jordbloc',m,0)';
N=spdiags(ones(m,1),-1,m,m);
M0 = (4*speye(m) + N + N')/6;
M1 = N - N';
M2 = -(2*speye(m)- N - N');
M3 = M1;
M4 = -M2;

A0 = c(1)*M0;
A1 = c(3)*M1;
A2 = c(5)*M2;
A3 = c(7)*M3;
A4 = c(9)*M4;

coeffs = {A0,A1,A2,A3,A4};

fun = @(lam) nlevp_monomials(lam,4);
F =  @(lam) coeffs{1} + lam*(coeffs{2} + lam*(coeffs{3} + lam*(coeffs{4} + 0*lam*coeffs{5}))); 
end
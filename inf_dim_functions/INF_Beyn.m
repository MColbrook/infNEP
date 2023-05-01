function [V,D,Res,cn] = INF_Beyn(N,m,quadpts,quadwts,jacobianpts,nlevp_solver,varargin)
%%% Solve-then-discretize version of Beyn's method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% N: size of discretization
% m: Initial upper bound for number of eigenvalues
% quadpts: quadrature points for integral
% quadwts: quadrature weights for integral
% jacobianpts: Jacobian at quadrature points
% nlevp_solver: solver of the form @(z,RHS) nlevp_solver(z, RHS), used to compute resolvents 

% OPTIONAL LABELLED INPUTS
% nlevp_residual: function of the form @(V,D) nlevp_residual(V,D) that computes residuals
% TOL: tolerance - throw away eigenpairs withe residual>TOL
% DIM: 0 if system is 1D (default), -k if system is k coupled 1D systems,
% 2 if system is 2D (Kronecker structure)

% OUTPUTS
% mu: smoothed measure at points X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Collect the optional inputs
p = inputParser;
addRequired(p,'N',@(x) (x==floor(x))&(x>0));
addRequired(p,'m',@(x) (x==floor(x))&(x>0));
addRequired(p,'quadpts',@isnumeric);
addRequired(p,'quadwts',@isnumeric);
addRequired(p,'jacobianpts',@isnumeric);
addRequired(p,'nlevp_solver',@(f) isa(f,'function_handle'))

addParameter(p,'nlevp_residual',@(V,D) 0*diag(D)-1,@(f) isa(f,'function_handle'));
addParameter(p,'TOL',1i,@isnumeric);
addParameter(p,'svdTOL',10^(-10),@isnumeric);
addParameter(p,'DIM',1,@(x) x==floor(x));
addParameter(p,'shift',sum(quadpts.*quadwts)/sum(quadwts),@isnumeric);

p.CaseSensitive = false;
parse(p,N,m,quadpts,quadwts,jacobianpts,nlevp_solver,varargin{:});
nlevp_residual=p.Results.nlevp_residual;
TOL=p.Results.TOL;
svdTOL=p.Results.svdTOL;
DIM=p.Results.DIM;
Z0=p.Results.shift;


if DIM==1
    V = randn(N, m);   % random functions
    V(2*m:end,:) = 0;  % cutoff a little to make it resolvable. 
    A1 = zeros(N,m);
    A0 = zeros(N,m);
elseif DIM==2
    V1 = randn(N^2, m);   % random functions
    V = zeros(N^2,m);
    m2 = ceil(sqrt(m));
    I1 = zeros(1,N);
    I1(1:m2) = 1;
    I1 = kron(I1,I1);
    V(I1==1,:)=V1(I1==1,:);
    A1 = zeros(N^2,m);
    A0 = zeros(N^2,m);
elseif DIM==0
    V = randn(2*m+1, m);
    V = [zeros(N-m,m);V;zeros(N-m,m)];
    A1 = zeros(2*N+1,m);
    A0 = zeros(2*N+1,m);
elseif DIM<0
    V=[];
    for jj=1:(-DIM)
        Va = randn(N, m);   % random functions
        Va(m+1:end,:) = 0;  % cutoff a little to make it resolvable
        V=[V;Va];
    end
    A1 = zeros((-DIM)*N,m);
    A0 = zeros((-DIM)*N,m);
end

% Compute the two main matrices A0 and A1 in Beyn's method:
% A0 = (1/(2pi*i))*int_{contour} T(z)^{-1}V dz
% A1 = (1/(2pi*i))*int_{contour} zT(z)^{-1}V dz

n = length(quadwts(:));
for j = 1:n
    U = nlevp_solver(quadpts(j), V);
    A0 = A0 + quadwts(j)*jacobianpts(j)*U; % Don't assume T is vectorized.
    A1 = A1 + quadwts(j)*jacobianpts(j)*(quadpts(j)-Z0)*U;
end
A0 = A0/2i/pi;
A1 = A1/2i/pi;

% Reduced SVD of A0:
[V0, S0, W0] = svd( A0, 0);
if m>1
    idx = find(diag(S0)/S0(1,1)>svdTOL,1,'last');
else
    idx = 1;
end
V0 = V0(:,1:idx);  S0 = S0(1:idx,1:idx); W0 = W0(:,1:idx);

% Compute eigenvalues of V0'*A1*W0*S0^{-1}, which correspond to those of T inside the contour.
[T,B] = balance(V0'*A1*W0*diag(1./diag(S0)));
[VV, DD] = eig(B);
V = T*VV;
D = DD + Z0*eye(size(DD,1),size(DD,2));

% Compute eigenfunctions
V=V0*V;

% Compute condition numbers
cn=zeros(size(VV,2),1);
if nargout > 3
    for j=1:size(VV,2)
        lam=DD(j,j);
        AA=V0'*A1*W0-lam*S0;
        warning('off','all')
        [w,~,v]=svds(AA,1,'smallest');
        warning('on','all')
        BB=-S0;
        cn(j)=(norm(V0'*A1*W0)+abs(lam)*norm(BB))*norm(v)*norm(w)/(abs(lam)*abs(w'*BB*v));
    end
end

% Compute residuals
Res = nlevp_residual(V,D);
if sum(Res)<0
    Res=[];
elseif isreal(TOL)
    II = find(Res<TOL);
    D=D(II,II);
    V=V(:,II);
    cn=cn(II);
end


    

    
end


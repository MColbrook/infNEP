clear
close all
rng(1)

%% Damped beam example
%% NOTE: Extended precision is used ONLY to avoid stability issues caused by linearisation for the discretised problems.
%%       The mp package used can be found here: https://www.advanpix.com/
%%       This recovers the errors shown in the figure in the paper.

addpath(genpath('nlevp')); % NLEVP collection
addpath(genpath('inf_dim_functions')); % inf dim functions from M. Colbrook and A. Townsend, "Avoiding discretization issues for nonlinear eigenvalue problems"
addpath(genpath('computed_data')); % computed data

%% Parameters

M = 1;
kappa = 1;
sigma = kappa/M;
NN = 500; % discretization size for finite-dimensional NEPs
ne = 100; % number of evals for plot

%% FEM discretize-then-solve approach

[cfs,~,F,dF,~] = loaded_string2(NN); % finite element discretisation, dF is function handle @(lambda) F'(lambda)
A = cfs{1}; B = -cfs{2}; C = cfs{3}; % matrices for rational pencil

D0 = A; D1 = -C; D2 = -B; % CORK linearization
A2 = [D0, D1-D2; sparse(NN,NN), speye(NN,NN)];
B2 = [sparse(NN,NN), -D2; speye(NN,NN),speye(NN,NN)];
% lam_FEM = double(eig(mp(full(A2)),mp(full(B2)),'vector')); % uncomment for extended prec
lam_FEM = eig(full(A2),full(B2),'vector');

lam_FEM = lam_FEM(lam_FEM>2);
[~,I]=sort(real(lam_FEM),'ascend');
lam_FEM=lam_FEM(I(1:NN/2));

% compute condition numbers
cn_FEM=zeros(ne,1); % condition number from rational NEP

pf = parfor_progress(ne);
pfcleanup = onCleanup(@() delete(pf));
for j=1:ne
    lam=lam_FEM(j);
    AA=F(lam);
    [w,~,v]=svds(AA,1,'smallest');
    BB=dF(lam);
    cn_FEM(j)=norm(v)*norm(w)/(abs(lam)*abs(w'*BB*v));
    parfor_progress(pf);
end

cn_FEM= cn_FEM.*(norm(full(A))+norm(full(B))*abs(lam_FEM(1:ne))+norm(full(C))*abs(lam_FEM(1:ne)./(lam_FEM(1:ne)-sigma)));

%% Rectangular chebyshev collocation discretize-then-solve approach

A = diffmat([NN-2 NN], 2, [0,1], 'chebkind2', 'dirichlet', 'neumann');
B = diffmat([NN-2 NN], 0, [0,1], 'chebkind2'); B = [zeros(1,NN);B;zeros(1,NN)];
C = diffmat([NN-2 NN], 2, [0,1], 'chebkind2', 'dirichlet', 'dirichlet'); C(1:end-1,:) = 0;

D0 = A; D1 = -C; D2 = -B;
A2 = [D0, D1-D2; zeros(NN,NN), eye(NN,NN)];
B2 = [zeros(NN,NN), -D2; eye(NN,NN),eye(NN,NN)];
% lam_COL = double(eig(mp(full(A2)),mp(full(B2)),'vector')); % uncomment for extended prec
lam_COL = eig(full(A2),full(B2),'vector');

lam_COL = lam_COL(lam_COL>2);
[~,I]=sort(real(lam_COL),'ascend');
lam_COL=lam_COL(I(1:NN/2));

% compute condition numbers
cn_COL=zeros(ne,1); % condition number from rational NEP

pf = parfor_progress(ne);
pfcleanup = onCleanup(@() delete(pf));
for j=1:ne
    lam=lam_COL(j);
    AA=F(lam);
    [w,~,v]=svds(AA,1,'smallest');
    BB=dF(lam);
    cn_COL(j)=norm(v)*norm(w)/(abs(lam)*abs(w'*BB*v));
    parfor_progress(pf);
end

cn_COL= cn_COL.*(norm(full(A))+norm(full(B))*abs(lam_COL(1:ne))+norm(full(C))*abs(lam_COL(1:ne)./(lam_COL(1:ne)-sigma)));

%% Legendre Galerkin discretize-then-solve approach

D = leg_diffmat2(NN);
D2 = 4*full(D*D);

A = [(-1).^(0:NN-1); (0:NN-1).*(1:NN); D2(1:end-2,:)];
B = [zeros(2,NN); eye(NN-2,NN)];
C = [zeros(1,NN); ones(1,NN); zeros(NN-2,NN)];
D0 = A; D1 = -C; D2 = -B;
A2 = [D0, D1-D2; zeros(NN,NN), eye(NN,NN)];
B2 = [zeros(NN,NN), -D2; eye(NN,NN),eye(NN,NN)];
% lam_GAL = double(eig(mp(full(A2)),mp(full(B2)),'vector')); % uncomment for extended prec
lam_GAL = eig(full(A2),full(B2),'vector');

lam_GAL = lam_GAL(lam_GAL>2);
[~,I]=sort(real(lam_GAL),'ascend');
lam_GAL=lam_GAL(I(1:NN/2));

% compute condition numbers
cn_GAL=zeros(ne,1); % condition number from rational NEP

pf = parfor_progress(ne);
pfcleanup = onCleanup(@() delete(pf));
for j=1:ne
    lam=lam_GAL(j);
    AA=F(lam);
    [w,~,v]=svds(AA,1,'smallest');
    BB=dF(lam);
    cn_GAL(j)=norm(v)*norm(w)/(abs(lam)*abs(w'*BB*v));
    parfor_progress(pf);
end

cn_GAL= cn_GAL.*(norm(full(A))+norm(full(B))*abs(lam_GAL(1:ne))+norm(full(C))*abs(lam_GAL(1:ne)./(lam_GAL(1:ne)-sigma)));

%% Ultraspherical discretize-then-solve approach

D2 = ultraS.diffmat(NN, 2);
S02 = ultraS.convertmat(NN, 0, 1);

BC=[1.^(0:NN-1); (-1).^(0:NN-1);(0:NN-1).^2];

A = [BC(2,:); 2*BC(3,:); 4*D2(1:end-2,:)];
B = [sparse(2,NN); S02(1:end-2,:)];
C = [sparse(1,NN); BC(1,:); sparse(NN-2,NN)];
D0 = A; D1 = -C; D2 = -B;
A2 = [D0, D1-D2; sparse(NN,NN), speye(NN,NN)];
B2 = [sparse(NN,NN), -D2; speye(NN,NN),speye(NN,NN)];
% lam_US = double(eig(mp(full(A2)),mp(full(B2)),'vector')); % uncomment for extended prec
lam_US = eig(full(A2),full(B2),'vector');

lam_US = lam_US(lam_US>2);
[~,I]=sort(real(lam_US),'ascend');
lam_US=lam_US(I(1:NN/2));

% compute condition numbers
cn_US=zeros(ne,1); % condition number from rational NEP

pf = parfor_progress(ne);
pfcleanup = onCleanup(@() delete(pf));
for j=1:ne
    lam=lam_US(j);   
    AA=F(lam);
    [w,~,v]=svds(AA,1,'smallest');
    BB=dF(lam);
    cn_US(j)=norm(v)*norm(w)/(abs(lam)*abs(w'*BB*v));
    parfor_progress(pf);
end

cn_US= cn_US.*(norm(full(A))+norm(full(B))*abs(lam_US(1:ne))+norm(full(C))*abs(lam_US(1:ne)./(lam_US(1:ne)-sigma)));

%% Compute analytic eigenvalues for error plots

nn = (0:NN)';
ASYMP = pi^2*(nn+1/2).^2;

f = @(z) cos(sqrt(z)) + sqrt(z)./(z-kappa).*sin(sqrt(z));
df = @(z)  -sin(sqrt(z))./(2*sqrt(z)) + (sqrt(z).*(z-kappa).*cos(sqrt(z))-(kappa+z).*sin(sqrt(z)))./(2*sqrt(z).*(z-kappa).^2);
% lam_inf = Newtons_Simple(mp(ASYMP),f,df,10^(-20),10000); % uncomment for extended prec
lam_inf = Newtons_Simple(ASYMP,f,df,10^(-20),10000);   % Newton's method

%% Compute eigenvalues using INF_Beyn

spec = [];
cn_BEYN=[];
TOL = 10^(-4);

pf = parfor_progress(ne);
pfcleanup = onCleanup(@() delete(pf));
for j=1:ne
    LHS=double(max(lam_inf(j)-5,4));
    RHS=double(lam_inf(j+1)+5); % can easily alter for a larger or smaller cluster of eigenvalues - multiple eigenvalues (recommend > 5) needed for condition number to make sense
    cntr=(LHS+RHS)/2;
    L=(RHS-LHS)/2;
    n=max(2*ceil(L),100);

    contour = @(t) cntr + L*cos(t)+1i*sin(t);
    jacobian = @(t) -L*sin(t)+1i*cos(t);
    tpts = linspace(0,2*pi,n+1)+2*pi/(2*n); tpts(end)=[];  % trap rule
    quadpts = contour(tpts);
    quadwts = (2*pi)/n * ones(1,n);
    m = 5; % increase for slightly increased stability
    [V,D,~,cc] = INF_Beyn(NN,m,quadpts,quadwts,jacobian(tpts),@(z,RHS) nlevp_solver(z, RHS, M, kappa),...
        'nlevp_residual',@(V,D) nlevp_residual(V, D, M,kappa),'TOL',TOL);

    [D,I] = sort(real(diag(D)));
    spec=[spec;D];
    cn_BEYN=[cn_BEYN;cc(:)];
    parfor_progress(pf);
end
%% clean up computed spectra
spec2=1;
cn2=1;
for j=1:length(spec)
    if min(abs(spec(j)-spec2))>0.1
        I=find(abs(spec-spec(j))<0.1);
        ccc=cn_BEYN(I);
        spec2=[spec2(:);mean(spec(I(ccc==min(ccc))))];
        cn2=[cn2(:);mean(cn_BEYN(I))];
    end
end
spec2(1)=[];
cn2(1)=[];

%% Plot results

figure
loglog(cn_FEM,'.','markersize',20)
hold on
loglog(cn_COL,'.','markersize',20)
loglog(cn_GAL,'.','markersize',20)
loglog(cn_US,'.','markersize',20)
loglog(cn2,'k.','markersize',20)
xlim([1,ne])
ax = gca; ax.FontSize = 14;
%%

figure
loglog(abs(lam_FEM-lam_inf(1:length(lam_FEM)))./abs(lam_inf(1:length(lam_FEM))),'.','markersize',20)
hold on
semilogy(abs(lam_COL-lam_inf(1:length(lam_COL)))./abs(lam_inf(1:length(lam_COL))),'.','markersize',20)
semilogy(abs(lam_GAL-lam_inf(1:length(lam_GAL)))./abs(lam_inf(1:length(lam_GAL))),'.','markersize',20)
semilogy(abs(lam_US-lam_inf(1:length(lam_US)))./abs(lam_inf(1:length(lam_US))),'.','markersize',20)
semilogy(abs(spec2-lam_inf(1:length(spec2)))./abs(lam_inf(1:length(spec2))),'k.','markersize',20)
xlim([1,ne])
ax = gca; ax.FontSize = 14;




%% Solvers for infBEYN

function Sol = nlevp_solver(z, RHS, MM, kappa)
N = size(RHS,1);
D2 = ultraS.diffmat(N, 2);
S02 = ultraS.convertmat(N, 0, 1);

T = 4*D2 + z*S02;

B=[1.^(0:N-1); (-1).^(0:N-1)];
for ii=1
    b=(0:N-1).^2;
    for j=1:ii-1
        b=b.*(((0:N-1).^2-j^2)/(2*j+1));
    end
    B=[B; b; (-1).^(ii+(0:N-1)).*b];
end

BC = [B(2,:);
    2*B(3,:)+z*kappa*MM/(z-kappa)*B(1,:)];

T = [speye(2,N); T(1:end-2,:)];
            
RHSb = S02*RHS;
RHSb = [zeros(2,size(RHSb,2)); RHSb(1:end-2,:)];

U = eye(N,2); V = BC-eye(2,N);

Q1 = T\RHSb;
Q2 = T\U;
Sol = Q1-Q2*((eye(2,2)+V*Q2)\(V*Q1));
end

function Res = nlevp_residual(V, D, MM, kappa)
Z = diag(D);
Res = zeros(length(Z),1);

N = size(V,1);

D2 = ultraS.diffmat(N, 2);
S02 = ultraS.convertmat(N, 0, 1);

for j=1:length(Z)
    z = Z(j);

    T = 4*D2 + z*S02;

    B=[1.^(0:N-1); (-1).^(0:N-1)];
    for ii=1
        b=(0:N-1).^2;
        for jj=1:ii-1
            b=b.*(((0:N-1).^2-jj^2)/(2*jj+1));
        end
        B=[B; b; (-1).^(ii+(0:N-1)).*b];
    end

    BC = [B(2,:);
        2*B(3,:)+z*kappa*MM/(z-kappa)*B(1,:)];

    T = [BC; T(1:end-2,:)];

    u = T*V(:,j);
    Res(j) = norm(u)/norm(V(:,j));
end
end



function [percent, elapsed] = parfor_progress(N, varargin)
narginchk(1, 2);

if isnumeric(N) && N > 0
    if nargin > 1
        file = varargin{1};
    else
        file = [tempname '_parfor.txt'];
    end
    fid = fopen(file, 'w');
    if fid < 0, error('Could not open file for writing (perms?): %s', file); end
    % write N, start time (0.1s resolution), and iteration (0 for now)
    progress = [N floor(now*864000) 0];
    fprintf(fid, '%d\n%d\n0\n', progress);
    fclose(fid);
    computeprogress(progress, false, true);
    percent = file;
    elapsed = 0;
else
    file = N;
    if ~ischar(file) || ~exist(file, 'file')
        error('Not initialized. See HELP PARFOR_PROGRESS.');
    end
    % update (read and write) in one go
    fid = fopen(file, 'r+');
    progress = fscanf(fid, '%f');                   % read the 3 values
    progress(3) = progress(3) + 1;                  % update iteration number
    fseek(fid, 0, 'bof');
    fprintf(fid, '%d\n%d\n%d\n', round(progress));  % write back to file
    fclose(fid);
    [percent, elapsed] = computeprogress(progress, true, nargout == 0);
end

end

function [percent, elapsed] = computeprogress (progress, update, show)
    elapsed = (now - progress(2)/864000) * 86400;   % compute elapsed seconds
    percent = progress(3) / progress(1);
    if percent == 0
        duration = 0;
        remaining = 0;
    else
        duration = elapsed / percent;
        remaining = duration - elapsed;
    end
    if show
        r = humantime(remaining);
        e = humantime(elapsed);
        t = humantime(duration);
        s = sprintf('%8.2f%%, %s (el), %s (rem), %s (tot)\n', ...
                percent * 100, e, r, t);
        if update, back = repmat(char(8),1,length(s)); else back = ''; end
        fprintf('%s%s', back, s);
    end
end

function t = humantime (s)
    if s < 60, t = sprintf('%4.1fs', s);
    elseif s < 3600, t = sprintf('%4.1fm', s/60);
    elseif s < 86400, t = sprintf('%4.1fh', s/3600);
    elseif s < 604800, t = sprintf('%4.1fd', s/86400);      % 86400 = 1 day = 24 * 3600
    elseif s < 2629800, t = sprintf('%4.1fw', s/604800);    % 604800 = 1 week = 7 * 86400
    elseif s < 31557600, t = sprintf('%4.1fM', s/2629800);  % 2629800 = 1 average month = 365.25/12 * 86400
    else t = sprintf('%4.1fy', s/31557600);                 % 31557600 = 1 average year = 365.25 * 86400
    end
end
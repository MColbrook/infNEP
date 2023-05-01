clear
close all

%% Damped beam example

addpath(genpath('nlevp')); % NLEVP collection
addpath(genpath('inf_dim_functions')); % inf dim functions from M. Colbrook and A. Townsend, "Avoiding discretization issues for nonlinear eigenvalue problems"
addpath(genpath('computed_data')); % computed eigenvalues

%% Eigenvalues of discretised operator

n = 100; 
cfs = nlevp('damped_beam',n);

M = full(cfs{3}); D = full(cfs{2}); K = full(cfs{1}); 
m = norm(M); d = norm(D); k = norm(K);
gamma = sqrt(k/m); delta = 2/(k+d*gamma); 
tM = gamma^2*delta*M; 
tD = gamma*delta*D; 
tK = delta*K;

Z = 0*M; 
A = [ -tM Z ;  Z tK ];
B = [ Z tM ; tM tD ]; 

mu = eig(full(A),-full(B));
lambda = gamma*mu;

%% Physical parameters of the problem

E = 7e10; 
I = 0.05*0.005^3/12; 
L = 1; 
rhoAL = 0.674; 
c = 5;

beta=-c/(E*I);
alpha=-rhoAL/(E*I);

% %% Compute eigenvalues using
% spec=[];
% m = 5;
% N = 100;
% TOL = 10^(-4);
% 
% NN = 100;
% mu = @(x) alpha^(1/4)*L/2*sqrt(x);
% nn=(1:NN);
% cc=abs(8*1i*sqrt(abs(alpha))/(beta*L));
% ASYMP=(nn*pi-pi/2)*1i -1./((nn*pi)*cc);
% ASYMP=(ASYMP/mu(1)).^2;     % use asymptotics for contours
% 
% for jj=1:20
%     % construct contour
%     cntr=ASYMP(jj); LL=10; if jj>10; LL=1; end
%     n=50;
%     contour = @(t) cntr + LL*cos(t)+LL*1i*sin(t);
%     jacobian = @(t) -LL*sin(t)+LL*1i*cos(t);
%     tpts = linspace(0,2*pi,n+1)+2*pi/(2*n); tpts(end)=[];  % trap rule
%     quadpts = contour(tpts);
%     quadwts = (2*pi)/n * ones(1,n);
% 
%     [V,D,Res] = INF_Beyn(N,m,quadpts,quadwts,jacobian(tpts),@(z,RHS) nlevp_solver(z, RHS, alpha,beta),...
%         'nlevp_residual',@(V,D) nlevp_residual(V, D, alpha, beta),'DIM',-1,'TOL',TOL);
%     spec = [spec;diag(D);conj(diag(D))];
% end



%%
load('damped_beam_spec.mat')

figure(100);    lineH = plot(1,1,1,1,1,1,1,1,1,1);  C = get(lineH, 'Color');   close(100) % grab colours from matlab

close all
figure
plot(real(lambda),(imag(lambda)),'.','markersize',25,'color',C{2}), hold on,
plot(real(spec),(imag(spec)),'ko','markersize',3,'linewidth',2), hold on,
ax = gca; ax.FontSize = 14;
axis([-15,0,-10^6,10^6])

figure
plot(real(lambda),(imag(lambda)),'.','markersize',25,'color',C{2}), hold on,
plot(real(spec),(imag(spec)),'ko','markersize',3,'linewidth',2), hold on,
ax = gca; ax.FontSize = 14;
axis([-7.5,-7.4,-10^5/2,10^5/2])

%% Eigenfunction plot

cfs = nlevp('damped_beam',100);
[X,lambda] = polyeig(cfs{:});

I=find(imag(lambda)<0); X(:,I)=[];  lambda(I)=[];
I=find(real(lambda)<-0.1);  lambda1=lambda(I);  X1=X(:,I); % separate into the two sets of eigenvalues in UHP
I=find(real(lambda)>-0.1);  lambda2=lambda(I);  X2=X(:,I);

[~,I]=sort(imag(lambda1),'ascend'); lambda1=lambda1(I); X1=X1(:,I);
[~,I]=sort(imag(lambda2),'ascend'); lambda2=lambda2(I); X2=X2(:,I);

spec3=spec(201:2:end);

m = 5;  N = 100;    TOL = 10^(-4);

% lambda_10^1 eigenfunction
rng(6)

cntr=spec3(10); LL=1; n=50;
contour = @(t) cntr + LL*cos(t)+LL*1i*sin(t);   jacobian = @(t) -LL*sin(t)+LL*1i*cos(t);
tpts = linspace(0,2*pi,n+1)+2*pi/(2*n); tpts(end)=[];  % trap rule
quadpts = contour(tpts);    quadwts = (2*pi)/n * ones(1,n);

[V1,~,~] = INF_Beyn(N,m,quadpts,quadwts,jacobian(tpts),@(z,RHS) nlevp_solver(z, RHS, alpha,beta),...
    'nlevp_residual',@(V,D) nlevp_residual(V, D, alpha, beta),'DIM',-1,'TOL',TOL);

ua=chebfun(double(V1(1:N/2)),[0,1/2],'coeffs');
ub=chebfun(double(V1(N/2+1:end)),[1/2,1],'coeffs');
u=chebfun({ua,ub},[0 0.5 1]); % concatenate the domains
u=u/sqrt(sum(u*conj(u))); % normalise the eigenfunction

x1=chebpts(100, [0,1]);
XS=X2(:,10); % discretised eigenfunction
XS=[0,transpose(XS(1:end-1)),0,XS(end)];
PC = fnval(pwch((0:size(X,1)/2)/(size(X,1)/2),XS(1:2:end),XS(2:2:end)),x1); % reconstruct the discretised eigenfunction from cubic Hermite polynomials
PC = chebfun(PC,[0,1]);
PC = PC/sqrt(sum(PC*conj(PC))); % normalise the eigenfunction

[x2,w]=chebpts(200, [0,1]); % quadrature to compute subspace angle
U1=sqrt(w(:)).*u(x2(:));
U2=sqrt(w(:)).*PC(x2(:));

subspace(U1,U2)
abs(spec3(10)-lambda2(10))

PC = PC/PC(0.455)*u(0.455); % normalise the eigenfunction to match inf dim one

figure
xxx=0:0.005:1;
subplot(2,1,1)
plot(real(u),'linewidth',2)
hold on
plot(xxx,real(PC(xxx)),'.','markersize',12)
legend({'InfBeyn','Discretized'},'interpreter','latex','fontsize',14,'location','northeast','Orientation','horizontal');
ax = gca; ax.FontSize = 14;

subplot(2,1,2)
plot(imag(u),'linewidth',2)
hold on
plot(xxx,imag(PC(xxx)),'.','markersize',12)
ax = gca; ax.FontSize = 14;


% lambda_10^2 eigenfunction
rng(6) % this only effects the normalisation of the eigenfunctions for plots

cntr=spec2(10); LL=1;   n=50;
contour = @(t) cntr + LL*cos(t)+LL*1i*sin(t);   jacobian = @(t) -LL*sin(t)+LL*1i*cos(t);
tpts = linspace(0,2*pi,n+1)+2*pi/(2*n); tpts(end)=[];  % trap rule
quadpts = contour(tpts);    quadwts = (2*pi)/n * ones(1,n);

[V1,~,~] = INF_Beyn(N,m,quadpts,quadwts,jacobian(tpts),@(z,RHS) nlevp_solver(z, RHS, alpha,beta),...
    'nlevp_residual',@(V,D) nlevp_residual(V, D, alpha, beta),'DIM',-1,'TOL',TOL);

ua=chebfun(double(V1(1:N/2)),[0,1/2],'coeffs');
ub=chebfun(double(V1(N/2+1:end)),[1/2,1],'coeffs');
u=chebfun({ua,ub},[0 0.5 1]); % concatenate the domains
u=u/sqrt(sum(u*conj(u))); % normalise the eigenfunction

x1=chebpts(100, [0,1]);
XS=X1(:,10); % discretised eigenfunction
XS=[0,transpose(XS(1:end-1)),0,XS(end)];
PC = fnval(pwch((0:size(X,1)/2)/(size(X,1)/2),XS(1:2:end),XS(2:2:end)),x1); % reconstruct the discretised eigenfunction from cubic Hermite polynomials
PC = chebfun(PC,[0,1]);
PC = PC/sqrt(sum(PC*conj(PC))); % normalise the eigenfunction

[x2,w]=chebpts(200, [0,1]); % quadrature to compute subspace angle
U1=sqrt(w(:)).*u(x2(:));
U2=sqrt(w(:)).*PC(x2(:));

subspace(U1,U2)
abs(spec2(10)-lambda1(10))

PC = PC/PC(0.455)*u(0.455); % normalise the eigenfunction to match inf dim one

figure
xxx=0:0.005:1; % x points for plot
subplot(2,1,1)
plot(real(u),'linewidth',2)
hold on
plot(xxx,real(PC(xxx)),'.','markersize',12)
ax = gca; ax.FontSize = 14;
subplot(2,1,2)
plot(imag(u),'linewidth',2)
hold on
plot(xxx,imag(PC(xxx)),'.','markersize',12)
ax = gca; ax.FontSize = 14;


%% Solvers for infBEYN

function Sol = nlevp_solver(z, RHS, alpha, beta)
N = size(RHS,1)/2;
D4 = ultraS.diffmat(N, 4); 
S04 = ultraS.convertmat(N, 0, 3); 
T = (- D4 + alpha/(4^4)*z^2*S04);

B=[1.^(0:N-1); (-1).^(0:N-1)];
for ii=1:3
    b=(0:N-1).^2;
    for j=1:ii-1
        b=b.*(((0:N-1).^2-j^2)/(2*j+1));
    end
    B=[B; b; (-1).^(ii+(0:N-1)).*b];
end

BC = [B(2,:),zeros(1,N); B(6,:),zeros(1,N); zeros(1,N),B(1,:); zeros(1,N),B(5,:);
        B(1,:),-B(2,:); B(3,:),-B(4,:); B(5,:),-B(6,:); B(7,:)+beta/(4^3)*z*B(1,:),-B(8,:)];

T = [speye(4,N),sparse(4,N); sparse(4,N),speye(4,N); T(1:end-4,:), sparse(N-4,N); sparse(N-4,N),T(1:end-4,:) ];

RHSb = [S04*RHS(1:N,:);S04*RHS(N+1:end,:)];
RHSb = [zeros(8,size(RHSb,2)); RHSb([1:N-4,N+1:end-4],:)];

U = eye(N*2,8); V = BC-[eye(4,N),zeros(4,N);zeros(4,N),eye(4,N)];
Q1 = (T\RHSb);
Q2 = (T\U);
Sol = (Q1-Q2*((eye(8,8)+V*Q2)\(V*Q1)));
end


function Res = nlevp_residual(V, D, alpha, beta)
Z = diag(D);
Res = zeros(length(Z),1);

N = size(V,1)/2;

D4 = ultraS.diffmat(N, 4); 
S04 = ultraS.convertmat(N, 0, 3); 

B=[1.^(0:N-1); (-1).^(0:N-1)];
for ii=1:3
    b=(0:N-1).^2;
    for j=1:ii-1
        b=b.*(((0:N-1).^2-j^2)/(2*j+1));
    end
    B=[B; b; (-1).^(ii+(0:N-1)).*b];
end

for j=1:length(Z)
    z = Z(j);
    T = (- D4 + alpha/(4^4)*z^2*S04);

    BC = [B(2,:),sparse(1,N); B(6,:),sparse(1,N); sparse(1,N),B(1,:); sparse(1,N),B(5,:);
            B(1,:),-B(2,:); B(3,:),-B(4,:); B(5,:),-B(6,:); B(7,:)+beta/(4^3)*z*B(1,:),-B(8,:)];

    T = [BC; T(1:end-4,:), sparse(N-4,N); sparse(N-4,N),T(1:end-4,:) ];
   
    u = (T)*V(:,j);
    Res(j) = (norm(u)/norm(V(:,j)))/abs(z);
end
end

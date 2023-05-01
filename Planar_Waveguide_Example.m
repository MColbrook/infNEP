clear
close all

%% Planar waveguide example

addpath(genpath('nlevp')); % NLEVP collection
addpath(genpath('inf_dim_functions')); % inf dim functions from M. Colbrook and A. Townsend, "Avoiding discretization issues for nonlinear eigenvalue problems"
addpath(genpath('computed_data')); % computed eigenvalues and pseudospectra

cfs = nlevp('planar_waveguide');
lambda = polyeig(cfs{:}); % eigenvalues of discretised operator

load('waveguide_evals.mat') % load computed eigenvalues (these were found using infBEYN below), rrr are residuals

% physical parameters of problem
n0 = 1.5; n1 = 1.66; n2 = 1.6; n3 = 1.53; n4 = 1.66; n5 = 1.0;
ksq = (2*pi/0.6328)^2; deltasq = ksq*(n0^2-n5^2);
betasq2 = ksq*(n0^2+n5^2)/2 + deltasq.^2./(16*spec.^2) + spec.^2;
mu = betasq2/ksq;

alpha0 = 1i*(deltasq./(4*lambda)-lambda); alpha5 = -1i*(deltasq./(4*lambda)+lambda);
betasq = ksq*(n0^2+n5^2)/2 - (alpha0.^2+alpha5.^2)/2;
mu_fin = betasq/ksq; % mu for discretised problem

%% Example code to compute eigenvalues. In practice, use a larger contour to locate eigenvalues then hone in with smaller contours.
m = 20; 
N = 500;
TOL = 10^(-3);

cntr=0.8i; L=0.2; n=100;
contour = @(t) cntr + L*cos(t)+L*1i*sin(t);
jacobian = @(t) -L*sin(t)+L*1i*cos(t);
tpts = linspace(0,2*pi,n+1)+2*pi/(2*n); tpts(end)=[];  % trap rule
quadpts = contour(tpts);
quadwts = (2*pi)/n * ones(1,n);

[~,D,~] = INF_Beyn(N,m,quadpts,quadwts,jacobian(tpts),@(z,RHS) nlevp_solver(z, RHS),'DIM',-4,...
    'nlevp_residual',@(V,D) nlevp_residual(V, D),'TOL',TOL);
spec = diag(D);

%% Code to compute essential spectrum of underlying problem on [-inf,inf]

avec = [0:0.001:60,60.001:0.0001:70,70.001:0.001:100,10.^(2.001:0.001:30)];
avec = avec - ksq*1.5^2;
dp = ksq*(n0^2+n5^2)/2;
dm = ksq*(n0^2-n5^2)/2;
avec = avec+dp;

zv1 = sqrt(-avec/2+sqrt(avec.^2-dm/2)/2);
zv2 = sqrt(-avec/2-sqrt(avec.^2-dm/2)/2);
zv1 = zv1(abs(imag(zv1))<300);
zv2 = zv2(abs(imag(zv2))<300);

I=find(abs(zv2(2:end)-zv2(1:end-1))>2);
zv2=[zv2(1:I),imag(zv2(I))*1i,NaN,zv2(I+1:end)];

%% Plot results

figure(100);    lineH = plot(1,1,1,1,1,1,1,1,1,1,1,1,1,1);  C = get(lineH, 'Color');   close(100) % grab colours from matlab

% lambda plane

figure
plot(real(lambda),imag(lambda),'.','markersize',25,'color',C{2})
hold on
plot(real(spec),imag(spec),'ko','markersize',3,'linewidth',2), hold on,
plot(real(spec),-imag(spec),'ko','markersize',3,'linewidth',2), hold on,
plot(real(zv1),imag(zv1),'-','linewidth',2,'color',C{6})
plot(-real(zv1),-imag(zv1),'-','linewidth',2,'color',C{6})
plot(real(zv2),imag(zv2),'-','linewidth',2,'color',C{6})
plot(-real(zv2),-imag(zv2),'-','linewidth',2,'color',C{6})
ax = gca; ax.FontSize = 14;
axis([-6,6,-5,5])

load('waveguide_evals.mat','spec') % load full evals


figure
plot(real(lambda),imag(lambda),'.','markersize',25,'color',C{2})
hold on
plot(real(spec),imag(spec),'ko','markersize',3,'linewidth',2), hold on,
plot(real(spec),-imag(spec),'ko','markersize',3,'linewidth',2), hold on,
plot(real(zv1),imag(zv1),'-','linewidth',2,'color',C{6})
plot(-real(zv1),-imag(zv1),'-','linewidth',2,'color',C{6})
plot(real(zv2),imag(zv2),'-','linewidth',2,'color',C{6})
plot(-real(zv2),-imag(zv2),'-','linewidth',2,'color',C{6})
ax = gca; ax.FontSize = 14;
axis([-.002,.004,0,0.6])

figure
plot(real(lambda),imag(lambda),'.','markersize',25,'color',C{2})
hold on
plot(real(spec),imag(spec),'ko','markersize',3,'linewidth',2), hold on,
plot(real(spec),-imag(spec),'ko','markersize',3,'linewidth',2), hold on,
plot(real(zv1),imag(zv1),'-','linewidth',2,'color',C{6})
plot(-real(zv1),-imag(zv1),'-','linewidth',2,'color',C{6})
plot(real(zv2),imag(zv2),'-','linewidth',2,'color',C{6})
plot(-real(zv2),-imag(zv2),'-','linewidth',2,'color',C{6})
ax = gca; ax.FontSize = 14;
ylim([-225,225])
xlim([-10,10])




% mu plane

mu1 = dp/ksq + dm./(8*ksq*zv1.^2) + zv1.^2/ksq;

figure
plot(real(mu_fin),imag(mu_fin),'.','markersize',25,'color',C{2})
hold on
plot(real(mu),imag(mu),'ko','markersize',3,'linewidth',2)
plot(real(mu),-imag(mu),'ko','markersize',3,'linewidth',2)
plot(real(mu1),imag(mu1),'-','linewidth',2,'color',C{6})
ax = gca; ax.FontSize = 14;
axis([-500,50,-6,6])

figure
plot(real(mu_fin),imag(mu_fin),'.','markersize',25,'color',C{2})
hold on
plot(real(mu),imag(mu),'ko','markersize',3,'linewidth',2)
plot(real(mu),-imag(mu),'ko','markersize',3,'linewidth',2)
plot(real(mu1),imag(mu1),'-','linewidth',2,'color',C{6})
ax = gca; ax.FontSize = 14;
axis([-80,40,-1*0,0.6])

a2 = axes();
a2.Position = [0.71 0.75 0.18 0.15]; % xlocation, ylocation, xsize, ysize
plot(real(mu_fin),imag(mu_fin),'.','markersize',25,'color',C{2})
hold on
plot(real(mu),imag(mu),'ko','markersize',3,'linewidth',2)
plot(real(mu),-imag(mu),'ko','markersize',3,'linewidth',2)
plot(real(mu1),imag(mu1),'-','linewidth',2,'color',C{6})
axis([2.1,2.7,-1*0,0.025])



%% Convergence plots

figure
load('spec_conv3.mat')
Err=abs(spec2-spec);

betasq2 = ksq*(n0^2+n5^2)/2 + deltasq.^2./(16*spec.^2) + spec.^2;
spec2 = sqrt(betasq2/ksq);

semilogy(nvec(5:5:end),Err(5:5:end)/abs(spec2),'-o','markersize',10,'linewidth',1)
hold on
load('spec_conv2.mat')

Err=abs(spec2-spec);
betasq2 = ksq*(n0^2+n5^2)/2 + deltasq.^2./(16*spec.^2) + spec.^2;
spec2 = sqrt(betasq2/ksq);
semilogy(nvec(5:5:end),Err(5:5:end)/abs(spec2),'-o','markersize',10,'linewidth',1)
load('spec_conv1.mat')
Err=abs(spec2-spec);
betasq2 = ksq*(n0^2+n5^2)/2 + deltasq.^2./(16*spec.^2) + spec.^2;
spec2 = sqrt(betasq2/ksq);
semilogy(nvec(5:5:end),Err(5:5:end)/abs(spec2),'-o','markersize',10,'linewidth',1)
axis([0,140,10^(-17),0.01])
ax = gca; ax.FontSize = 14;

legend({'$\lambda=3.0163\times 10^{-3} + 0.4967i$','$\lambda=1.1987\times 10^{-4} + 0.1001i$','$\lambda=8.2287\times 10^{-7} + 0.0100i$'},'location','northeast','fontsize',14,'interpreter','latex')




%% Solvers for infBEYN


function Sol = nlevp_solver(z, RHS)
N = size(RHS,1)/4;
D2 = ultraS.diffmat(N, 2); 
S02 = ultraS.convertmat(N, 0, 1);

% parameters of problem
n0 = 1.5; n1 = 1.66; n2 = 1.6; n3 = 1.53; n4 = 1.66; n5 = 1.0;
ksq = (2*pi/0.6328)^2; deltasq = ksq*(n0^2-n5^2);

alpha0 = 1i*(deltasq/(4*z)-z); alpha5 = -1i*(deltasq/(4*z)+z);
betasq = ksq*(n0^2+n5^2)/2 - (alpha0^2+alpha5^2)/2;

T1 = 4^2*D2 + (ksq*n1^2-betasq)*S02;
T2 = 4^2*D2 + (ksq*n2^2-betasq)*S02;
T3 = 4^2*D2 + (ksq*n3^2-betasq)*S02;
T4 = 4^2*D2 + (ksq*n4^2-betasq)*S02;

B=[1.^(0:N-1); (-1).^(0:N-1)];
for ii=1:1
    b=(0:N-1).^2;
    for j=1:ii-1
        b=b.*(((0:N-1).^2-j^2)/(2*j+1));
    end
    B=[B; b; (-1).^(ii+(0:N-1)).*b];
end

ZZ = zeros(1,N);
BC = [4*B(4,:)-1i*alpha0*B(2,:),ZZ, ZZ, ZZ;
        ZZ,ZZ,ZZ, 4*B(3,:)+1i*alpha5*B(1,:);
        B(1,:),-B(2,:), ZZ,ZZ;
        B(3,:),-B(4,:), ZZ,ZZ;
        ZZ,B(1,:),-B(2,:), ZZ;
        ZZ,B(3,:),-B(4,:), ZZ;
        ZZ,ZZ,B(1,:),-B(2,:);
        ZZ,ZZ,B(3,:),-B(4,:)];

T = [speye(2,N),sparse(2,3*N);
                sparse(2,N),speye(2,N),sparse(2,2*N);
                sparse(2,2*N),speye(2,N),sparse(2,N);
                sparse(2,3*N),speye(2,N);
                T1(1:end-2,:), sparse(N-2,3*N);
                sparse(N-2,N),T2(1:end-2,:), sparse(N-2,2*N);
                sparse(N-2,2*N),T3(1:end-2,:), sparse(N-2,N);
                sparse(N-2,3*N),T4(1:end-2,:)];
            
RHSb = [S02*RHS(1:N,:);S02*RHS(N+1:2*N,:);S02*RHS(2*N+1:3*N,:);S02*RHS(3*N+1:end,:)];
RHSb = [zeros(8,size(RHSb,2));
            RHSb([1:(N-2),(N+1):(2*N-2),(2*N+1):(3*N-2),(3*N+1):(4*N-2)],:)];

U = eye(N*4,8); V = BC-[eye(2,N),zeros(2,3*N);
                        zeros(2,N),eye(2,N),zeros(2,2*N);
                        zeros(2,2*N),eye(2,N),zeros(2,N);
                        zeros(2,3*N),eye(2,N)];

Q1 = T\RHSb;
Q2 = T\U;
Sol = Q1-Q2*((eye(8,8)+V*Q2)\(V*Q1));
end


function Res = nlevp_residual(V, D)

Z = diag(D);
Res = zeros(length(Z),1);

N = size(V,1)/4;

D2 = ultraS.diffmat(N, 2); 
S02 = ultraS.convertmat(N, 0, 1);

% parameters of problem
n0 = 1.5; n1 = 1.66; n2 = 1.6; n3 = 1.53; n4 = 1.66; n5 = 1.0;
ksq = (2*pi/0.6328)^2; deltasq = ksq*(n0^2-n5^2);

B=[1.^(0:N-1); (-1).^(0:N-1)];
for ii=1:1
    b=(0:N-1).^2;
    for j=1:ii-1
        b=b.*(((0:N-1).^2-j^2)/(2*j+1));
    end
    B=[B; b; (-1).^(ii+(0:N-1)).*b];
end

for j=1:length(Z)
    z = Z(j);
    alpha0 = 1i*(deltasq/(4*z)-z);
    alpha5 = -1i*(deltasq/(4*z)+z);
    betasq = ksq*(n0^2+n5^2)/2 - (alpha0^2+alpha5^2)/2;

    T1 = 4^2*D2 + (ksq*n1^2-betasq)*S02;
    T2 = 4^2*D2 + (ksq*n2^2-betasq)*S02;
    T3 = 4^2*D2 + (ksq*n3^2-betasq)*S02;
    T4 = 4^2*D2 + (ksq*n4^2-betasq)*S02;

    ZZ = zeros(1,N);
    BC = [4*B(4,:)-1i*alpha0*B(2,:),ZZ, ZZ, ZZ;
            ZZ,ZZ,ZZ, 4*B(3,:)+1i*alpha5*B(1,:);
            B(1,:),-B(2,:), ZZ,ZZ;
            B(3,:),-B(4,:), ZZ,ZZ;
            ZZ,B(1,:),-B(2,:), ZZ;
            ZZ,B(3,:),-B(4,:), ZZ;
            ZZ,ZZ,B(1,:),-B(2,:);
            ZZ,ZZ,B(3,:),-B(4,:)];

    T = [BC;
                    T1(1:end-2,:), sparse(N-2,3*N);
                sparse(N-2,N),T2(1:end-2,:), sparse(N-2,2*N);
                sparse(N-2,2*N),T3(1:end-2,:), sparse(N-2,N);
                sparse(N-2,3*N),T4(1:end-2,:)];
    
    u = T*V(:,j);
    Res(j) = norm(u)/norm(V(:,j));
end
end
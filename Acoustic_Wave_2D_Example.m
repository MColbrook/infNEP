clear
close all

%% Acoustic Wave 2d example

addpath(genpath('nlevp')); % NLEVP collection
addpath(genpath('inf_dim_functions')); % inf dim functions from M. Colbrook and A. Townsend, "Avoiding discretization issues for nonlinear eigenvalue problems"
impedance=1;

%% Eigenvalues of discretised operators in region where there is no spectrum

DS=zeros(3,1);      % actual discretisation sizes
lamFD=cell(3,1);    % eigenvalues for each discretisation size
nvec=[50,250,500];

for j=1:length(nvec)
    [cfs,fun,F,xcoeffs] = nlevp('acoustic_wave_2d',nvec(j),impedance);
    lamFD{j} = polyeig(cfs{1},cfs{2},cfs{3});
    DS(j)=size(cfs{1},1);
end

%% Plot the results

figure(100);    lineH = plot(1,1,1,1,1,1,1,1,1,1);  CC = get(lineH, 'Color');   close(100) % grab colours from matlab

figure
plot(real(lamFD{1}),imag(lamFD{1}),'.','markersize',30);
hold on
plot(real(lamFD{2}),imag(lamFD{2}),'.','markersize',30);
plot(real(lamFD{3}),imag(lamFD{3}),'.','markersize',30);
c = get(gca,'Color');
axis([-10,10,0.6,3.5])
ax = gca; ax.FontSize = 14;

%% Compute inf dim eigenvalues in region where there is spectrum

% contour parameters
n=200;                                                  % number of quad points
L=7;cntr =7;
contour = @(t) cntr + L*cos(t)+3i*sin(t);               % contour
jacobian = @(t) -L*sin(t)+3i*cos(t);
tpts = linspace(0,2*pi,n+1)+2*pi/(2*n); tpts(end)=[];   % trap rule
quadpts = contour(tpts);
quadwts = (2*pi)/n * ones(1,n);

mm = 80;            % max no. of eigenvalues inside contour.
N = 1000;           % max discretisation size
TOL = 10^(-5);      % tolerance for inf dim residual
spec=[];

for m=1:19
    [~,D,~] = INF_Beyn(N,mm,quadpts,quadwts,jacobian(tpts),@(z,RHS) nlevp_solver(z, RHS, impedance,m),...
        'nlevp_residual',@(V,D) nlevp_residual(V,D,impedance,m),'TOL',TOL);
    spec = [spec(:);diag(D)];
end

%% Plot results

figure
semilogy(real(lamFD{3}),imag(lamFD{3}),'.','markersize',30,'color',CC{3})
hold on
plot(real(spec),imag(spec),'ko','markersize',4,'linewidth',2)
ax = gca; ax.FontSize = 14;
axis([0,10,0,1])

%% Slow convergence for purely imaginary eigenvalues

impedance=0.8;
spec2=[];

tfact=sign(real(1i*sqrt(1/(impedance^2-1))))/sqrt(1-1/impedance^2); % from asymptotic form

for m=1:10
    cntr = m*tfact/2;    L1 = abs(tfact/2);
    contour1 = @(t) cntr + L1*(cos(t)+1i*sin(t));    jacobian1 = @(t) L1*(-sin(t)+1i*cos(t));
    tpts1 = linspace(0,2*pi,n+1)+2*pi/(2*n); tpts1(end)=[];
    quadpts1 = contour1(tpts1);    quadwts1 = (2*pi)/n * ones(1,n);

    [~,D,~] = INF_Beyn(N,mm,quadpts1,quadwts1,jacobian1(tpts1),@(z,RHS) nlevp_solver(z, RHS, impedance,m),...
        'nlevp_residual',@(V,D) nlevp_residual(V,D,impedance,m),'TOL',TOL);
    spec2 = [spec2(:);diag(D)];
end

[~,I] = sort(imag(spec2),'ascend');
spec2=spec2(I);

lamFD2=cell(3,1);
nvec=[30,42,156];
DS2=zeros(3,1);     % check actual discretisation sizes

for j=1:length(nvec)
    [cfs,fun,F,xcoeffs] = nlevp('acoustic_wave_2d',nvec(j),impedance);
    lamFD2{j} = polyeig(cfs{1},cfs{2},cfs{3});
    DS2(j)=size(cfs{1},1);
end



%%

for j=1:3
    figure
    plot(real(lamFD2{j}),imag(lamFD2{j}),'.','markersize',30,'color',CC{5})
    hold on
    plot(real(spec2),imag(spec2),'ko','markersize',4,'linewidth',2)
    ax = gca; ax.FontSize = 16;
    axis([-5,5,0.5,3])
end


%% Solvers for infBEYN

function Sol = nlevp_solver(z,RHS,impedance,m)
N = size(RHS,1);
D2 = ultraS.diffmat(N, 2); 
S02 = ultraS.convertmat(N, 0, 1); 
T = - 4*D2 - (4*pi^2*z^2-m^2*pi^2)*S02; 
nn = 0:N-1;

BC = [  (-1).^nn ; 2*(nn.^2) + 2*pi*1i*z/impedance ];

T = [speye(2,N); T(1:(end-2),:)];
RHS=S02*RHS;
RHS = [zeros(2,size(RHS,2)); RHS(1:end-2,:)];
U = eye(N,2); V = BC-eye(2,N);
Q1 = T\RHS;
Q2 = T\U;
Sol = Q1-Q2*((eye(2,2)+V*Q2)\(V*Q1));
end


function Res = nlevp_residual(V,D,impedance,m)
Z = diag(D);
Res = zeros(length(Z),1);

N = size(V,1);
D2 = ultraS.diffmat(N, 2); 
S02 = ultraS.convertmat(N, 0, 1);
nn = 0:N-1;

for j=1:length(Z)
    z = Z(j);
    T = - 4*D2 - (4*pi^2*z^2-m^2*pi^2)*S02;
    T = [[  (-1).^nn ; 2*(nn.^2) + 2*pi*1i*z/impedance ]; T(1:(end-2),:)];
    u = T*V(:,j);
    Res(j) = norm(u);
end

end


clear
close all

%% Orr-Sommerfeld example

addpath(genpath('nlevp')); % NLEVP collection
addpath(genpath('inf_dim_functions')); % inf dim functions from M. Colbrook and A. Townsend, "Avoiding discretization issues for nonlinear eigenvalue problems"
addpath(genpath('computed_data')); % computed eigenvalues and pseudospectra

%% Set parameters

R = 5772.22;
omega = 0.264002;
N = 64;

x_pts = (0:0.1:1.51)-0.001;    y_pts = -0.5:0.1:2; % increase resolution for plot
z_pts = kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts = z_pts(:);		% complex points where we compute pseudospectra

v = (10.^(-10:0.5:10)); % contours for pseudospectra plot

%% Discretise-then-solve approach

RESfin = zeros(size(z_pts));

[cfs,~,F,G,H,W] = nlevp('orr_sommerfeld2',N,R,omega);
lambda = polyeig(cfs{:});

pf = parfor_progress(length(z_pts));
pfcleanup = onCleanup(@() delete(pf));
for kk=1:length(z_pts)    
    [V,D] = eig(W*H(1));
    W1=(V*diag(sqrt(diag(D))))/V; % weight matrix for inner product
    W2=(V*diag(1./sqrt(diag(D))))/V; % weight matrix for inner product
    RESfin(kk)=min(svd(W1*(G(z_pts(kk))\F(z_pts(kk)))*W2));
    parfor_progress(pf);
end

figure
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(reshape(RESfin,length(y_pts),length(x_pts)))),log10(v));
cbh=colorbar;
set(gca,'YDir','normal')
colormap gray
ax=gca; ax.FontSize=14; axis tight;  
hold on
plot(real(lambda),imag(lambda),'.r','markersize',12)
axis([0,1.5,-0.5,2])
clim([-7,0])

%% Infinite-dimensional approach

RES=PseudoSpec(z_pts,N+37,R,omega);   % +37 is due to rectangular truncation - the number of columns will be N

figure
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(reshape(RES,length(y_pts),length(x_pts)))),log10(v));
cbh=colorbar;
set(gca,'YDir','normal')
colormap gray
ax=gca; ax.FontSize=14; axis tight;  
axis([0,1.5,-0.5,2])
clim([-7,0])

%% Load high res data
clear

v = (10.^(-10:0.5:10));

load('OS_fin1.mat')
figure
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(reshape(RESfin,length(y_pts),length(x_pts)))),log10(v));
cbh=colorbar;
set(gca,'YDir','normal')
colormap gray
cbh.Ticks=log10(10.^(-7:1:0));
cbh.TickLabels=["1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1"];
hold on
plot(real(lambda),imag(lambda),'.r','markersize',14)
ax=gca; ax.FontSize=20; axis tight;
axis([0,1.5,-0.5,2])
clim([-7,0])

%%

clear

v = (10.^(-10:0.5:10));

load('OS_inf1.mat')
figure
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(reshape(RES,length(y_pts),length(x_pts)))),log10(v));
cbh=colorbar;
set(gca,'YDir','normal')
colormap gray
cbh.Ticks=log10(10.^(-7:1:0));
cbh.TickLabels=["1e-7","1e-6","1e-5","1e-4","1e-3","1e-2","1e-1","1"];
ax=gca; ax.FontSize=20; axis tight;
axis([0,1.5,-0.5,2])
clim([-7,0])





function r = PseudoSpec(z_pts,N,R,omega)
nn=0:N-1;
r=zeros(size(z_pts));

S1 = leg_normalize(N, 0.5);
D = leg_diffmat2(N);
D2 = D*D;

A = spdiags(((2*nn+7)./(2*nn+3))',0,N,N)...
    -2*[sparse(2,N);spdiags(((2*nn(1:end-2)+5)./(2*nn(1:end-2)+3))',0,N-2,N-2), sparse(N-2,2)]...
    +[sparse(4,N);speye(N-4,N)];
A = S1*A;
[Q,~] = qr(A);

pf = parfor_progress(length(z_pts));
pfcleanup = onCleanup(@() delete(pf));
for jj=1:length(z_pts)
    z = z_pts(jj);
    
    L0 = (D2-z^2*speye(N,N))*(D2-z^2*speye(N,N))/R - 1i*(  leg_multmat(N,chebfun(@(x) z*(1-x^2)-omega),0.5)*(D2 - z^2*speye(N,N))  +2*z*speye(N,N)); 
    
    B0 = -(D2-z^2*speye(N,N));
    B0 = [ones(1,N);  (-1).^(1:N);  B0(1:end-2,:)];
    B1 = -(D2-abs(1)^2*speye(N,N));
    
    L = (L0/S1)*Q;
    [ii,~,~] = find(L(:,1));
    
    L = L(:,1:end-(max(ii)));
    L=[zeros(2,size(L,2));L(1:end-2,:)];
    L=S1*(B0\L);

    L=full(L);
    L(abs(L)<10^(-40))=0;
    [lower,~] = bandwidth(L);

    L = L(:,1:end-lower);
    
    L=L'*(S1*(B1/S1))*L;

    W2=full(Q'*S1*B1*(S1\Q));
    W2=full(W2(1:end-2,1:end-2));
    W2=(W2+W2')/2;
    [VG,DG]=eig(W2);
    W3=VG*diag(sqrt(1./diag(DG)))*VG';
    W3=W3(1:size(L,2),1:size(L,2));
    
    L=W3*L*W3;
    
    r(jj) = sqrt(min(svd(L)));
    parfor_progress(pf);
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




clear
close all

%% Butterfly example

addpath(genpath('nlevp')); % NLEVP collection
addpath(genpath('inf_dim_functions')); % inf dim functions from M. Colbrook and A. Townsend, "Avoiding discretization issues for nonlinear eigenvalue problems"
addpath(genpath('computed_data')); % computed eigenvalues and pseudospectra


% %% Compute eigenvalues of discretisation
% 
% n = 200;
% c = [0.2i 0 1.3 0 0.1 0 1.0 0 0 0];
% 
% cfs = butterfly2(n,c);
% lambda = polyeig(cfs{1:4}); % eigenvalues of discretised operator
% 
% 
% %% Compute pseudospectra
% 
% x_pts=-3:0.1:3;    y_pts=-3:0.1:3; % increase resolution for plot
% z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);		% complex points where we compute pseudospectra
% RES=reshape( PseudoSpec(z_pts,n,c),length(y_pts),length(x_pts));
% RESfin=reshape( PseudoSpecF(z_pts,n,c),length(y_pts),length(x_pts));


%% Plot results
load('buttefly_data.mat') % load data
z_pts=kron(x_pts,ones(length(y_pts),1))+1i*kron(ones(1,length(x_pts)),y_pts(:));    z_pts=z_pts(:);

close all
v=(10.^(-16:1:0));
vv=[-16,-1];

figure
RR=max(reshape(RES,length(y_pts),length(x_pts)),min(v));
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RR)),log10(v));
cbh=colorbar;
set(gca,'YDir','normal')
colormap gray
ax=gca; ax.FontSize=14; axis tight;
cbh.Ticks=log10(10.^(-16:2:-2));
cbh.TickLabels=["1e-16","1e-14","1e-12","1e-10","1e-08","1e-06","1e-04","1e-02","1e-02"];
axis([-3,3,-3,3])
clim(vv)

figure

RR=max(reshape(RESfin,length(y_pts),length(x_pts)),min(v));
contourf(reshape(real(z_pts),length(y_pts),length(x_pts)),reshape(imag(z_pts),length(y_pts),length(x_pts)),log10(real(RR)),log10(v));
cbh=colorbar;
set(gca,'YDir','normal')
colormap gray
ax=gca; ax.FontSize=14; axis tight;
cbh.Ticks=log10(10.^(-16:2:-2));
cbh.TickLabels=["1e-16","1e-14","1e-12","1e-10","1e-08","1e-06","1e-04","1e-02","1e-02"];
hold on
plot(real(lambda),imag(lambda),'.r','markersize',12)
axis([-3,3,-3,3])
clim(vv)



%% Pseudospectra routines



function r = PseudoSpec(z_pts,N,c)
[~,~,F] = butterfly2(N+2,c); % N+2 due to rectangular truncation
r=zeros(size(z_pts));

pf = parfor_progress(length(z_pts));
pfcleanup = onCleanup(@() delete(pf));
for jj=1:length(z_pts) % this can be trivially parallelised
    B=F(z_pts(jj));
    B2=B';
    B=B(:,2:end-1); % rectangular truncations
    B2=B2(:,2:end-1); % rectangular truncations
    r(jj) = min(svds(B,1,'smallest'),svds(B2,1,'smallest'));
    parfor_progress(pf); % keep track of computation time
end
end


function r = PseudoSpecF(z_pts,N,c)
[~,~,F] = butterfly2(N,c); 
r=zeros(size(z_pts));

pf = parfor_progress(length(z_pts));
pfcleanup = onCleanup(@() delete(pf));
for jj=1:length(z_pts) % this can be trivially parallelised
    B=F(z_pts(jj));
    r(jj) = svds(B,1,'smallest');
    parfor_progress(pf); % keep track of computation time
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

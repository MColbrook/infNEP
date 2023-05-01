clear
close all

%% Acoustic Wave 1d example

addpath(genpath('nlevp'));
impedance = 1;%1.0001;

%% Compute eigenvalues of discretisation

nvec = [10,100,500];
F = zeros(4,length(nvec));
ct = 1;
lamFD=cell(length(nvec),1);

for n = nvec
    cfs = nlevp('acoustic_wave_1d',n,impedance); 
    lamFD{ct} = polyeig(cfs{:});
    ct = ct+1;
end

%% Plot the results

figure
plot(real(lamFD{1}),imag(lamFD{1}),'.','markersize',30)
hold on
plot(real(lamFD{2}),imag(lamFD{2}),'.','markersize',30)
plot(real(lamFD{3}),imag(lamFD{3}),'.','markersize',30)

if impedance~=1
    exact = (atan(1i*impedance) + (-n*2:n*2)*pi)/2/pi;
    plot(real(exact),imag(exact),'ko','markersize',4,'linewidth',2)
end

ax = gca; ax.FontSize = 14;
axis([-10,10,0,1])


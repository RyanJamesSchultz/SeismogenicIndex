% Simple script to plot the fit of Seismogenic Index data.
clear;

% Load in and predefine some variables.
load('SI_data.mat');
b=1.05;
Mc=1.3;
Vs=0;
Ve=0;
RED=[247,96,96]/256;

% Compute Seismogenic Index.
[VcumEQs, sigma,S,S_err, R2_s, Vcum_fits,Ncum_fits, Vss]=SeismogenicIndex( Top,Vcum, T,M, b,Mc, Vs,Ve);
N=length(VcumEQs);

% Output fit statistics to user prompt.
fprintf('Sigma: %0.3f ± %0.3f (R² %0.3f)\n', S, S_err, R2_s);


% Plot cumulative event count vs. volume.
figure(1); clf;
subplot(121);
loglog(VcumEQs, 1:N, 'o','Color', 'k', 'MarkerFaceColor', RED, 'MarkerSize', 6 ); hold on;
loglog(Vcum_fits,Ncum_fits,'-','Color','k');
xlabel('Injection Volume (m^3)'); ylabel('Cumulative Seismic Event Count');
ylim([0.7 1.3*N]); xlim([0.7*min(VcumEQs) 1.3*max(VcumEQs)]);

% Plot Seismogenic Index vs. volume.
subplot(122);
plot(VcumEQs, sigma, '-o', 'Color',RED); hold on;
plot(Vcum_fits, repmat(S, size(Vcum_fits)), '-','Color', 'k');
plot(Vcum_fits, repmat(S+S_err, size(Vcum_fits)), '--','Color', 'k');
plot(Vcum_fits, repmat(S-S_err, size(Vcum_fits)), '--','Color', 'k');
xlabel('Injection Volume (m^3)'); ylabel('Seismogenic Index');
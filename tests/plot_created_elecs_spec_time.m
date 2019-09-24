clear all
close all
clc

global bins;
bins=logspace(1,log10(40000),32);

yy = importdata('created_electrons');

energies = yy(:,1);
times = yy(:,2);



%% energy
[n_simu_drm,xout] = make_spectrum(energies);


histogram('BinEdges',xout,'BinCounts',n_simu_drm,'DisplayStyle','stairs','LineWidth',2);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

xlabel('kinetic energy (keV)')
ylabel('de/dn spectrum (a.u.)')
title('created electron energy spectrum')

grid on


%% time
figure
[counts,time_bins] = histcounts(times,64);
histogram('BinEdges',time_bins,'BinCounts',counts,'DisplayStyle','stairs','LineWidth',2);
set(gca, 'YScale', 'log')
xlabel('Time (microsecond)')
ylabel('Counts per bin')
title('created electron time histogram')
grid on

%%
function [n_simu,xout] = make_spectrum(ener_list)
global bins;

[n_simu,xout] = histcounts(ener_list,bins);
n_simu = n_simu ./ diff(bins);

end
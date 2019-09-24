clear all
close all
clc

global bins;
bins=logspace(1,log10(40000),64);

yy = importdata('created_photons');

energies = yy(:,1);
times = yy(:,2);




[n_simu_drm,xout] = make_spectrum(energies);


histogram('BinEdges',xout,'BinCounts',n_simu_drm,'DisplayStyle','stairs','LineWidth',2,'DisplayName','Simulation');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

xlabel('energy (keV)')
ylabel('de/dn spectrum (a.u.)')

grid on


function [n_simu,xout] = make_spectrum(ener_list)
global bins;

[n_simu,xout] = histcounts(ener_list,bins);
n_simu = n_simu ./ diff(bins);

end
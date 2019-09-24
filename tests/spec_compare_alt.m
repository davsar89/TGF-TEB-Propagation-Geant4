clear all
close all
clc


load('/Home/siv29/dsa030/Desktop/GEANT4/ASIM-TEB/analysis/spectrum_TGF/spec_TGF_phot_10km.mat');
histogram('BinEdges',spec_TGF_phot.bins,'BinCounts',spec_TGF_phot.counts_per_ener/mean(spec_TGF_phot.counts_per_ener),'DisplayStyle','stairs','LineWidth',1.5);
hold on

load('/Home/siv29/dsa030/Desktop/GEANT4/ASIM-TEB/analysis/spectrum_TGF/spec_TGF_phot_12km.mat');
histogram('BinEdges',spec_TGF_phot.bins,'BinCounts',spec_TGF_phot.counts_per_ener/mean(spec_TGF_phot.counts_per_ener),'DisplayStyle','stairs','LineWidth',1.5);

load('/Home/siv29/dsa030/Desktop/GEANT4/ASIM-TEB/analysis/spectrum_TGF/spec_TGF_phot_15km.mat');
histogram('BinEdges',spec_TGF_phot.bins,'BinCounts',spec_TGF_phot.counts_per_ener/mean(spec_TGF_phot.counts_per_ener),'DisplayStyle','stairs','LineWidth',1.5);

load('/Home/siv29/dsa030/Desktop/GEANT4/ASIM-TEB/analysis/spectrum_TGF/spec_TGF_phot_17km.mat');
histogram('BinEdges',spec_TGF_phot.bins,'BinCounts',spec_TGF_phot.counts_per_ener/mean(spec_TGF_phot.counts_per_ener),'DisplayStyle','stairs','LineWidth',1.5);

load('/Home/siv29/dsa030/Desktop/GEANT4/ASIM-TEB/analysis/spectrum_TGF/spec_TGF_phot_19km.mat');
histogram('BinEdges',spec_TGF_phot.bins,'BinCounts',spec_TGF_phot.counts_per_ener/mean(spec_TGF_phot.counts_per_ener),'DisplayStyle','stairs','LineWidth',1.5);


xlabel('Energy (keV)')
ylabel('dn/de spectrum (keV^{-1})')
grid on
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

legend('source at 10 km','source at 12 km','source at 15 km','source at 17 km','source at 19 km')
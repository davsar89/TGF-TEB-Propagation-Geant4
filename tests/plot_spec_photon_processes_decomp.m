clearvars
close all
clc

global part_select;
part_select = 22; % 22 is gamma
global ener_thres; % keV
ener_thres = 10;

global bins;
bins = logspace(log10(10),log10(38000),512);

%% reading data

[yy,process_list]=load_g4_file('recorded_photon_data.txt');

ID = yy(:,3);
ener_list = yy(:,6);

%% making spectrum
[n_simu,xout] = make_spectrum(ener_list);
histogram('BinEdges',xout,'BinCounts',n_simu,'DisplayStyle','stairs','LineWidth',1.5);
hold on

% smoothed one
n_simu2=n_simu;
n_simu2(1:239) = smooth(n_simu2(1:239),0.05);
n_simu2(245:end) = smooth(n_simu2(245:end),0.010);

histogram('BinEdges',xout,'BinCounts',n_simu2,'DisplayStyle','stairs','LineWidth',1.5);

title({'no B-field', ...
    'source at 12 km altitude',...
    'Initial beaming : isotropic with 45^o half angle','Initial photon spec : exp(E/7.3 MeV)/E','record at 500 km'})
xlabel('Energy (keV)')
ylabel('dn/de spectrum (keV^{-1})')
grid on
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')


%% spectrum with creator process decomposition
figure

plot_spec_for_process(process_list,ener_list)

[n_simu,xout] = make_spectrum(ener_list);
histogram('BinEdges',xout,'BinCounts',n_simu,'DisplayStyle','stairs','LineWidth',1.5,'LineStyle',':','DisplayName','sum');

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

legend('show')
title({'no B-field', ...
    'source at 12 km altitude',...
    'Initial beaming : isotropic with 45^o half angle','Initial photon spec : exp(E/7.3 MeV)/E','record at 500 km','decomposed by creation process'})
xlabel('Energy (keV)')
ylabel('dn/de spectrum (keV^{-1})')
grid on








%%
function plot_spec_for_process(process_list,ener_list)
unique_process = unique(process_list);

for ii=1:length(unique_process)
    
    list_e = ener_list(strcmp(process_list, unique_process{ii}));
    
    %  'annihil'    'compt'    'conv'    'eBrem'    'eIoni'    'phot'
    [n_simu,xout] = make_spectrum(list_e);
    histogram('BinEdges',xout,'BinCounts',n_simu,'DisplayStyle','stairs','LineWidth',1.5,'DisplayName',unique_process{ii});
    hold on
end

end

%%
function [n_simu,xout] = make_spectrum(ener_list)
global bins
[n_simu,xout] = histcounts(ener_list,bins);
n_simu = n_simu ./ diff(bins);
end

%%
function [n_simu,xout] = make_rd_dist(rd_list)
bins = logspace(log10(min(rd_list)),log10(max(rd_list)),64);
size_rings = pi .* (bins(2:end)-bins(1:end-1)).^2;
[n_simu,xout] = histcounts(rd_list,bins);
n_simu = n_simu ./ size_rings;
end


%% can probably be optimized more
function [yy,process_list]=load_g4_file(filename)
global part_select;
global ener_thres;

yy0 = importdata(filename);
yy=[];

DD = regexp(yy0, ' ', 'split');
DD = vertcat(DD{:});
FF = DD(:,1:8);
yy = cellfun(@str2num,FF);

process_list0 = DD(:,9);

PDG_list = yy(:,4);
ener_list = yy(:,6);
radial_dist = yy(:,6);

% to_keep = rad_dist<200 & PDG_list==11 & ener_list>10;
to_keep = PDG_list==part_select & ener_list>ener_thres & radial_dist<300000 & ~(process_list0=="annihil" & ener_list>=512);

yy = yy(to_keep,:);
process_list=process_list0(to_keep);


end

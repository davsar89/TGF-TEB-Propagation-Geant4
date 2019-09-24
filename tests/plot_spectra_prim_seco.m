clearvars
close all
clc

%%
yy0 = load_data('./');
%%

for ii = 1:length(yy0)
    out = sscanf(yy0{ii},'%f',[1 7]);
    yy(ii,:) = out;
    chr = yy0{ii};
    [AA,~,~,nextindex] = sscanf(chr,'%s ',[1 7]);
    process_list0{ii} = chr(nextindex:end);
end


PDG_list = yy(:,1);
rad_dist = yy(:,5);
ener_list = yy(:,3);

to_keep = rad_dist<200 & PDG_list==22 & ener_list>10;
yy = yy(to_keep,:);

process_list_tmp={};
for ii=1:length(process_list0)
    if (to_keep(ii))
        process_list_tmp{end+1} = process_list0{ii};
    end
end

process_list = process_list_tmp;

ID = yy(:,7);
ener_list = yy(:,3);

[n_simu,xout] = make_spectrum(ener_list(ID==1));
histogram('BinEdges',xout,'BinCounts',n_simu,'DisplayStyle','stairs','LineWidth',1.5);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

hold on
[n_simu,xout] = make_spectrum(ener_list(ID>1));
histogram('BinEdges',xout,'BinCounts',n_simu,'DisplayStyle','stairs','LineWidth',1.5);

[n_simu,xout] = make_spectrum(ener_list);
histogram('BinEdges',xout,'BinCounts',n_simu,'DisplayStyle','stairs','LineWidth',1.5);

legend('primaries','secondaries','total')
title({'Photons only','no B-field', 'Record radial distance : between 0 and 200 km',...
    'source at 12 km altitude',...
    'Initial beaming : isotropic with 45^o half angle','Initial photon spec : exp(E/7.6 MeV)/E'})
xlabel('Energy (keV)')
ylabel('dn/de spectrum (keV^{-1})')
grid on

%% creator process decomposition for secondaries
figure

second_eners = ener_list(ID>1)';

process_list2={};
for ii=1:length(process_list)
    if (ID(ii)>1)
        process_list2{end+1} = process_list{ii};
    end
end

process_list = process_list2;

eBrem = second_eners(strcmp(process_list, 'eBrem'));
annihil = second_eners(strcmp(process_list, 'annihil'));

%  'annihil'    'compt'    'conv'    'eBrem'    'eIoni'    'phot'
[n_simu,xout] = make_spectrum(eBrem);
histogram('BinEdges',xout,'BinCounts',n_simu,'DisplayStyle','stairs','LineWidth',1.5);
hold on

[n_simu,xout] = make_spectrum(annihil);
histogram('BinEdges',xout,'BinCounts',n_simu,'DisplayStyle','stairs','LineWidth',1.5);

[n_simu,xout] = make_spectrum(ener_list(ID>1));
histogram('BinEdges',xout,'BinCounts',n_simu,'DisplayStyle','stairs','LineWidth',1.5,'LineStyle',':');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

legend('Bremsstrahlung','annihilation','all secondaries')
title({'Photons only','no B-field', 'Record radial distance : between 0 and 200 km',...
    'source at 12 km altitude',...
    'Initial beaming : isotropic with 45^o half angle','Initial photon spec : exp(E/7.6 MeV)/E'})
xlabel('Energy (keV)')
ylabel('dn/de spectrum (keV^{-1})')


%% radial distance distribution
rad_dist = rad_dist(PDG_list==22 & rad_dist>5);
figure
[n_simu1,xout1] = make_rd_dist(rad_dist);
histogram('BinEdges',xout1,'BinCounts',n_simu1,'DisplayStyle','stairs','LineWidth',1.5);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('Distance (km)')
ylabel('particle counts per area per bin size')
xlabel('Energy (keV)')
ylabel('dn/de spectrum (keV^{-1})')

grid on


%%
function [n_simu,xout] = make_spectrum(ener_list)
bins = logspace(log10(10),log10(32000),64);
[n_simu,xout] = histcounts(ener_list,bins);
n_simu = n_simu ./ diff(bins);
end

function yy=load_data(folder)

files = dir([folder 'detParticles_*.out']);
yy = [];
for ii=1:length(files)
    files(ii).name;
    yy = [yy;importdata([folder files(ii).name])];
end

end

function [n_simu,xout] = make_rd_dist(rd_list)
bins = logspace(log10(min(rd_list)),log10(max(rd_list)),64);
size_rings = pi .* (bins(2:end)-bins(1:end-1)).^2;
[n_simu,xout] = histcounts(rd_list,bins);
n_simu = n_simu ./ size_rings;
end

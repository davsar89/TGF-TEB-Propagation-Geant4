clear all
close all
clc

MIN_ENER = 100; % keV

hemisphere_to_keep = 0 ; % 0 for lower hemisphere and 1 for higher hemisphere

energy_grid = logspace(log10(MIN_ENER),log10(40000),32);

%% gather data

folder = './541_15_30_Gaussian_0/';

% Define the command to concatenate files
command = ['cat ' folder '* > fused.out'];
% Run the command
[status, cmdout] = system(command);
data = better_importdata('fused.out');

%% process

seeds = data(:,1);
event_nb = data(:,5);

nb_sampled_initial_source = get_nb_photons_sampled_total(seeds,event_nb);

disp(['Number of initialy sampled TGF photons : ' num2str(nb_sampled_initial_source)])

SOURCE_TGF_LAT = unique(data(:,21));
SOURCE_TGF_LONG = unique(data(:,22));

type = data(:,7);

data = data(type==11 | type==-11 ,:); % remove photons (if any)

type = data(:,7);
energies = data(:,9);
times = data(:,8)/1000.0; % us to ms

momentum_x = data(:,17); % normalized momentum i.e. vx^2+vy^2+vz^2 == 1
momentum_y = data(:,18);
momentum_z = data(:,19);

mom_norm = sqrt(momentum_x.^2+momentum_y.^2+momentum_z.^2);
momentum_x = momentum_x./mom_norm;
momentum_y = momentum_y./mom_norm;
momentum_z = momentum_z./mom_norm;
mom_norm = sqrt(momentum_x.^2+momentum_y.^2+momentum_z.^2);

nb_e = sum(type==11);
nb_p = sum(type==-11);

positrons_electron_ratio = nb_p / (nb_e + nb_p);

disp(['positrons-to-electrons ratio : ' num2str(positrons_electron_ratio)])

figure(1) % plot on map
% electrons
lat = data(:,11);
lat_e = data(type==11,11);
lon_e = data(type==11,12);
lat_p = data(type==-11,11);
lon_p = data(type==-11,12);

plot(lon_e,lat_e,'b+','DisplayName',['electrons']);
xlabel('long')
ylabel('lat')
hold on
plot(lon_p,lat_p,'r+','DisplayName',['positrons']);
grid on
legend('show', 'Location', 'northoutside')

plot(SOURCE_TGF_LONG,SOURCE_TGF_LAT,'+m','DisplayName',['source TGF'], 'MarkerSize', 10)

% find latitude separation between the thow hemispheres of the TEB (somethimes it is not simply 0)
middle_lat = find_separation_latitude(lat_e);
yline(middle_lat,'--','DisplayName','TEB hemisphere separator');
title('lat long map')

saveas(gcf,'map.png')

% filter omly one hemisphere

if hemisphere_to_keep==0
    tk = lat<middle_lat;
elseif hemisphere_to_keep==1
    tk = lat>middle_lat;
end


figure(2) % plot energy spectra

% electrons
[N_e,~] = histcounts(energies(type==11 & tk),energy_grid);
nde_e = N_e./diff(energy_grid);
histogram('BinEdges',energy_grid,'BinCounts',nde_e,'DisplayStyle','stairs','LineWidth',1.5,'DisplayName',['electrons'],'edgecolor','b');
grid on
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Energy (keV)')
ylabel('n/de spectrum (keV^{-1})')
title('energy spectra')
% positrons
hold on
[N_p,~] = histcounts(energies(type==-11 & tk),energy_grid);
nde_p = N_p./diff(energy_grid);
histogram('BinEdges',energy_grid,'BinCounts',nde_p,'DisplayStyle','stairs','LineWidth',1.5,'DisplayName',['positrons'],'edgecolor','r');
legend('show')
saveas(gcf,'energy_spectrum.png')

% lightcurves
figure(3)
time_grid = min(times):1:max(times);
[lc_e,~] = histcounts(times(type==11 & tk),time_grid);
[lc_p,~] = histcounts(times(type==-11 & tk),time_grid);

histogram('BinEdges',time_grid,'BinCounts',lc_e,'DisplayStyle','stairs','LineWidth',1.5,'edgecolor','b','DisplayName',['electrons']);
hold on
histogram('BinEdges',time_grid,'BinCounts',lc_p,'DisplayStyle','stairs','LineWidth',1.5,'edgecolor','r','DisplayName',['positrons']);
grid on
xlabel('time (ms)')
ylabel('counts per time bin')
legend('show')
title('lightcurve')
axis tight
saveas(gcf,'lightcurve.png')


%%
function data = better_importdata(filename)

% read first line to get the number of columns
fid = fopen(filename,'r');
your_text = fgetl(fid);
fclose(fid);
number_of_columns = length(str2num(your_text));

% read the full file
fileID = fopen(filename);
C = textscan(fileID, [repmat('%f ',[1 number_of_columns])],'CommentStyle','#');
fclose(fileID);

data=[C{:}];

end

%% find latitude separation between the thow hemispheres of the TEB (somethimes it is not simply 0)
function middle_x = find_separation_latitude(x)

% Sort the x data
x_sorted = sort(x);

% Find gaps in the sorted x data
gaps = diff(x_sorted);

% Find the index of the largest gap
[~, index] = max(gaps);

% Calculate the middle x value between the largest gap
middle_x = (x_sorted(index) + x_sorted(index + 1)) / 2;

end

% get the number of initial TGF photons rhat have been sampled (could be usefull to recover source brightness in the future)
function nb_sampled_initial_source = get_nb_photons_sampled_total(seed_nb,event_nb)

uniq_seed = unique(seed_nb);

nb_sampled_initial_source = 0;

for ii = 1:length(uniq_seed)
    nb_sampled_initial_source = nb_sampled_initial_source + max(event_nb(seed_nb==uniq_seed(ii)));
    
end

end

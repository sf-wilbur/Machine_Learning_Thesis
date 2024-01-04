clear; close all; clc;

% NEED TO MERGE THE SUM DATA locations into this...the raw locations are
% crap because they just come from the associator.

% dbstop if error

months = {'04','05','06','07','08','09','10','11','12'};
counts = {'1','2', '3','4', '5', '6', '7', '8', '9'};
% Build the individual files
for ii = 1:numel(months)

phase_file = sprintf('Y2000-2020%s.phs',months{ii});
fprintf('Processing %s\n', phase_file);

summary_file = sprintf('hypo%s.sum',counts{ii});
fprintf('Supplementing with %s\n', summary_file);

c = read_PHS_eqt( phase_file, summary_file );

[~,output_file,~] = fileparts( phase_file );

save( fullfile('PICK_DATA', output_file), 'c' );
clear c

end

%% Build a single files with all EQT events

clear; clc;

Y2000_files = dir('PICK_DATA/Y2000*');

% load the first month
file = fullfile( Y2000_files(1).folder, Y2000_files(1).name );
load( file );
S = c; % save a copy   

for ii = 2:numel(Y2000_files)
    file = fullfile( Y2000_files(ii).folder, Y2000_files(ii).name );
    load( file );
    S = [S, c]; % append along first dimension
end
save( 'PICK_DATA/eqt-all', 'S' );

%% Do some plotting

clear; clc;

load( 'PICK_DATA/eqt-all', 'S' );
c = S; clear S; % change the name

%%


clc;

X(:,1) = [c.lon];
X(:,2) = [c.lat];

% lon_limits = [min(X(:,1)), max(X(:,1))];
lon_limits = [-115.6, -114.6];
% lat_limits = [min(X(:,2)), max(X(:,2))];
lat_limits = [43.0, 44.8];

kill_idx = X(:,2) > lat_limits(2) & X(:,2) < lat_limits(1);
X(kill_idx,:) = []; % remove events outside of the lat_limits
% kill_idx = X(:,1) > lon_limits(2) & X(:,1) < lon_limits(1);
% X(kill_idx,:) = []; % remove events outside of the lon_limits

lon_Ctrs = linspace( lon_limits(1), lon_limits(2), 51);
lat_Ctrs = linspace( lat_limits(1), lat_limits(2), 51);

h = figure;
subplot(1,2,2)
hist3(X, 'Ctrs', {lon_Ctrs lat_Ctrs}, 'CdataMode', 'auto');
xlabel('Longitude [deg]');
colorbar;
view(2); axis('square');

subplot(1,2,1)
geoplot(X(:,2),X(:,1),'ok');
geobasemap('colorterrain');
% plot([c.lon],[c.lat],'ok');
% axis([lon_limits, lat_limits]);
geolimits(lat_limits, lon_limits);
% axis('square'); 
% grid on;
% xlabel('Longitude [deg]');
% ylabel('Latitude [deg]');

set(h,'color','w');
set( findall( h, '-property', 'FontSize' ), 'FontSize', 18 );
set( findall( h, '-property', 'FontName' ), 'FontName', 'Helvetica' );
set( findall( h, '-property', 'FontWeight' ), 'FontWeight', 'Bold' );
set( h, 'Position', [100 100 1500 1500] );
set( h, 'PaperPositionMode', 'auto' );

% print( h, '-dpng', 'M2p5_heatmap.png', '-r300');



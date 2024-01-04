%Load Sulphur EQ databases
% clear
min_lat=42.45;
max_lat=42.7;
min_lon=-111.6;
max_lon=-111.30;
 EQT_Sulphur_AK=readtable('EQ_Catalogs/EQT_Sulphur_ak135.csv'); %AK135 model
% EQT_Sulphur_AK=readtable('EQT_Sulphur_grad_vel.csv'); %gradient model



A_B_AK=[EQT_Sulphur_AK.Verr<=100 & EQT_Sulphur_AK.Herr<=100 ];
% A_B_AK=[1:length(EQT_Sulphur_AK.lat)];


EQT_Sulphur=readtable('EQ_Catalogs/eqt_Sulphur_sulvel.csv'); % local velocity model


A_B=[EQT_Sulphur.Verr<=5 & EQT_Sulphur.Herr<=5];

USGS_Sulphur=readtable('EQ_Catalogs/Sulphur_2017_2018_EQ_catalog.csv'); %USGS 2020 catalog
no_10km_idx=find(USGS_Sulphur.depth~=10&USGS_Sulphur.depthError<3); % remove EQ with default depth or large depth error

% p_model=readtable('sawtooth_p_model.txt');
stations=readtable('Station_Data/UU_Sulphur_Stations.csv');
USGS_mainshock=[42.6474	-111.4492]; %M5.3 mainshock, 9/2/2017
UU_mainshock=[42.651 -111.441];
USGS_aftershock=[42.563 -111.416]; %M5 aftershock, 9/10/2017
UU_aftershock=[42.569 -111.419];
Qfaults = shaperead('Shape_Files/QuaternaryFaults.shp', 'UseGeoCoords', true,'BoundingBox',[min_lon min_lat; max_lon max_lat]);
Qfaults2 = shaperead('Shape_Files/Qfaults_US_Database.shp', 'UseGeoCoords', true,'BoundingBox',[min_lon min_lat; max_lon max_lat]);

%%  bin aftershocks with step size defined above

j=1;
 step_size_lat=.00227;%This is 0.25 km latitude
 step_size_lon=.00357;%This is 0.25 km longitude
%   step_size_lat=.009;%This is 1 km latitude
%  step_size_lon=.014;%This is 1 km longitude
 circle_size=12;%for plotting with geoscatter
 
for lat=min_lat:step_size_lat:max_lat
    for lon=min_lon:step_size_lon:max_lon
        %AK135 EQT data
        depth_mean=mean(EQT_Sulphur_AK.depth_FSL(EQT_Sulphur_AK.lat(A_B_AK)<lat+step_size_lat&EQT_Sulphur_AK.lat(A_B_AK)>lat-step_size_lat&EQT_Sulphur_AK.lon(A_B_AK)>lon-step_size_lon&EQT_Sulphur_AK.lon(A_B_AK)<lon+step_size_lon));
        depth_std=std(EQT_Sulphur_AK.depth_FSL(EQT_Sulphur_AK.lat(A_B_AK)<lat+step_size_lat&EQT_Sulphur_AK.lat(A_B_AK)>lat-step_size_lat&EQT_Sulphur_AK.lon(A_B_AK)>lon-step_size_lon&EQT_Sulphur_AK.lon(A_B_AK)<lon+step_size_lon));
        depth_bin_size=length(EQT_Sulphur_AK.depth_FSL(EQT_Sulphur_AK.lat(A_B_AK)<lat+step_size_lat&EQT_Sulphur_AK.lat(A_B_AK)>lat-step_size_lat&EQT_Sulphur_AK.lon(A_B_AK)>lon-step_size_lon&EQT_Sulphur_AK.lon(A_B_AK)<lon+step_size_lon));
        output_depth_AK(j,:)=[depth_mean depth_std depth_bin_size lat lon];
        %local velocity model EQT data
        depth_mean=mean(EQT_Sulphur.depth_FSL(EQT_Sulphur.lat(A_B)<lat+step_size_lat&EQT_Sulphur.lat(A_B)>lat-step_size_lat&EQT_Sulphur.lon(A_B)>lon-step_size_lon&EQT_Sulphur.lon(A_B)<lon+step_size_lon));
        depth_std=std(EQT_Sulphur.depth_FSL(EQT_Sulphur.lat(A_B)<lat+step_size_lat&EQT_Sulphur.lat(A_B)>lat-step_size_lat&EQT_Sulphur.lon(A_B)>lon-step_size_lon&EQT_Sulphur.lon(A_B)<lon+step_size_lon));
        depth_bin_size=length(EQT_Sulphur.depth_FSL(EQT_Sulphur.lat(A_B)<lat+step_size_lat&EQT_Sulphur.lat(A_B)>lat-step_size_lat&EQT_Sulphur.lon(A_B)>lon-step_size_lon&EQT_Sulphur.lon(A_B)<lon+step_size_lon));
        output_depth(j,:)=[depth_mean depth_std depth_bin_size lat lon];
        % USGS data
        depth_mean=mean(USGS_Sulphur.depth(USGS_Sulphur.latitude(no_10km_idx)<lat+step_size_lat&USGS_Sulphur.latitude(no_10km_idx)>lat-step_size_lat&USGS_Sulphur.longitude(no_10km_idx)>lon-step_size_lon&USGS_Sulphur.longitude(no_10km_idx)<lon+step_size_lon));
        depth_std=std(USGS_Sulphur.depth(USGS_Sulphur.latitude(no_10km_idx)<lat+step_size_lat&USGS_Sulphur.latitude(no_10km_idx)>lat-step_size_lat&USGS_Sulphur.longitude(no_10km_idx)>lon-step_size_lon&USGS_Sulphur.longitude(no_10km_idx)<lon+step_size_lon));
        depth_bin_size=length(USGS_Sulphur.depth(USGS_Sulphur.latitude(no_10km_idx)<lat+step_size_lat&USGS_Sulphur.latitude(no_10km_idx)>lat-step_size_lat&USGS_Sulphur.longitude(no_10km_idx)>lon-step_size_lon&USGS_Sulphur.longitude(no_10km_idx)<lon+step_size_lon));
        output_depth_USGS(j,:)=[depth_mean depth_std depth_bin_size lat lon];
        j=j+1;
    end
end
output_depth_AK(output_depth_AK(:,3)<2,3)=NaN;%remove cells with fewer than 5 EQ's
output_depth_AK=output_depth_AK(~isnan(output_depth_AK(:,1)),:);
output_depth_AK=output_depth_AK(~isnan(output_depth_AK(:,3)),:);

output_depth(output_depth(:,3)<2,3)=NaN;%remove cells with fewer than 5 EQ's
output_depth=output_depth(~isnan(output_depth(:,1)),:);
output_depth=output_depth(~isnan(output_depth(:,3)),:);

output_depthUSGS(output_depth_USGS(:,3)<1,3)=NaN; %remove cells with fewer than 2 EQ's
output_depth_USGS=output_depth_USGS(~isnan(output_depth_USGS(:,1)),:);
output_depth_USGS=output_depth_USGS(~isnan(output_depth_USGS(:,3)),:);
%%
figure
%Plot AK135 EQ density map
subplot(332);
geoscatter(output_depth_AK(:,4),output_depth_AK(:,5),circle_size,output_depth_AK(:,3),'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)

%geoscatter(output_depth_AK(:,4),output_depth_AK(:,5),circle_size,output_depth_AK(:,1),'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
geobasemap streets
geolimits([min_lat max_lat],[min_lon max_lon])
set([get(gca(),'LatitudeLabel') get(gca(),'LongitudeLabel')],'FontSize',16)
% caxis([10 1000])
set(gcf,'DefaultTextFontSize',14); hold on
set(gcf, 'DefaultAxesFontSize',14);
colorbar; 
h=colorbar;
h.Label.String = 'Number of Earthquakes';
h.FontSize = 14;
caxis([10 80])
hold on
colormap('jet')

% for i=1:2
%     geoplot(Qfaults(i).Lat,Qfaults(i).Lon,'r')
%     geoplot(Qfaults2(i).Lat,Qfaults2(i).Lon,'k')
% end

geoplot(stations.Latitude,stations.Longitude,'kv')
% geoplot(USGS_mainshock(1),USGS_mainshock(2),'rp','MarkerSize',12)
geoplot(UU_mainshock(1),UU_mainshock(2),'rp','MarkerSize',12)
% geoplot(USGS_aftershock(1),USGS_aftershock(2),'gp','MarkerSize',12)
geoplot(UU_aftershock(1),UU_aftershock(2),'gp','MarkerSize',12)
title('EQT Earthquake Density - AK135f', 'FontSize',14)

% plot histogram of Sulphur earthquake depths with AK135 model
subplot(335);histogram(EQT_Sulphur_AK.depth_FSL(A_B_AK),'BinWidth',.5)
hold on
% histogram(EQT_SUlphur_AK.depth_FSL_FSL_FSL(A_B_AK),'BinWidth',.5)
xlim([-1 20])
xlabel('Depth (km)')
ylabel('Count')
title('Sulphur Aftershocks ( AK135f EQT)')
text(15,1500,['Mean: ', num2str(round(mean(EQT_Sulphur_AK.depth_FSL(A_B_AK)),1)),' (' num2str(round(std(EQT_Sulphur_AK.depth_FSL(A_B_AK)),1)) ')'])
% text(17,15,['Mean: ', num2str(round(mean(EQT_SUlphur_AK.depth_FSL(B_event_idx_AK)),1)),' (' num2str(round(std(EQT_SUlphur_AK.depth_FSL(B_event_idx_AK)),1)) ')'])
text(19,1500,['count=' num2str(length(EQT_Sulphur_AK.depth_FSL(A_B_AK)))])
view(90,90)

%plot histogram of Sulphur earthquake RMSE with AK135 velocity model
subplot(338);histogram(EQT_Sulphur_AK.rmse(A_B_AK),'BinWidth',.02)
hold on
% histogram(EQT_Sulphur_AK.rmse(A_event_idx_AK),'BinWidth',.02)
xlim([0 1])
xlabel('RMSE (sec)')
ylabel('Count')
title('EQT RMSE - AK135f Velocity Model')
text(.5,60,['Mean: ', num2str(round(mean(EQT_Sulphur_AK.rmse(A_B_AK)),2)),' (' num2str(round(std(EQT_Sulphur_AK.rmse(A_B_AK)),2)) ')'])
% text(.7,200,['Mean: ', num2str(round(mean(EQT_Sulphur_AK.rmse(B_event_idx_AK)),2)),' (' num2str(round(std(EQT_Sulphur_AK.rmse(B_event_idx_AK)),2)) ')'])
view(90,90)

% Plot EQ density map with local velocity model
subplot(331);geoscatter(output_depth(:,4),output_depth(:,5),circle_size,output_depth(:,3),'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
geobasemap streets
geolimits([min_lat max_lat],[min_lon max_lon])
% caxis([10 1000])
colorbar; 
h=colorbar;
h.Label.String = 'Number of Earthquakes';
h.FontSize = 14;
caxis([10 80])
hold on
% for i=1:2
%     geoplot(Qfaults(i).Lat,Qfaults(i).Lon,'r')
%     geoplot(Qfaults2(i).Lat,Qfaults2(i).Lon,'k')
% 
% end

geoplot(stations.Latitude,stations.Longitude,'kv')
% geoplot(USGS_mainshock(1),USGS_mainshock(2),'rp','MarkerSize',12)
geoplot(UU_mainshock(1),UU_mainshock(2),'rp','MarkerSize',12)
% geoplot(USGS_aftershock(1),USGS_aftershock(2),'gp','MarkerSize',12)
geoplot(UU_aftershock(1),UU_aftershock(2),'gp','MarkerSize',12)
title('EQT Earthquake Density - Local', 'FontSize',14)
set([get(gca(),'LatitudeLabel') get(gca(),'LongitudeLabel')],'FontSize',16)

% plot histogram of Sulphur earthquake depths with local velocity model
subplot(334);histogram(EQT_Sulphur.depth_FSL(A_B),'BinWidth',.5)
hold on
% histogram(EQT_Sulphur.depth(A_event_idx),'BinWidth',.5)
xlim([-1 20])
xlabel('Depth (km)')
ylabel('Count')
title('Sulphur Aftershocks (Local EQT)', 'FontSize',14)
text(15,150,['Mean: ', num2str(round(mean(EQT_Sulphur.depth_FSL(A_B)),1)),' (' num2str(round(std(EQT_Sulphur.depth_FSL(A_B)),1)) ')'])
% text(17,150,['Mean: ', num2str(round(mean(EQT_Sulphur.depth_FSL(B_event_idx)),1)),' (' num2str(round(std(EQT_Sulphur.depth(B_event_idx)),1)) ')'])
text(19,150,['count=' num2str(length(EQT_Sulphur.depth_FSL(A_B)))])
view(90,90)

%plot histogram of Sulphur earthquake RMSE with local velocity model
subplot(337);histogram(EQT_Sulphur.rmse(A_B),'BinWidth',.02)
hold on
% histogram(EQT_Sulphur.rmse(A_event_idx),'BinWidth',.02)
xlim([0 1])
xlabel('RMSE (sec)')
ylabel('Count')
title('EQT RMSE - Local Velocity Model')
text(.5,200,['Mean: ', num2str(round(mean(EQT_Sulphur.rmse(A_B)),2)),' (' num2str(round(std(EQT_Sulphur.rmse(A_B)),2)) ')'])
% text(.7,200,['Mean: ', num2str(round(mean(EQT_Sulphur.rmse(B_event_idx)),2)),' (' num2str(round(std(EQT_Sulphur.rmse(B_event_idx)),2)) ')'])
view(90,90)

% Plot USGS catalog EQ density
subplot(333);
geoscatter(output_depth_USGS(:,4),output_depth_USGS(:,5),circle_size,output_depth_USGS(:,3),'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
geobasemap streets
geolimits([min_lat max_lat],[min_lon max_lon])
set([get(gca(),'LatitudeLabel') get(gca(),'LongitudeLabel')],'FontSize',16)
colorbar; 
h=colorbar;
h.Label.String = 'Number of Earthquakes';
h.FontSize = 14;
caxis([10 80])
hold on
% for i=1:2
%     geoplot(Qfaults(i).Lat,Qfaults(i).Lon,'r')
%     geoplot(Qfaults2(i).Lat,Qfaults2(i).Lon,'k')
% 
% end

geoplot(stations.Latitude,stations.Longitude,'kv')
geoplot(USGS_mainshock(1),USGS_mainshock(2),'rp','MarkerSize',12)
% geoplot(UU_mainshock(1),UU_mainshock(2),'rp','MarkerSize',12)
geoplot(USGS_aftershock(1),USGS_aftershock(2),'gp','MarkerSize',12)
% geoplot(UU_aftershock(1),UU_aftershock(2),'gp','MarkerSize',12)
title('USGS Earthquake Density', 'FontSize',14)
% geoplot([min_lat min_lat max_lat max_lat min_lat],[min_lon max_lon max_lon min_lon min_lon],'linew',6)
%Plot depth distributions for USGS catalog
subplot(336);
hold on
histogram(USGS_Sulphur.depth(no_10km_idx),'BinWidth',.5)
xlim([-3.5 20])
xlabel('Depth (km)')
ylabel('Count')
title('Sulphur Aftershocks (USGS)', 'FontSize',14)
text(15,20,['Mean: ', num2str(round(mean(USGS_Sulphur.depth(no_10km_idx)),1)),...
    ' (' num2str(round(std(USGS_Sulphur.depth(no_10km_idx)),1)) ')'])
text(19,20,['count=' num2str(length(USGS_Sulphur.depth(no_10km_idx)))])
% text(13,40,num2str(round(mean(USGS_Sulphur.depth(no_10km_idx))),1))
% text(13,30,'Mean:')

view(90,90)

%plot RMS distribution for USGS catalog
subplot(339);
hold on
histogram(USGS_Sulphur.rms(no_10km_idx),'BinWidth',.02)
xlim([0 1])
xlabel('RMSE (sec)')
ylabel('Count')
title('USGS Aftershocks')
text(.5,20,['Mean: ', num2str(round(mean(USGS_Sulphur.rms(no_10km_idx)),2)),...
    ' (' num2str(round(std(USGS_Sulphur.rms(no_10km_idx)),2)) ')'])
view(90,90)


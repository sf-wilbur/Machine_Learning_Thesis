%Load Stanley EQ databases
% clear
EQT_Stanley_AK= readtable('stanley_ak135_EQT.csv'); %local model sawvel
A_event_idx_AK=find(strcmp('A',EQT_Stanley_AK.quality));
B_event_idx_AK=find(strcmp('B',EQT_Stanley_AK.quality));
C_event_idx_AK=find(strcmp('C',EQT_Stanley_AK.quality));
A_B_AK=[EQT_Stanley_AK.Verr<=5 & EQT_Stanley_AK.Herr <=5];
% 
EQT_Stanley=readtable('EQT_Stanley_Gradient.csv'); % local velocity model with new
% depths 
A_event_idx=find(strcmp('A',EQT_Stanley.quality));
B_event_idx=find(strcmp('B',EQT_Stanley.quality));
C_event_idx=find(strcmp('C',EQT_Stanley.quality));
A_B=[EQT_Stanley.Verr<=5 & EQT_Stanley.Herr <=5];

USGS_Stanley=readtable('Stanley_2020_EQ_catalog.csv'); %USGS 2020 catalog
no_10km_idx=find(USGS_Stanley.depth~=10&USGS_Stanley.depthError<3); % remove EQ with default depth or large depth error

% p_model=readtable('sawtooth_p_model.txt');
stations=readtable('BSU_Stations.csv');
USGS_mainshock=[44.465 -115.118];
Montana_mainshock=[44.401 -115.239];
Qfaults = shaperead('QuaternaryFaults.shp', 'UseGeoCoords', true,'BoundingBox',[-115.8 44; -114.7 44.75]);
Qfaults2 = shaperead('Qfaults_US_Database.shp', 'UseGeoCoords', true,'BoundingBox',[-115.8 44; -114.7 44.75]);
%%  bin aftershocks with step size defined above

j=1;
 step_size_lat=.0045;%This is 0.5 km latitude
 step_size_lon=.007;%This is 0.5 km longitude
 circle_size=12;%for plotting with geoscatter
 
for lat=44:step_size_lat:44.7
    for lon=-115.4:step_size_lon:-114.75
        %AK135 EQT data
        depth_mean=mean(EQT_Stanley_AK.depth(EQT_Stanley_AK.lat(A_B_AK)<lat+step_size_lat&EQT_Stanley_AK.lat(A_B_AK)>lat-step_size_lat&EQT_Stanley_AK.lon(A_B_AK)>lon-step_size_lon&EQT_Stanley_AK.lon(A_B_AK)<lon+step_size_lon));
        depth_std=std(EQT_Stanley_AK.depth(EQT_Stanley_AK.lat(A_B_AK)<lat+step_size_lat&EQT_Stanley_AK.lat(A_B_AK)>lat-step_size_lat&EQT_Stanley_AK.lon(A_B_AK)>lon-step_size_lon&EQT_Stanley_AK.lon(A_B_AK)<lon+step_size_lon));
        depth_bin_size=length(EQT_Stanley_AK.depth(EQT_Stanley_AK.lat(A_B_AK)<lat+step_size_lat&EQT_Stanley_AK.lat(A_B_AK)>lat-step_size_lat&EQT_Stanley_AK.lon(A_B_AK)>lon-step_size_lon&EQT_Stanley_AK.lon(A_B_AK)<lon+step_size_lon));
        output_depth_AK(j,:)=[depth_mean depth_std depth_bin_size lat lon];
        %local velocity model EQT data
        depth_mean=mean(EQT_Stanley.depth_FSL(EQT_Stanley.lat(A_B)<lat+step_size_lat&EQT_Stanley.lat(A_B)>lat-step_size_lat&EQT_Stanley.lon(A_B)>lon-step_size_lon&EQT_Stanley.lon(A_B)<lon+step_size_lon));
        depth_std=std(EQT_Stanley.depth_FSL(EQT_Stanley.lat(A_B)<lat+step_size_lat&EQT_Stanley.lat(A_B)>lat-step_size_lat&EQT_Stanley.lon(A_B)>lon-step_size_lon&EQT_Stanley.lon(A_B)<lon+step_size_lon));
        depth_bin_size=length(EQT_Stanley.depth_FSL(EQT_Stanley.lat(A_B)<lat+step_size_lat&EQT_Stanley.lat(A_B)>lat-step_size_lat&EQT_Stanley.lon(A_B)>lon-step_size_lon&EQT_Stanley.lon(A_B)<lon+step_size_lon));
        output_depth(j,:)=[depth_mean depth_std depth_bin_size lat lon];
        % USGS data
        depth_mean=mean(USGS_Stanley.depth(USGS_Stanley.latitude(no_10km_idx)<lat+step_size_lat&USGS_Stanley.latitude(no_10km_idx)>lat-step_size_lat&USGS_Stanley.longitude(no_10km_idx)>lon-step_size_lon&USGS_Stanley.longitude(no_10km_idx)<lon+step_size_lon));
        depth_std=std(USGS_Stanley.depth(USGS_Stanley.latitude(no_10km_idx)<lat+step_size_lat&USGS_Stanley.latitude(no_10km_idx)>lat-step_size_lat&USGS_Stanley.longitude(no_10km_idx)>lon-step_size_lon&USGS_Stanley.longitude(no_10km_idx)<lon+step_size_lon));
        depth_bin_size=length(USGS_Stanley.depth(USGS_Stanley.latitude(no_10km_idx)<lat+step_size_lat&USGS_Stanley.latitude(no_10km_idx)>lat-step_size_lat&USGS_Stanley.longitude(no_10km_idx)>lon-step_size_lon&USGS_Stanley.longitude(no_10km_idx)<lon+step_size_lon));
        output_depth_USGS(j,:)=[depth_mean depth_std depth_bin_size lat lon];
        j=j+1;
    end
end
output_depth_AK(output_depth_AK(:,3)<5,3)=NaN;%remove cells with fewer than 5 EQ's
output_depth_AK=output_depth_AK(~isnan(output_depth_AK(:,1)),:);
output_depth_AK=output_depth_AK(~isnan(output_depth_AK(:,3)),:);

output_depth(output_depth(:,3)<5,3)=NaN;%remove cells with fewer than 5 EQ's
output_depth=output_depth(~isnan(output_depth(:,1)),:);
output_depth=output_depth(~isnan(output_depth(:,3)),:);

output_depthUSGS(output_depth_USGS(:,3)<2,3)=NaN; %remove cells with fewer than 2 EQ's
output_depth_USGS=output_depth_USGS(~isnan(output_depth_USGS(:,1)),:);
output_depth_USGS=output_depth_USGS(~isnan(output_depth_USGS(:,3)),:);
%%
figure
%Plot AK135 EQ density map
subplot(331);geoscatter(output_depth_AK(:,4),output_depth_AK(:,5),circle_size,output_depth_AK(:,3),'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
geobasemap streets
geolimits([43.9 44.7],[-115.7 -114.6])
% caxis([10 1000])
caxis([10 150])
colorbar
hold on

for i=1:11
    geoplot(Qfaults(i).Lat,Qfaults(i).Lon,'r')
    geoplot(Qfaults2(i).Lat,Qfaults2(i).Lon,'k')
end
geoplot(stations.Latitude,stations.Longitude,'kv')
% geoplot(USGS_mainshock(1),USGS_mainshock(2),'rp','MarkerSize',12)
geoplot(Montana_mainshock(1),Montana_mainshock(2),'rp','MarkerSize',12)
title('EQT Earthquake Density - Local Velocity Model')

% plot histogram of Stanley earthquake depths with AK135 model
subplot(334);
hold on
histogram(EQT_Stanley_AK.depth_FSL(A_B_AK),'BinWidth',.5)
xlim([-1 20])
xlabel('Depth (km)')
ylabel('Count')
title('Stanley Aftershocks (Local EQT)')

text(17,10,['Mean: ', num2str(round(mean(EQT_Stanley_AK.depth_FSL(A_B_AK)),1)),' (' num2str(round(std(EQT_Stanley_AK.depth(A_B_AK)),1)) ')'])
text(19, 10, ['count= ' num2str(length(EQT_Stanley_AK.depth_FSL(A_B_AK)))])
view(90,90)

%plot histogram of Stanley earthquake RMSE with AK135 velocity model
subplot(337);
hold on
histogram(EQT_Stanley_AK.rmse(A_B_AK),'BinWidth',.02)
xlim([0 1])
xlabel('RMSE (km)')
ylabel('Count')
title('EQT RMSE - Local Velocity Model ')

text(.7,10,['Mean: ', num2str(round(nanmean(EQT_Stanley_AK.rmse(A_B_AK)),2)),' (' num2str(round(nanstd(EQT_Stanley_AK.rmse(A_B_AK)),2)) ')'])
view(90,90)

% Plot EQ density map with local velocity model
subplot(332);geoscatter(output_depth(:,4),output_depth(:,5),circle_size,output_depth(:,3),'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
geobasemap streets
geolimits([43.9 44.7],[-115.7 -114.6])
% caxis([10 1000])
caxis([10 150])
colorbar
hold on

for i=1:11
    geoplot(Qfaults(i).Lat,Qfaults(i).Lon,'r')
    geoplot(Qfaults2(i).Lat,Qfaults2(i).Lon,'k')
end
geoplot(stations.Latitude,stations.Longitude,'kv')
% geoplot(USGS_mainshock(1),USGS_mainshock(2),'rp','MarkerSize',12)
geoplot(Montana_mainshock(1),Montana_mainshock(2),'rp','MarkerSize',12)
title('EQT Earthquake Density - Gradient Velocity Model ')

% plot histogram of Stanley earthquake depths with local velocity model
subplot(335);
hold on
histogram(EQT_Stanley.depth_FSL(A_B),'BinWidth',.5)
xlim([-2.5 20])
xlabel('Depth (km)')
ylabel('Count')
title('Stanley Aftershocks (Gradient EQT)')

text(17,10,['Mean: ', num2str(round(nanmean(EQT_Stanley.depth_FSL(A_B)),1)),' (' num2str(round(nanstd(EQT_Stanley.depth_FSL(A_B)),1)) ')'])
text(19, 10, ['count=' num2str(length(EQT_Stanley.depth_FSL(A_B)))])
view(90,90)

%plot histogram of Stanley earthquake RMSE with local velocity model
subplot(338);
hold on
histogram(EQT_Stanley.rmse(A_B),'BinWidth',.02)
xlim([0 1])
xlabel('RMSE (km)')
ylabel('Count')
title('EQT RMSE - Gradient Velocity Model')
text(.7,10,['Mean: ', num2str(round(nanmean(EQT_Stanley.rmse(A_B)),2)),' (' num2str(round(nanstd(EQT_Stanley.rmse(A_B)),2)) ')'])
view(90,90)

% Plot USGS catalog EQ density
subplot(333);
geoscatter(output_depth_USGS(:,4),output_depth_USGS(:,5),circle_size,output_depth_USGS(:,3),'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
geobasemap streets
geolimits([43.9 44.7],[-115.7 -114.6])
caxis([2 10])
colorbar
hold on

for i=1:11
    geoplot(Qfaults(i).Lat,Qfaults(i).Lon,'r')
    geoplot(Qfaults2(i).Lat,Qfaults2(i).Lon,'k')
end

geoplot(stations.Latitude,stations.Longitude,'kv')
geoplot(USGS_mainshock(1),USGS_mainshock(2),'rp','MarkerSize',12)
% geoplot(Montana_mainshock(1),Montana_mainshock(2),'rp','MarkerSize',12)
title('USGS Earthquake Density')

%Plot depth distributions for USGS catalog
subplot(336);
hold on
histogram(USGS_Stanley.depth(no_10km_idx),'BinWidth',.5)
xlim([-3 20])
xlabel('Depth (km)')
ylabel('Count')
title('2020 Stanley Aftershocks (USGS)')
text(17,20,['Mean: ', num2str(round(mean(USGS_Stanley.depth(no_10km_idx)),1)),...
    ' (' num2str(round(std(USGS_Stanley.depth(no_10km_idx)),1)) ')'])
text(19, 20, ['count= ' num2str(length(USGS_Stanley.depth(no_10km_idx)))])
% text(13,40,num2str(round(mean(USGS_Stanley.depth(no_10km_idx))),1))
% text(13,30,'Mean:')

view(90,90)

%plot RMS distribution for USGS catalog
subplot(339);
hold on
histogram(USGS_Stanley.rms(no_10km_idx),'BinWidth',.02)
xlim([0 1])
xlabel('RMSE (km)')
ylabel('Count')
title('USGS Aftershocks')
text(.7,20,['Mean: ', num2str(round(mean(USGS_Stanley.rms(no_10km_idx)),2)),...
    ' (' num2str(round(std(USGS_Stanley.rms(no_10km_idx)),2)) ')'])


view(90,90)

% %% Plot the raw events 
% 
% 
% figure
% %Plot AK135 EQ density map
% colordata = categorical(EQT_Stanley_AK.depth_FSL(A_B_AK));
% subplot(331);
% geobubble(EQT_Stanley_AK.lat(A_B_AK),EQT_Stanley_AK.lon(A_B_AK),EQT_Stanley_AK.depth_FSL(A_B_AK), '.')
% geobasemap streets
% geolimits([43.9 44.7],[-115.7 -114.6])
% % caxis([10 1000])
% caxis([10 150])
% colorbar
% hold on
% 
% for i=1:11
%     geoplot(Qfaults(i).Lat,Qfaults(i).Lon,'r')
%     geoplot(Qfaults2(i).Lat,Qfaults2(i).Lon,'k')
% end
% geoplot(stations.Latitude,stations.Longitude,'kv')
% % geoplot(USGS_mainshock(1),USGS_mainshock(2),'rp','MarkerSize',12)
% geoplot(Montana_mainshock(1),Montana_mainshock(2),'rp','MarkerSize',12)
% title('EQT Earthquakes - AK135 Velocity Model')


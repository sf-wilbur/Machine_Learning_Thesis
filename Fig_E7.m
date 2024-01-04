%Load Challis EQ databases
clear all; 

%%
min_lat=44;
max_lat=44.9;
min_lon=-114.8;
max_lon=-113.4;

main_lat = 44.62;
main_lon = -114.33;
main_lat2 = 44.6;	
main_lon2 = -114.33;
main_lat3 =44.5074;	
main_lon3 = -114.112;
EQT_Challis_AK=readtable('EQ_Catalogs/EQT_Challis_chalvel.csv'); %davenport  model

A_B_AK=[EQT_Challis_AK.Herr<=5&EQT_Challis_AK.Verr<=5];

EQT_Challis=readtable('EQ_Catalogs/challis_ak135_EQT.csv'); % gradient velocity model new depths


A_B=[EQT_Challis.Herr<=5 &EQT_Challis.Verr<=5];

USGS_Challis=readtable('EQ_Catalogs/challis_usgs.csv'); %USGS 2020 catalog
no_10km_idx=find(USGS_Challis.depth~=10); % remove EQ with default depth no depth error available for these picks 

% p_model=readtable('sawtooth_p_model.txt');
stations=readtable('Station_Data/challis_local_stations.csv');
USGS_mainshock=[42.6474	-111.4492]; %M5.3 mainshock, 9/2/2017
UU_mainshock=[42.651 -111.441];
USGS_aftershock=[42.563 -111.416]; %M5 aftershock, 9/10/2017
UU_aftershock=[42.569 -111.419];
Qfaults = shaperead('Shape_Files/QuaternaryFaults.shp', 'UseGeoCoords', true,'BoundingBox',[min_lon min_lat; max_lon max_lat]);
Qfaults2 = shaperead('Shape_Files/Qfaults_US_Database.shp', 'UseGeoCoords', true,'BoundingBox',[min_lon min_lat; max_lon max_lat]);
%% bin aftershocks with step size defined above

j=1;
 step_size_lat=.00454;%This is 0.5 km latitude
 step_size_lon=.00714;%This is 0.5 km longitude
%   step_size_lat=.009;%This is 1 km latitude
%  step_size_lon=.014;%This is 1 km longitude
 circle_size=12;%for plotting with geoscatter
 output_depth_AK=zeros(length(min_lat:step_size_lat:max_lat)*(length(min_lon:step_size_lon:max_lon)),5);
 output_depth=zeros(length(min_lat:step_size_lat:max_lat)*(length(min_lon:step_size_lon:max_lon)),5);
 output_depth_USGS=zeros(length(min_lat:step_size_lat:max_lat)*(length(min_lon:step_size_lon:max_lon)),5);
 output_depth_OD=zeros(length(min_lat:step_size_lat:max_lat)*(length(min_lon:step_size_lon:max_lon)),5);
for lat=min_lat:step_size_lat:max_lat
    for lon=min_lon:step_size_lon:max_lon
        %AK135 EQT data
        tmp1=EQT_Challis_AK.depth(find(EQT_Challis_AK.lat(A_B_AK)<lat+step_size_lat&EQT_Challis_AK.lat(A_B_AK)>lat-step_size_lat&EQT_Challis_AK.lon(A_B_AK)>lon-step_size_lon&EQT_Challis_AK.lon(A_B_AK)<lon+step_size_lon));
        depth_mean=mean(tmp1);
        depth_std=std(tmp1);
        depth_bin_size=length(tmp1);
        output_depth_AK(j,:)=[depth_mean depth_std depth_bin_size lat lon];
        %local velocity model EQT data
        tmp2=EQT_Challis.depth_FSL(EQT_Challis.lat(A_B)<lat+step_size_lat&EQT_Challis.lat(A_B)>lat-step_size_lat&EQT_Challis.lon(A_B)>lon-step_size_lon&EQT_Challis.lon(A_B)<lon+step_size_lon);
        depth_mean=mean(tmp2);
        depth_std=std(tmp2);
        depth_bin_size=length(tmp2);
        output_depth(j,:)=[depth_mean depth_std depth_bin_size lat lon];
        % USGS data
        tmp3=USGS_Challis.depth(USGS_Challis.lat(no_10km_idx)<lat+step_size_lat&USGS_Challis.lat(no_10km_idx)>lat-step_size_lat&USGS_Challis.lon(no_10km_idx)>lon-step_size_lon&USGS_Challis.lon(no_10km_idx)<lon+step_size_lon);
        depth_mean=mean(tmp3);
        depth_std=std(tmp3);
        depth_bin_size=length(tmp3);
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


output_depth_USGS(output_depth_USGS(:,3)<2,3)=NaN; %remove cells with fewer than 2 EQ's
output_depth_USGS=output_depth_USGS(~isnan(output_depth_USGS(:,1)),:);
output_depth_USGS=output_depth_USGS(~isnan(output_depth_USGS(:,3)),:);

%%
h = figure;


%Plot davenport model or local velocity model EQ density map
subplot(331);
% output_depth_AK(output_depth_AK(:,3) <= 2, 1) = NaN;
geoscatter(output_depth_AK(:,4),output_depth_AK(:,5),circle_size,output_depth_AK(:,3),'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);

geobasemap streets
set(gcf,'DefaultTextFontSize',14); hold on
set(gcf, 'DefaultAxesFontSize',14);
%geoscatter(output_depth_AK(:,4),output_depth_AK(:,5),circle_size,output_depth_AK(:,1),'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
geolimits([min_lat max_lat],[min_lon max_lon])
set([get(gca(),'LatitudeLabel') get(gca(),'LongitudeLabel')],'FontSize',16)

colorbar; 
h=colorbar;
h.Label.String = 'Number of Earthquakes';
h.FontSize = 14;
caxis([1 50])
hold on
colormap('jet')


% for i=1:32
%     geoplot(Qfaults(i).Lat,Qfaults(i).Lon,'r')
%     geoplot(Qfaults2(i).Lat,Qfaults2(i).Lon,'k')
% end

geoplot(stations.lat,stations.lon,'kv')
% geoplot(USGS_mainshock(1),USGS_mainshock(2),'rp','MarkerSize',12)

geoplot(main_lat,main_lon,'rp','MarkerSize',12)
geoplot(main_lat2,main_lon2,'rp','MarkerSize',12)
geoplot(main_lat3,main_lon3,'rp','MarkerSize',12)
% geoplot(USGS_aftershock(1),USGS_aftershock(2),'gp','MarkerSize',12)
% geoplot(UU_aftershock(1),UU_aftershock(2),'gp','MarkerSize',12)
title('EQT Earthquake Density - Local', 'FontSize',14)


% plot histogram of Challis earthquake depths with AK135 model
subplot(334);histogram(EQT_Challis_AK.depth(A_B_AK),'BinWidth',.5);
hold on
% histogram(EQT_Challis_AK.depth_FSL(A_B_AK),'BinWidth',.5)
xlim([-1 20])
xlabel('Depth (km)')
ylabel('Count')
title('Challis aftershocks (Local EQT)')
text(15,600,['Mean: ', num2str(round(mean(EQT_Challis_AK.depth_FSL(A_B_AK)),1)),' (' num2str(round(std(EQT_Challis_AK.depth_FSL(A_B_AK)),1)) ')'])
% text(17,15,['Mean: ', num2str(round(mean(EQT_Challis_AK.depth(B_event_idx_AK)),1)),' (' num2str(round(std(EQT_Challis_AK.depth(B_event_idx_AK)),1)) ')'])
text(19,600,['count=' num2str(length(EQT_Challis_AK.depth_FSL(A_B_AK)))])
view(90,90)

%plot histogram of Challis earthquake RMSE with dav/local velocity model
subplot(337);histogram(EQT_Challis_AK.rmse(A_B_AK),'BinWidth',.02);
hold on
% histogram(EQT_Challis_AK.rmse(A_event_idx_AK),'BinWidth',.02)
xlim([0 1])
xlabel('RMSE (sec)')
ylabel('Count')
title('EQT RMSE - Local Velocity Model')
text(.5,100,['Mean: ', num2str(round(mean(EQT_Challis_AK.rmse(A_B_AK)),2)),' (' num2str(round(std(EQT_Challis_AK.rmse(A_B_AK)),2)) ')'])
% text(.7,200,['Mean: ', num2str(round(mean(EQT_Challis_AK.rmse(B_event_idx_AK)),2)),' (' num2str(round(std(EQT_Challis_AK.rmse(B_event_idx_AK)),2)) ')'])
view(90,90)

% Plot EQ density map with local velocity model
subplot(332);
% output_depth(output_depth(:,3) <= 2, 1) = NaN;
geoscatter(output_depth(:,4),output_depth(:,5),circle_size,output_depth(:,3),'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
geobasemap streets

geolimits([min_lat max_lat],[min_lon max_lon])
set([get(gca(),'LatitudeLabel') get(gca(),'LongitudeLabel')],'FontSize',16)
% caxis([10 1000])
colorbar; 
h=colorbar;
h.Label.String = 'Number of Earthquakes';
h.FontSize = 14;
caxis([1 50])
hold on
% for i=1:32
%     geoplot(Qfaults(i).Lat,Qfaults(i).Lon,'r')
%     geoplot(Qfaults2(i).Lat,Qfaults2(i).Lon,'k')
% 
% end
geoplot(main_lat,main_lon,'rp','MarkerSize',12)
geoplot(main_lat2,main_lon2,'rp','MarkerSize',12)
geoplot(main_lat3,main_lon3,'rp','MarkerSize',12)
geoplot(stations.lat,stations.lon,'kv')
% geoplot(USGS_mainshock(1),USGS_mainshock(2),'rp','MarkerSize',12)
% geoplot(UU_mainshock(1),UU_mainshock(2),'rp','MarkerSize',12)
% geoplot(USGS_aftershock(1),USGS_aftershock(2),'gp','MarkerSize',12)
% geoplot(UU_aftershock(1),UU_aftershock(2),'gp','MarkerSize',12)
title('EQT Earthquake Density - AK135f','FontSize',14)
% title('EQT earthquake density - Gradient velocity model')
% plot histogram of Challis earthquake depths with gradient velocity model
subplot(335);histogram(EQT_Challis.depth_FSL(A_B),'BinWidth',.5)
hold on
% histogram(EQT_Challis.depth_FSL(A_event_idx),'BinWidth',.5)
xlim([-2.5 20])
xlabel('Depth (km)')
ylabel('Count')
title('Challis Aftershocks (AK135f EQT)')
text(15,100,['Mean: ', num2str(round(mean(EQT_Challis.depth(A_B)),1)),' (' num2str(round(std(EQT_Challis.depth(A_B)),1)) ')'])
% text(17,150,['Mean: ', num2str(round(mean(EQT_Challis.depth(B_event_idx)),1)),' (' num2str(round(std(EQT_Challis.depth(B_event_idx)),1)) ')'])
text(19,100,['count=' num2str(length(EQT_Challis.depth(A_B)))])
% text(15,80,['Mean: ', num2str(round(mean(EQT_Challis.depth(A_B)),1)),' (' num2str(round(std(EQT_Challis.depth(A_B)),1)) ')'])
% text(19,80,['count=' num2str(length(EQT_Challis.depth(A_B)))])
view(90,90)


%plot histogram of Challis earthquake RMSE with gradient velocity model
subplot(338);histogram(EQT_Challis.rmse(A_B),'BinWidth',.02)
hold on
% histogram(EQT_Challis.rmse(A_event_idx),'BinWidth',.02)
xlim([0 1])
xlabel('RMSE (sec)')
ylabel('Count')
% title('EQT RMSE - Gradient velocity model')
title('EQT RMSE - AK135f velocity model')
text(.5,40,['Mean: ', num2str(round(mean(EQT_Challis.rmse(A_B)),2)),' (' num2str(round(std(EQT_Challis.rmse(A_B)),2)) ')'])
% text(.7,200,['Mean: ', num2str(round(mean(EQT_Challis.rmse(B_event_idx)),2)),' (' num2str(round(std(EQT_Challis.rmse(B_event_idx)),2)) ')'])
% text(.5,45,['Mean: ', num2str(round(mean(EQT_Challis.rmse(A_B)),2)),' (' num2str(round(std(EQT_Challis.rmse(A_B)),2)) ')'])
view(90,90)

% Plot USGS catalog EQ density
subplot(333);
geoscatter(output_depth_USGS(:,4),output_depth_USGS(:,5),circle_size,output_depth_USGS(:,3),'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
geobasemap streets
geolimits([min_lat max_lat],[min_lon max_lon])
set([get(gca(),'LatitudeLabel') get(gca(),'LongitudeLabel')],'FontSize',16)
% caxis([2 10])
colorbar; 
h=colorbar;
h.Label.String = 'Number of Earthquakes';
h.FontSize = 14;
caxis([1 50])
hold on
% for i=1:32
%     geoplot(Qfaults(i).Lat,Qfaults(i).Lon,'r')
%     geoplot(Qfaults2(i).Lat,Qfaults2(i).Lon,'k')
% 
% end
geoplot(main_lat,main_lon,'rp','MarkerSize',12)
geoplot(main_lat2,main_lon2,'rp','MarkerSize',12)
geoplot(main_lat3,main_lon3,'rp','MarkerSize',12)
geoplot(stations.lat,stations.lon,'kv')
% geoplot(USGS_mainshock(1),USGS_mainshock(2),'rp','MarkerSize',12)
% geoplot(UU_mainshock(1),UU_mainshock(2),'rp','MarkerSize',12)
% geoplot(USGS_aftershock(1),USGS_aftershock(2),'gp','MarkerSize',12)
% geoplot(UU_aftershock(1),UU_aftershock(2),'gp','MarkerSize',12)
title('USGS earthquake density','FontSize',14)
% geoplot([min_lat min_lat max_lat max_lat min_lat],[min_lon max_lon max_lon min_lon min_lon],'linew',6)
%Plot depth distributions for USGS catalog
subplot(336);
hold on
histogram(USGS_Challis.depth(no_10km_idx),'BinWidth',.5)
xlim([-1 20])
xlabel('Depth (km)')
ylabel('Count')
title('Challis Aftershocks (USGS)')
text(15,15,['Mean: ', num2str(round(mean(USGS_Challis.depth(no_10km_idx)),1)),...
    ' (' num2str(round(std(USGS_Challis.depth(no_10km_idx)),1)) ')'])
text(19,15, ['count=' num2str(length(USGS_Challis.depth(no_10km_idx)))])
% text(13,40,num2str(round(mean(EQT_Challis.depth(no_10km_idx))),1))
% text(13,30,'Mean:')

view(90,90)

%plot RMS distribution for USGS catalog
subplot(339);
hold on
histogram(USGS_Challis.rms(no_10km_idx),'BinWidth',.02)
xlim([0 1])
xlabel('RMSE (sec)')
ylabel('Count')
title('USGS Aftershocks')
text(.5,6,['Mean: ', num2str(round(mean(USGS_Challis.rms(no_10km_idx)),2)),...
    ' (' num2str(round(std(USGS_Challis.rms(no_10km_idx)),2)) ')'])
view(90,90)



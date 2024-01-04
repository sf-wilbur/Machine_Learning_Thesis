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

EQT_Challis=readtable('EQ_Catalogs/EQT_Challis_Gradient.csv'); % gradient velocity model new depths


A_B=[EQT_Challis.Herr<=5 &EQT_Challis.Verr<=5];


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
 
 output_depth=zeros(length(min_lat:step_size_lat:max_lat)*(length(min_lon:step_size_lon:max_lon)),5);
 output_depth_OD=zeros(length(min_lat:step_size_lat:max_lat)*(length(min_lon:step_size_lon:max_lon)),5);
for lat=min_lat:step_size_lat:max_lat
    for lon=min_lon:step_size_lon:max_lon
        %AK135 EQT data

        %local velocity model EQT data
        tmp2=EQT_Challis.depth_FSL(EQT_Challis.lat(A_B)<lat+step_size_lat&EQT_Challis.lat(A_B)>lat-step_size_lat&EQT_Challis.lon(A_B)>lon-step_size_lon&EQT_Challis.lon(A_B)<lon+step_size_lon);
        depth_mean=mean(tmp2);
        depth_std=std(tmp2);
        depth_bin_size=length(tmp2);
        output_depth(j,:)=[depth_mean depth_std depth_bin_size lat lon];
  
        j=j+1;
        
    end
end


output_depth(output_depth(:,3)<2,3)=NaN;%remove cells with fewer than 5 EQ's
output_depth=output_depth(~isnan(output_depth(:,1)),:);
output_depth=output_depth(~isnan(output_depth(:,3)),:);

%%
figure
% Plot EQ density map with local velocity model
subplot(3,2,1:4);
% output_depth(output_depth(:,3) <= 2, 1) = NaN;
geoscatter(output_depth(:,4),output_depth(:,5),circle_size,output_depth(:,3),'filled','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5)
geobasemap streets

geolimits([min_lat max_lat],[min_lon max_lon])
% caxis([10 1000])
colorbar; 
h=colorbar;
h.Label.String = 'Number of Earthquakes';
caxis([1 50])

hold on
Qfaults3 = shaperead('Shape_Files/Geologic_Map_Idaho_Faults.shp', 'UseGeoCoords', true,'BoundingBox',[-115.75 44; -114.7 44.7]);
for i=8:90
    geoplot(Qfaults3(i).Lat,Qfaults3(i).Lon,'k')
end
for i=1:32
    geoplot(Qfaults(i).Lat,Qfaults(i).Lon,'r')
    geoplot(Qfaults2(i).Lat,Qfaults2(i).Lon,'k')
end

geoplot(main_lat,main_lon,'rp','MarkerSize',12)
geoplot(main_lat2,main_lon2,'rp','MarkerSize',12)
geoplot(main_lat3,main_lon3,'rp','MarkerSize',12)
geoplot(stations.lat,stations.lon,'kv')
% geoplot(USGS_mainshock(1),USGS_mainshock(2),'rp','MarkerSize',12)
% geoplot(UU_mainshock(1),UU_mainshock(2),'rp','MarkerSize',12)
% geoplot(USGS_aftershock(1),USGS_aftershock(2),'gp','MarkerSize',12)
% geoplot(UU_aftershock(1),UU_aftershock(2),'gp','MarkerSize',12)
title('EQT Earthquake Density - Gradient Velocity Model')
% title('EQT earthquake density - Gradient velocity model')
% plot histogram of Challis earthquake depths with gradient velocity model
subplot(3,2,5);histogram(EQT_Challis.depth_FSL(A_B),'BinWidth',.5)
hold on
% histogram(EQT_Challis.depth_FSL(A_event_idx),'BinWidth',.5)
xlim([-2.5 20])
xlabel('Depth (km)')
ylabel('Count')
title('Challis Aftershocks (Gradient EQT)')
text(15,150,['Mean: ', num2str(round(mean(EQT_Challis.depth(A_B)),1)),' (' num2str(round(std(EQT_Challis.depth(A_B)),1)) ')'])
% text(17,150,['Mean: ', num2str(round(mean(EQT_Challis.depth(B_event_idx)),1)),' (' num2str(round(std(EQT_Challis.depth(B_event_idx)),1)) ')'])
text(19,150,['count=' num2str(length(EQT_Challis.depth(A_B)))])
% text(15,80,['Mean: ', num2str(round(mean(EQT_Challis.depth(A_B)),1)),' (' num2str(round(std(EQT_Challis.depth(A_B)),1)) ')'])
% text(19,80,['count=' num2str(length(EQT_Challis.depth(A_B)))])
view(90,90)


%plot histogram of Challis earthquake RMSE with gradient velocity model
subplot(3,2,6);histogram(EQT_Challis.rmse(A_B),'BinWidth',.02)
hold on
% histogram(EQT_Challis.rmse(A_event_idx),'BinWidth',.02)
xlim([0 1])
xlabel('RMSE (sec)')
ylabel('Count')
% title('EQT RMSE - Gradient velocity model')
title('EQT RMSE - Gradient velocity model')
text(.5,200,['Mean: ', num2str(round(mean(EQT_Challis.rmse(A_B)),2)),' (' num2str(round(std(EQT_Challis.rmse(A_B)),2)) ')'])
% text(.7,200,['Mean: ', num2str(round(mean(EQT_Challis.rmse(B_event_idx)),2)),' (' num2str(round(std(EQT_Challis.rmse(B_event_idx)),2)) ')'])
% text(.5,45,['Mean: ', num2str(round(mean(EQT_Challis.rmse(A_B)),2)),' (' num2str(round(std(EQT_Challis.rmse(A_B)),2)) ')'])
view(90,90)
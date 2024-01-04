 clear;
 %%
% Read downloaded catalog 
 data=readtable('EQ_Catalogs/EQT_Stanley_QE.csv');
T= readtable("1970_present_EQ_westernUS_M2.5.xlsx"); %read in USGS catalog 
m_t = readtable("EQ_Catalogs/Mom_Tens_EQT_Locations.csv");
Vert_err=1;
Hor_err=.5;
max_depth=20;
USGS_mainshock=[44.465 -115.118];
Montana_mainshock=[44.401 -115.239];
idx =find(data.depth_FSL<max_depth&data.lon<-114.9&data.lon>-115.4&data.Herr<Hor_err&data.Verr<Vert_err);
%% change time format in USGS catalog 

formatIn = 'yyyy-mm-ddTHH:MM:SS'; %tell matlab the format for time in the catalog T
t1 = T.time(:);
t1 = datenum(t1, formatIn); %Get the correct format to compare the times between arrays 
T.otime = t1;

%% STANLEY create table
%for loop for horizontal error
T2 = height(T);
St = array2table(zeros(length(T2),7));
St.Properties.VariableNames = {'otime','lon', ...
           'lat','depth', 'horErr', 'rms', 'mag'};
for i = 1:T2
    if i>49630 && i< 52643 && 45 >= T.latitude(i)... 
        && T.latitude(i) >= 43.5 && T.longitude(i) >= -115.5 && T.longitude(i) <= -114 
        main_shock = find(T.mag == 6.5);
        St.otime(i) = T.otime(i);
        St.lat(i) = T.latitude(i);
        St.lon(i) = T.longitude(i);
        St.depth(i) = T.depth(i);
        St.horErr(i) = T.horizontalError(i);
        St.rms(i) = T.rms(i); 
        St.mag(i) = T.mag(i);
     
        %quality 
        
    end
    
end 

St(~St.otime,:) = []; % deletes the rows where otime is 0 and save the table as C 

%% MATCH THE USGS EVENTS TO EQT for Stanley
clear dt dist val ids counter ;
% Start_date = {'31-March-2020 00:00:00'};
Start_date = datetime(T.otime(main_shock(1)), 'ConvertFrom','datenum');
for x = 1:height(data) % create a loop that finds time difference since the mainshock on March 31st 
  days(x) = daysdif(Start_date, data.otime(x));
end 
data.days = days'; 
 %Change the format to math the sum files and compare the times
E = referenceEllipsoid('Earth'); % reference ellipse [m] for distance calculation
counter = 0;
t1 = St.otime(:,1); % USGS origin times since mainshock 
% t2_st = zeros(height(data),1);
t2_st = data.otime(idx) ;                  %t2 = inl(k).otime; for structure

for ii = 1 : numel(t1) % loop through eqt events CHANGE TO T1
    
    % find the USGS event that is closest in time to eqt_event(i)
    t_diff  = t2_st - t1(ii); 
    
    [val, ids(ii)] = min( abs(t_diff) ); % this lets them be negative 

    % Display some useful information to the user
    
  
    

    fprintf('EQT: %s, USGS: %s\n', datestr(t2_st(ids(ii))), datestr(t1(ii)));
    
    dt(ii) = val * 24*3600; % [s] convert time difference to seconds
%      dt_nm(i) = dt(i);
    if dt(ii) > 10
        counter = counter +1 ;
          
 

    else 

                % compute distance in meters
       St.latnew(ii) = data.lat(idx(ids(ii))); 
       St.lonnew(ii) = data.lon(idx(ids(ii))); 
     
       [arclen(ii), az(ii)]  = distance( St.lat((ii)), St.lon((ii)),...  % right now this is set to run for A events 
        data.lat(ids(ii)), data.lon(ids(ii)), E) ; %CHANGE CATALOGS HERE 
        data.mag(idx(ids(ii)))= St.mag(ii);
    end 
   % Append mag values from USGS to Stanl
end

t1_m = m_t.otime; % moment origin times
% t2_st = zeros(height(data),1);
                %t2 = inl(k).otime; for structure
counterm = 0;
for i = 1 : numel(t1_m) % loop through eqt events CHANGE TO T1
    
    % find the USGS event that is closest in time to eqt_event(i)
    t_diffm  = t2_st - t1_m(i); 
    
    [val, idi(i)] = min( abs(t_diffm) ); % this lets them be negative 

    % Display some useful information to the user
    
    fprintf('EQT: %s, Moment Tensors: %s\n', datestr(t2_st(idi(i))), datestr(t1_m(i)));
    
    dtm(i) = val * 24*3600; % [s] convert time difference to seconds
%      dt_nm(i) = dt(i);
     if dtm(i) > 5
        counterm = counterm +1 ;
        mt.lat(i) = nan;
        mt.lon(i) = nan; 
     else 
  
       m_t.lat(i) = data.lat(idx(idi(i))); 
       m_t.lon(i) = data.lon(idx(idi(i))); 
       m_t.her(i) = data.Herr(idx(idi(i)));
       m_t.ver(i) = data.Verr(idx(idi(i)));
     end
end
%% Replace lat and lon coord in moment tensor table with locations from EQT
m_t.Latitude(i) = data.lat(idi(i)); 
m_t.Longitude(i) = data.lon(idi(i));
%% convert EQ catalog to UTM
% find group that meets criteria
% set cross section angle
% rotate coordinate system to rotation angle
[EQ_UTM(:,1),EQ_UTM(:,2),utmzone] = deg2utm(data.lat,data.lon); %deg2utm is a matlab function in the current directory
[USGS_UTM(:,1),USGS_UTM(:,2),utmzone] = deg2utm(St.lat,St.lon);

groupA=find(data.depth_FSL<max_depth&data.lon<-114.9&data.lon>-115.4&data.Herr<Hor_err&data.Verr<Vert_err);
EQ_count=length(groupA);
groupM = find(data.mag>0);
group1=find(data.lat>=44.4&data.lat<44.49&data.depth_FSL<max_depth&data.lon<-115.15&data.lon>-115.25&data.Herr<Hor_err&data.Verr<Vert_err);
% group1=find(data.lat>=44.3&data.lat<44.33&data.depth_FSL<max_depth&data.lon<-115.05&data.lon>-115.25&data.Herr<Hor_err&data.Verr<Vert_err);
group_USGS2=find(St.lat>=44.4&St.lat<44.49&St.depth<max_depth&St.lon<-115.15&St.lon>-115.25);

group2=find(data.lat<44.3&data.lat>44.26&data.depth_FSL<max_depth&data.lon<-115.045&data.lon>-115.09&data.Herr<Hor_err&data.Verr<Vert_err);
group3=find(data.lat<44.16&data.depth_FSL<25&data.lon<-115.0&data.lon>-115.15&data.Herr<Hor_err&data.Verr<Vert_err);
group4=find(data.lat<=44.27&data.lat>44.2&data.depth_FSL<max_depth&data.lon<-114.92&data.lon>-115.22&data.Herr<Hor_err&data.Verr<Vert_err);

%calculate best fit UTM linear trend to group1 values
p1a=polyfit(EQ_UTM(group1,2),EQ_UTM(group1,1),1);
f1a = polyval(p1a,EQ_UTM(group1,2));

theta1=p1a(1)*100; %rotation angle from vertical (e.g., -10 = N10W)
 theta1=0;
R1=[cosd(theta1) -sind(theta1); sind(theta1) cosd(theta1)];
rotcoord1=EQ_UTM(:,:)*R1'; %find values near Stanley

%calculate best fit depth linear trend to group1
p1b=polyfit(data.depth_FSL(group1),rotcoord1(group1,1)/1000,1);
f1b = polyval(p1b,data.depth_FSL(group1));
dip1=p1b(1)*100+90;

%calculate best fit UTM linear trend to group2 values
p2a=polyfit(EQ_UTM(group2,2),EQ_UTM(group2,1),1);

f2a = polyval(p2a,EQ_UTM(group2,2));

theta2=p2a(1)*100; %rotation angle from vertical (e.g., -10 = N10W)
% theta2=0;
R2=[cosd(theta2) -sind(theta2); sind(theta2) cosd(theta2)];
rotcoord2=EQ_UTM(:,:)*R2'; %find values near Stanley

%calculate best fit depth linear trend to group2
p2b=polyfit(data.depth_FSL(group2),rotcoord2(group2,1)/1000,1);
f2b = polyval(p2b,data.depth_FSL(group2));
dip2=p2b(1)*100+90;

%calculate best fit UTM linear trend to group3 values
p3a=polyfit(EQ_UTM(group3,2),EQ_UTM(group3,1),1);
f3a = polyval(p3a,EQ_UTM(group3,2));

theta3=p3a(1)*100; %rotation angle from vertical (e.g., -10 = N10W)
R3=[cosd(theta3) -sind(theta3); sind(theta3) cosd(theta3)];
rotcoord3=EQ_UTM(:,:)*R3'; %find values near Stanley

%calculate best fit depth linear trend to group3
p3b=polyfit(data.depth_FSL(group3),rotcoord2(group3,1)/1000,1);
f3b = polyval(p3b,data.depth_FSL(group3));
dip3=p3b(1)*100+90;

%calculate best fit UTM linear trend to group4 values
p4a=polyfit(EQ_UTM(group4,2),EQ_UTM(group4,1),1);
f4a = polyval(p4a,EQ_UTM(group4,2));

theta4=p4a(1)*100; %rotation angle from vertical (e.g., -10 = N10W)
 theta4=-20;
R4=[cosd(theta4) -sind(theta4); sind(theta4) cosd(theta4)];
rotcoord4=EQ_UTM(:,:)*R4'; %find values near Stanley

%calculate best fit depth linear trend to group4
p4b=polyfit(data.depth_FSL(group4),rotcoord2(group4,1)/1000,1);
f4b = polyval(p4b,data.depth_FSL(group4));
dip4=p4b(1)*100+90;
%% 
figure;
%plot all data in lat/long
% subplot(2,4,1:2);
% plot(data.lon,data.lat,'bo')
% hold on
%plot group1 data in lat/long
Qfaults = shaperead('Shape_Files/QuaternaryFaults.shp', 'UseGeoCoords', true,'BoundingBox',[-115.75 44; -114.7 44.7]);

Qfaults2 = shaperead('Shape_Files/Geologic_Map_Idaho_Faults.shp', 'UseGeoCoords', true,'BoundingBox',[-115.75 44; -114.7 44.7]);

geoscatter(data.lat(groupA),data.lon(groupA),5,data.depth_FSL(groupA),'filled'); % color coded by depth 
% geoscatter(data.lat(groupA),data.lon(groupA),5,data.timeM(groupA),'filled'); % color coded by origin time
hold on
mdx = find(St.mag>=4);
geoscatter(St.latnew(mdx),St.lonnew(mdx), exp(St.mag(mdx)),'filled', 'o','MarkerFaceColor','k')
% mdx = find(m_t.Magnitude>=3.5);
% geoscatter(m_t.lat(mdx),m_t.lon(mdx), exp(m_t.Magnitude(mdx)),'filled', 'o','MarkerFaceColor','k')
hold on
for i=8:90
    geoplot(Qfaults2(i).Lat,Qfaults2(i).Lon,'k')
end

for i = 8:11
    geoplot(Qfaults(i).Lat,Qfaults(i).Lon,'r','linew',4) %Sawtooth fault
%     geoplot(Qfaults(1).Lat,Qfaults(1).Lon,'r') % Deadwood fault
%     geoplot(Qfaults2(1).Lat,Qfaults2(1).Lon,'k')
end 
% baseURL = "https://basemap.nationalmap.gov/ArcGIS/rest/services";
% usgsURL = baseURL + "/USGSTopo/MapServer/tile/${z}/${y}/${x}";
% basemap="USGSShadedReliefOnly";
geobasemap('topographic')
geoplot(44.475216, -115.142315,'kp', 'MarkerSize', 16,'MarkerFaceColor','blue'  )
geoplot(44.2155, -114.9352,'rp','Markersize',20)
geoplot(USGS_mainshock(1),USGS_mainshock(2), 'kp','MarkerSize',12, 'MarkerFaceColor','auto')
geoplot(Montana_mainshock(1),Montana_mainshock(2), 'kp','MarkerSize',12,'MarkerFaceColor','auto')
% geoplot(44.4573,-115.0980, 'kp','MarkerSize',12,'MarkerFaceColor','g')
% geoplot([44.41 44.41 44.49 44.49 44.41],[-115.25 -115.15 -115.15 -115.25 -115.25],'r-','linew',2)
% % geoplot([44.3 44.3 44.33 44.33 44.3],[-115.25 -115.05 -115.05 -115.25 -115.25],'r-','linew',2)
% bounding boxes
% geoplot([44.3 44.3 44.26 44.26 44.3],[-115.09 -115.045 -115.045 -115.09 -115.09],'k-','linew',2)
% geoplot([44.08 44.08 44.15 44.15 44.08],[-115. -115.15 -115.15 -115. -115.],'k-','linew',2)
% geoplot([44.2 44.2 44.27 44.27 44.2],[-115.22 -114.92 -114.92 -115.22 -115.22],'k-','linew',2)

% geoplot([44.12 44.125],[-115.25 -115],'k-','linew',2)
% geoplot([44.22 44.28],[-115.22 -114.95],'k-','linew',2)
% group3=find(data.lat<44.15&data.depth_FSL<25&data.lon<-115.0&data.lon>-115.15&data.Herr<Hor_err&data.Verr<Vert_err);
% group4=find(data.lat<=44.28&data.lat>44.23&data.depth_FSL<max_depth&data.lon<-114.95&data.lon>-115.22&data.Herr<Hor_err&data.Verr<Vert_err);

% 
caxis([4 16]) %turn this on for depth axis 
colormap('jet')
title(['EQ count = ', num2str(EQ_count)])
% plot([-114 -114],[min(data.lat) max(data.lat)],'k-') %UTM zone divider
% grid
% axis([min(data.lon(group1)) max(data.lon(group1)) min(data.lat(group1)) max(data.lat(group1))])
h1=colorbar
ylabel(h1, 'Depth (km)')
geolimits([44.05 44.55],[-115.25 -114.80])
%%
% subplot(2,4,3:4)
% plot(EQ_UTM(groupA,1),EQ_UTM(groupA,2),'b.')
% hold on
% plot(EQ_UTM(group1,1),EQ_UTM(group1,2),'r.')
% plot(f1a,EQ_UTM(group1,2),'k-','linew',4)
% plot(EQ_UTM(group2,1),EQ_UTM(group2,2),'g.')
% plot(f2a,EQ_UTM(group2,2),'k-','linew',4)
% plot(EQ_UTM(group3,1),EQ_UTM(group3,2),'m.')
% plot(f3a,EQ_UTM(group3,2),'k-','linew',4)
% plot(EQ_UTM(group4,1),EQ_UTM(group4,2),'c.')
% plot(f4a,EQ_UTM(group4,2),'k-','linew',4)
% axis('equal')
% grid minor
% 
% % dip=25;  %actual dip is 90-dip
% % disp(['fault dip=' num2str(90-dip)])
% title(['trend1= ',num2str(theta1) '  trend2= ',num2str(theta2)])
%%
figure
ax3=subplot(2,2,1);
% scatter(rotcoord1(group1,1)/1000-min(rotcoord1(group1,1)/1000),data.depth_FSL(group1),1,data.depth_FSL(group1)) % color code by depth 
scatter(rotcoord1(group1,1)/1000-min(rotcoord1(group1,1)/1000),data.depth_FSL(group1),1,data.days(group1)) % color code by time 
hold on
% plot([mean(rotcoord1(group2,1)/1000)+dip/2*tan(deg2rad(dip))/2 mean(rotcoord1(group2,1)/1000)-dip/2*tan(deg2rad(dip))/2],[0 20],'linew',4)
plot(f1b-min(rotcoord1(group1,1)/1000),data.depth_FSL(group1),'r-','linew',2)
axis('equal')
grid
view(0,-90)
title(['Dip= ',num2str(dip1,'%10.1f\n'),' trend= ',num2str(theta1,'%10.1f\n'),' n= ',num2str(length(group1)) ])
xlabel('Distance (km)')
ylabel('Depth (km)')
ylim([0 16])
% colormap(ax3,'copper')
% colorbar
colormap('jet')
% caxis([4 16]) % turn on for color coded by depth 
h2=colorbar;
ylabel(h2, ' Days Since Mainshock ')

ax4=subplot(222);
% scatter(rotcoord2(group2,1)/1000-min(rotcoord2(group2,1)/1000),data.depth_FSL(group2),2,data.depth_FSL(group2))%
% for color coded by depth
scatter(rotcoord2(group2,1)/1000-min(rotcoord2(group2,1)/1000),data.depth_FSL(group2),2,data.days(group2)) % color coded by origin time
hold on
% plot([mean(rotcoord1(group2,1)/1000)+dip/2*tan(deg2rad(dip))/2 mean(rotcoord1(group2,1)/1000)-dip/2*tan(deg2rad(dip))/2],[0 20],'linew',4)
% plot(f2b-min(rotcoord2(group2,1)/1000),data.depth_FSL(group2),'m-','linew',4)
axis('equal')
grid
view(0,-90)
title(['Dip= ',num2str(dip2,'%10.1f\n'),' trend= ',num2str(theta2,'%10.1f\n'),' n= ',num2str(length(group2)) ])
xlabel('Distance (km)')
ylabel('Depth (km)')
ylim([0 16])
% caxis([4 16]) uncommonet for depth

% colormap(ax4,'copper')
% colorbar

ax5=subplot(223);
% scatter(rotcoord2(group3,1)/1000-min(rotcoord2(group3,1)/1000),data.depth_FSL(group3),2,data.depth_FSL(group3))%color coded by depth 
scatter(rotcoord2(group3,1)/1000-min(rotcoord2(group3,1)/1000),data.depth_FSL(group3),2,data.days(group3)) % color coded by time
hold on
% plot([mean(rotcoord1(group2,1)/1000)+dip/2*tan(deg2rad(dip))/2 mean(rotcoord1(group2,1)/1000)-dip/2*tan(deg2rad(dip))/2],[0 20],'linew',4)
% plot(f2b-min(rotcoord2(group2,1)/1000),data.depth_FSL(group2),'g-','linew',4)
axis('equal')
grid
view(0,-90)
ylim([0 16])
caxis([4 16])

title([' trend= ',num2str(theta3,'%10.1f\n'),' n= ',num2str(length(group3)) ])
xlabel('Distance (km)')
ylabel('Depth (km)')
% colormap(ax5,'copper')
% colorbar

ax6=subplot(224);
% scatter(rotcoord4(group4,1)/1000-min(rotcoord4(group4,1)/1000),data.depth_FSL(group4),2,data.depth_FSL(group4))% color coded by depth
scatter(rotcoord4(group4,1)/1000-min(rotcoord4(group4,1)/1000),data.depth_FSL(group4),2,data.days(group4)) %color coded by time 
hold on
% plot([mean(rotcoord1(group2,1)/1000)+dip/2*tan(deg2rad(dip))/2 mean(rotcoord1(group2,1)/1000)-dip/2*tan(deg2rad(dip))/2],[0 20],'linew',4)
% plot(f4b-min(rotcoord2(group4,1)/1000),data.depth_FSL(group4),'c-','linew',4)
axis('equal')
grid
view(0,-90)
ylim([0 16])
% caxis([4 16]) uncomment for depth

title([' trend= ',num2str(theta4,'%10.1f\n'),' n= ',num2str(length(group4)) ])
xlabel('Distance (km)')
ylabel('Depth (km)')
% colormap(ax6,'copper')
% colorbar

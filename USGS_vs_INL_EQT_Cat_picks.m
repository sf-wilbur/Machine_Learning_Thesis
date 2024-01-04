clear all 
%%
Chal = readtable('EQT_Challis_Gradient.csv');
Stanl = readtable('EQt_Stanley_Gradient.csv');
sul = readtable('EQT_Sulphur_sulvel.csv');
sul_inl = readtable("inl_sulphur.csv");
chal_inl = readtable("inl_challis.csv");
stan_inl =load('norhin-seisan.mat'); 
stan_inl = struct2table(stan_inl.c);

T = readtable("1970_present_EQ_westernUS_M2.5.xlsx"); %read in USGS catalog 

%% change time format in USGS catalog 
formatIn = 'yyyy-mm-ddTHH:MM:SS'; %tell matlab the format for time in the catalog T
t1 = T.time(:);
t1 = datenum(t1, formatIn); %Get the correct format to compare the times between arrays 
T.time = t1;
%% USGS Challis
clear  i  idx chal
T2 = height(T);
chal = array2table(zeros(length(T2),7));
chal.Properties.VariableNames = {'otime','lon', ...
           'lat','depth', 'horErr', 'rms', 'mag'};


%for loop for horizontal error
for i = 1:T2
    if i>42494 && i< 46312 && 45 >= T.latitude(i)... 
        && T.latitude(i) >= 44 && T.longitude(i) >= -114.8 && T.longitude(i) <= -113.4 
        
        chal.otime(i) = T.time(i);
        chal.lat(i) = T.latitude(i);
        chal.lon(i) = T.longitude(i);
        chal.depth(i) = T.depth(i);
        chal.horErr(i) = T.horizontalError(i);
        chal.rms(i) = T.rms(i); 
        %quality 
        
    end
    
end 
chal(~chal.otime,:) = [];
%% USGS Stanley 
%% STANLEY create table
clear  i  idx St
T2 = height(T);
St = array2table(zeros(length(T2),7));
St.Properties.VariableNames = {'otime','lon', ...
           'lat','depth', 'horErr', 'rms', 'mag'};

%for loop for horizontal error
for i = 1:T2
    if i>49630 && i< 52643 && 45 >= T.latitude(i)... 
        && T.latitude(i) >= 43.5 && T.longitude(i) >= -115.5 && T.longitude(i) <= -114 
       
        St.otime(i) = T.time(i);
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
%% USGS create Table sulphur




% S = Catalog.retrieve('iris', 'minimumLatitude', '41.8','maximumLatitude', '43.56', 'minimumLongitude', '-113.15', 'maximumLongitude', '-110.1','starttime', '2017-09-01', 'endtime', '2017-11-01')
clear  i  idx S; 
S = array2table(zeros(length(T2),6));
S.Properties.VariableNames = {'otime','lon', ...
           'lat','depth', 'horErr', 'rms'}

%for loop for horizontal error
for i = 1:T2
    if i>46694 && i< 47421 && 43.5 >= T.latitude(i)... 
        && T.latitude(i) >= 42 && T.longitude(i) >= -112 && T.longitude(i) <= -111 
        
        S.otime(i) = T.time(i);
        S.lat(i) = T.latitude(i);
        S.lon(i) = T.longitude(i);
        S.depth(i) = T.depth(i);
        S.horErr(i) = T.horizontalError(i); 
        S.rms(i) = T.rms(i); 
        %quality 
        
    end
    
end 
S(~S.otime,:) = [];

%% Set Up Time arrays for USGS AND stanley EQT COMPARISON
clear t1 t2_n t2 ii;
t1 = St.otime(:,1); % USGS origin times
t2_st = zeros(height(Stanl),1);
t2_st = Stanl.otime ;                  %t2 = inl(k).otime; for structure


t3 = stan_inl.otime(:,1); %inl orogin times 

% MATCH THE USGS EVENTS TO EQT for Stanley
clear dt dist val ids counter ; 
 %Change the format to math the sum files and compare the times
E = referenceEllipsoid('Earth'); % reference ellipse [m] for distance calculation
counter = 0;

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
        lat_diff(ii) = Stanl.lat(ids(ii))-St.lat(ii); 
        lon_diff(ii) = Stanl.lon(ids(ii))-St.lon(ii); 
                % compute distance in meters
        dist(ii) = distance( Stanl.lat(ids(ii)), Stanl.lon(ids(ii)),...  % right now this is set to run for A events 
        St.lat(ii), St.lon(ii), E) ; %CHANGE CATALOGS HERE 
        
    end 
   % Append mag values from USGS to Stanley cat 
    Stanl.mag(ids(ii))= St.mag(ii);
end
counter_inl = 0;

for ii = 1 : numel(t3) % loop through eqt events CHANGE TO T1
    
    % find the USGS event that is closest in time to eqt_event(i)
    t_diff_inl  = t2_st - t3(ii); 
    
    [val, ids(ii)] = min( abs(t_diff_inl) ); % this lets them be negative 

    % Display some useful information to the user
    
  
    

    fprintf('EQT: %s, USGS: %s\n', datestr(t2_st(ids(ii))), datestr(t3(ii)));
    
    dt_inl(ii) = val * 24*3600; % [s] convert time difference to seconds
%      dt_nm(i) = dt(i);
    if dt_inl(ii) > 10
        counter_inl = counter_inl +1 ;
          
 

    else 

                % compute distance in meters
        dist_inll(ii) = distance( Stanl.lat(ids(ii)), Stanl.lon(ids(ii)),...  % right now this is set to run for A events 
        stan_inl.lat(ii), stan_inl.lon(ii), E) ; %CHANGE CATALOGS HERE 
   
% Plots for Stanley Sawtooth Vel (EDIT) 
    end
end


clear h

h = figure; 
set(gcf,'DefaultTextFontSize',18); hold on
set(gcf, 'DefaultAxesFontSize',18); 
% time differences: use 0.1s bins from 0 to 2 seconds
subplot(1,2,1); 
title_text = sprintf('Stanley Sequence n=%d', length(t1)-counter );
histogram(dt, 0:0.1:10); grid on;
xlabel('Origin time difference [s] EQT vs USGS'); 
ylabel('No. events'); title(title_text);

subplot(1,2,2);

title_text_inl = sprintf('Stanley Sequence n=%d', length(t3)-counter_inl );
histogram(dt_inl, 0:0.1:10); grid on;
xlabel('Origin time difference [s] EQT vs INL'); 
ylabel('No. events'); title(title_text_inl);
%%
% epicenter distance difference: use 1km bins from 0 to 15 kilometers
H = figure;
set(gcf,'DefaultTextFontSize',14); hold on
set(gcf, 'DefaultAxesFontSize',18);

subplot(1,2,1); 
histogram(lat_diff); grid on;
m_stv = mean(lat_diff);
st_stv = std(lat_diff);
MAD = mad(lat_diff); 
mdist = mean(dist/1000);
str5 = sprintf('MAD = %.3f', MAD);
str = sprintf('Mean = %.3f', m_stv);
str2 = sprintf('Std = %.3f', st_stv);
str4 = 'North';
str6 = 'South' ;
text(-1.5,400,str)
text(-1.5,350, str2)
text(-1.5, 300, str5)
text(-1.5, 470, str6)
text(0.75, 470, str4)
xlabel('Latitude difference [Degrees]'); 
ylabel('No. events'); title(title_text);dim = [0.414186507936508,0.960567053676594,0.202604166666667,0.031213191990577];
str3 = sprintf('Mean Difference in Epicenters (km) = %.3f', mdist);
annotation('textbox',dim,'String',str3,'FontSize', 14,'FitBoxToText','on');
subplot(1,2,2);
m_stv = mean(lon_diff);
st_stv = std(lon_diff);
MAD = mad(lon_diff); 
str5 = sprintf('MAD = %.3f', MAD);
str = sprintf('Mean = %.3f', m_stv);
str2 = sprintf('Std = %.3f', st_stv);
str4 = 'East' ;
str6 = 'West' ;
histogram(lon_diff); grid on;
text(-0.8, 422, str4);
text(1.5, 422, str6);
text(0.3,350,str)
text(0.3,300, str2)
text(0.3, 250, str5)
xlabel('Longitude difference [Degrees]'); 
ylabel('No. events'); title(title_text);

%% set up time arrays for USGS and Challis Comparison
clear t1 t2_n t2;
t1 = chal.otime(:,1); % USGS origin times
t2_st = Chal.otime;

t3 = chal_inl.otime(:,1); %inl origin times 
% set up time arrays for USGS and Challis Comparison
counter = 0;
E = referenceEllipsoid('Earth');
for ii = 1:numel(t1)
         % find the USGS event that is closest in time to eqt_event(i)
        t_diff  = t2_st - t1(ii); 
        [val, ids(ii)] = min( abs(t_diff) ); % this lets them be negative 

        fprintf('EQT: %s, USGS: %s\n', datestr(t2_st(ids(ii))), datestr(t1(ii)));
    
         dt(ii) = val * 24*3600; % [s] convert time difference to seconds
        if dt(ii) > 10
            counter = counter +1 ;
          
        end
        % compute distance in meters
        dist(ii) = distance( Chal.lat(ids(ii)), Chal.lon(ids(ii)),...  % right now this is set to run for A events 
        chal.lat(ii), chal.lon(ii), E) ; %CHANGE CATALOGS HERE 
        dist_nm(ii) = dist(ii); %create new vectors that keep track of these values
end

counter_inl = 0;

for ii = 1 : numel(t3) % loop through eqt events CHANGE TO T1
    
    % find the USGS event that is closest in time to eqt_event(i)
    t_diff_inl  = t2_st - t3(ii); 
    
    [val, ids(ii)] = min( abs(t_diff_inl) ); % this lets them be negative 

    % Display some useful information to the user
    
  
    

    fprintf('EQT: %s, USGS: %s\n', datestr(t2_st(ids(ii))), datestr(t3(ii)));
    
    dt_inl(ii) = val * 24*3600; % [s] convert time difference to seconds
%      dt_nm(i) = dt(i);
    if dt_inl(ii) > 10
        counter_inl = counter_inl +1 ;
          
 

    else 

                % compute distance in meters
        dist_inll(ii) = distance( Chal.lat(ids(ii)), Chal.lon(ids(ii)),...  % right now this is set to run for A events 
        chal_inl.lat(ii), chal_inl.lon(ii), E) ; %CHANGE CATALOGS HERE 
% Plots for Stanley Sawtooth Vel (EDIT) 
    end
end

clear h
title_text = sprintf('Challis Sequence n=%d', length(t1)-counter );
h = figure; 
set(gcf,'DefaultTextFontSize',18); hold on
set(gcf, 'DefaultAxesFontSize',18); 
% time differences: use 0.1s bins from 0 to 2 seconds
subplot(1,2,1); 
histogram(dt, 0:0.1:10); grid on;
xlabel('Origin time difference [s] EQT vs USGS'); 
ylabel('No. events'); title(title_text);

subplot(1,2,2);
title_text_inl = sprintf('Challis Sequence n=%d', length(t3)-counter_inl );
histogram(dt_inl, 0:0.1:10); grid on;
xlabel('Origin time difference [s] EQT vs INL'); 
ylabel('No. events'); title(title_text_inl);

%% set up time arrays for USGS and sulphur Comparison
clear t1 t2_n t2 t2_st t_diff_inl t_diff counter counter_inl;
t1 = S.otime(:,1); % USGS origin times
t2_st = sul.otime;
tmin = min(sul.otime);
t1_usgs=find(S.otime>=tmin);
t1 = S.otime(t1_usgs);
% set up time arrays for USGS and Sulphur Comparison
counter = 0;
E = referenceEllipsoid('Earth');
for ii = 1:numel(t1)
         % find the USGS event that is closest in time to eqt_event(i)
        t_diff  = t2_st - t1(ii); 
        [val, ids(ii)] = min( abs(t_diff) ); % this lets them be negative 

        fprintf('EQT: %s, USGS: %s\n', datestr(t2_st(ids(ii))), datestr(t1(ii)));
    
         dt(ii) = val * 24*3600; % [s] convert time difference to seconds
        if dt(ii) > 10
            counter = counter +1 ;
          
        end
        % compute distance in meters
        dist(ii) = distance( sul.lat(ids(ii)), sul.lon(ids(ii)),...  % right now this is set to run for A events 
        S.lat(ii), S.lon(ii), E) ; %CHANGE CATALOGS HERE 
        dist_nm(ii) = dist(ii); %create new vectors that keep track of these values
end


counter_inl = 0;
t1_inl=find(sul_inl.otime>=tmin);
t3 = sul_inl.otime(t1_inl); 
for ii = 1 : numel(t3) % loop through eqt events CHANGE TO T1
    
    % find the INL event that is closest in time to eqt_event(i)
    t_diff_inl  = t2_st - t3(ii); 
    
    [val, ids(ii)] = min( abs(t_diff_inl) ); % this lets them be negative 

    % Display some useful information to the user
    
  
    

    fprintf('EQT: %s, USGS: %s\n', datestr(t2_st(ids(ii))), datestr(t3(ii)));
    
    dt_inl(ii) = val * 24*3600; % [s] convert time difference to seconds
%      dt_nm(i) = dt(i);
    if dt_inl(ii) > 10
        counter_inl = counter_inl +1 ;
          
 

    else 

                % compute distance in meters
        dist_inll(ii) = distance( sul.lat(ids(ii)), sul.lon(ids(ii)),...  % right now this is set to run for A events 
        sul_inl.lat(ii), sul_inl.lon(ii), E) ; %CHANGE CATALOGS HERE 
% Plots for Stanley Sawtooth Vel (EDIT) 
    end
end

clear h
title_text = sprintf('Sulphur Sequence n=%d', length(t1)-counter );
h = figure; 
set(gcf,'DefaultTextFontSize',18); hold on
set(gcf, 'DefaultAxesFontSize',18); 
% time differences: use 0.1s bins from 0 to 2 seconds
subplot(1,2,1); 
histogram(dt, 0:0.1:10); grid on;
xlabel('Origin time difference [s] EQT vs USGS'); 
ylabel('No. events'); title(title_text);

subplot(1,2,2);
title_text_inl = sprintf('Sulphur Sequence n=%d', length(t3)-counter_inl );
histogram(dt_inl, 0:0.1:10); grid on;
xlabel('Origin time difference [s] EQT vs INL'); 
ylabel('No. events'); title(title_text_inl);

%% Add magnitude to stanley catalog events csv 
%Read in stanley local velocity model 
writetable(Stanl, 'EQT_Stanley_Gradient.csv','Delimiter',',')
type 'EQT_Stanley_Gradient.csv'

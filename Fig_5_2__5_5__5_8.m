close all;
clear;
clc;
%read in stanley Data
load EQ_Catalogs/EQT_ST.mat %load Stanley timetable 
% weeks = EQT_ST.Dates; % these are the x-axis values.
weeks = min(EQT_ST.Dates):calweeks(1):max(EQT_ST.Dates);
%Read in challis Data
EQT_Challis = readtable('EQ_Catalogs/EQT_Challis_QE.csv'); %read in table
dates = datetime(EQT_Challis.otime, 'convertfrom', 'datenum'); % create a Datetime column 
EQT_Challis.Dates = dates;
EQT_Challis_TT = table2timetable(EQT_Challis);
t1c = min(EQT_Challis_TT.Dates);EQT_Challis_TT.quality = [];
t2c = max(EQT_Challis_TT.Dates);EQT_Challis_TT.magtype = [];
EQT_CT = retime(EQT_Challis_TT,'weekly','mean');
t1c = min(EQT_Challis_TT.Dates);
t2c = max(EQT_Challis_TT.Dates);
months = t1c:calweeks(1):t2c;
%read in Sulphure
EQT_Sulphur = readtable('EQ_Catalogs/EQT_Sulphur_QE.csv');
EQT_Sul_TT = table2timetable(EQT_Sulphur);
EQT_Sul_TT.quality = [];
EQT_Sul_TT.magtype = [];
EQT_Sul = retime(EQT_Sul_TT,'weekly','mean');
t1s = min(EQT_Sul_TT.d_time);
t2s = max(EQT_Sul_TT.d_time);
week_s = t1s:calweeks(1):t2s;

load Station_Data/Challis_Network.mat
load Station_Data/XP_Network.mat
load Station_Data/Sulphur_StatList.mat
% There were two empty rows that I removed.
XPNetworkList(end-1:end,:) = [];

station_count = zeros(1, numel(weeks)); % an empty vector the length of your x-axis (i.e. weeks)
station_countc = zeros(1, numel(months)); 
station_counts = zeros(1,numel(week_s)); 
% for loop  for stanley through each station
for ii = 1 :size(XPNetworkList,1)
   t1 = XPNetworkList.StartTime(ii);
   t2 = XPNetworkList.EndTime(ii);
   tmp_count = (weeks >= t1) & (weeks <= t2);
   station_count = station_count + tmp_count; % add this station to the total
end

% For loop for challis 
for ic = 1 :size(challisstationlist,1)
   t1c = challisstationlist.StartTime(ic);
   t2c = challisstationlist.EndTime(ic);
   tmp_countc = (months >= t1c) & (months <= t2c);
   station_countc = station_countc + tmp_countc; % add this station to the total
end

% For loop for sulphur 
for is = 1 :size(sulphurstationlist,1)
   t1s = sulphurstationlist.StartTime(is);
   t2s = sulphurstationlist.EndTime(is);
   tmp_counts = (week_s >= t1s) & (week_s <= t2s);
   station_counts = station_counts + tmp_counts; % add this station to the total
end
station_counts= [1; station_counts'];
% For minimal lines of code to do this!! The above loop is easier to read
% tough. You can delete the commented code below. Just leaving it in so you
% can see it.
% 
% station_count = zeros(1, numel(weeks)); % an empty vector the length of your x-axis (i.e. weeks)
% for ii = 1 :size(XPNetworkList,1)
%    station_count = station_count +...
%        (weeks >= XPNetworkList.StartTime(ii)) & (weeks <= XPNetworkList.EndTime(ii)); % add this station to the total
% end

h = figure('Color','w');

yyaxis left
% plot horizontal error
l1 = plot(EQT_ST.otime, EQT_ST.Herr, '-'); hold on;
% plot vertical error with same color, but different line style
plot(EQT_ST.otime, EQT_ST.Verr, '--', 'Color', l1.Color); grid on;
axis([EQT_ST.otime(1), EQT_ST.otime(end), 0, 3]);
ylabel('Error (km)'); 

yyaxis right
plot(EQT_ST.otime, station_count, '-o');
ylim([0, 20]);
ylabel('Station Count');

legend({'\sigma_H','\sigma_Z','Temporary Network'},...
    'Orientation','Horizontal','Location','NorthOutside')
xlabel('Date'); datetick('x',6); 

print(h,'error_vs_station.png','-dpng');

%figure for challis
c = figure('Color','w');
yyaxis left
% plot horizontal error
l2 = plot(EQT_CT.otime, EQT_CT.Herr, '-'); hold on;
% plot vertical error with same color, but different line style
plot(EQT_CT.otime, EQT_CT.Verr, '--', 'Color', l2.Color); grid on;
axis([EQT_CT.otime(1), EQT_CT.otime(end), 0, 3]);
ylabel('Error (km)'); 

yyaxis right
plot(EQT_CT.otime, station_countc, '-o');
ylim([0, 20]);
ylabel('Station Count');

legend({'\sigma_H','\sigma_Z','Temporary Network'},...
    'Orientation','Horizontal','Location','NorthOutside')
xlabel('Date'); datetick('x',6); 

print(c,'error_vs_station.png','-dpng');


%figure for sulphur
s = figure('Color','w');
yyaxis left
% plot horizontal error
l2 = plot(EQT_Sul.otime, EQT_Sul.Herr, '-'); hold on;
% plot vertical error with same color, but different line style
plot(EQT_Sul.otime, EQT_Sul.Verr, '--', 'Color', l2.Color); grid on;
axis([EQT_Sul.otime(1), EQT_Sul.otime(end), 0, 3]);
ylabel('Error (km)'); 

yyaxis right
plot(EQT_Sul.otime, station_counts, '-o');
ylim([0, 10]);
ylabel('Station Count');
title('Sulphur Peak Residual Error vs Station Deployment')
legend({'\sigma_H','\sigma_Z','Temporary Network'},...
    'Orientation','Horizontal','Location','NorthOutside')
xlabel('Date'); datetick('x',6); 

print(c,'error_vs_station.png','-dpng');



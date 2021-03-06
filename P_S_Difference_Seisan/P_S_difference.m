clear; close all; clc;

% load the two catalogs
load('PICK_DATA/eqt-all.mat');
load('norhin-seisan.mat');

% Rename the catalogs like this
eqt = S; % ...then you dont' have 'eqt.S.otime', but 'eqt.otime' instead. 
% It makes using these variables easier.
inl = c;
clear S c; % delete the old variables

%% Extract earthquake origin times from the two catalogs

% Match events
t_inl = [inl.otime]; % INL origin times
t_eqt = [eqt.otime]; % EQT origin times

% Plot the times just to check that things overlap
h = figure('color','w');
plot(t_eqt, t_eqt, 'k.'); hold on;
plot(t_inl, t_inl, 'ro');
legend({'EQT event','INL event'},'location','northwest');
datetick('x'); datetick('y');

set( findall( h, '-property', 'FontSize' ), 'FontSize', 18 );
set( findall( h, '-property', 'FontName' ), 'FontName', 'Helvetica' );
set( findall( h, '-property', 'FontWeight' ), 'FontWeight', 'Bold' );
set( h, 'Position', [100 100 900 900] );
set( h, 'PaperPositionMode', 'auto' );

print(h,'event_overlap','-dpng');

%% Let's also do some time history still plot (i.e. how many events per day

dayS = datenum('2020-03-31'); % start date
dayE = datenum('2021-01-01'); % end date
day_vector = dayS : 7 : dayE; % a vector of days for the histogram

XP_removed = datenum('2020-10-28');

% day histogram -- dh_*
h = figure('color','w');
subplot(1,2,1);
dh_eqt = histogram(t_eqt,day_vector); title('EQT'); hold on; ax = axis;
datetick('x',3); xtickangle(45); ylabel('No. Events'); 
grid on; axis tight;
plot([XP_removed, XP_removed],[ax(3) ax(4)],'r','linewidth',3); 

subplot(1,2,2);
dh_inl = histogram(t_inl,day_vector); title('Manual, M\geq2.5'); hold on; ax = axis;
datetick('x',3); xtickangle(45); ylabel('No. Events'); 
grid on; axis tight;
plot([XP_removed, XP_removed],[ax(3) ax(4)],'r','linewidth',3);

set( findall( h, '-property', 'FontSize' ), 'FontSize', 18 );
set( findall( h, '-property', 'FontName' ), 'FontName', 'Helvetica' );
set( findall( h, '-property', 'FontWeight' ), 'FontWeight', 'Bold' );
set( h, 'Position', [100 100 1000 400] );
set( h, 'PaperPositionMode', 'auto' );

print(h,'event_number_differences','-dpng');

%% Find the EQT event that is closest in time to each INL event.

for ii = 1 : numel( t_inl )
    
    time_sec = abs( t_inl(ii) - t_eqt ) * 24 * 3600; % time difference in seconds
    
    [~, eqt_ID(ii)] = min( time_sec ); % find the eqt_event index
    
end

% Compute the real time difference (i.e. can be positive or negative)
time_diff = ( t_inl - t_eqt(eqt_ID) ) .* 24 .* 3600;

h = figure('Color','w');
histogram(time_diff, linspace(-2, 2, 51) ); grid on; grid minor; axis('tight');
title('Manual minus EQT time')
xlabel('Origin time difference [s]'); ylabel('No. Events');
set( findall( h, '-property', 'FontSize' ), 'FontSize', 18 );
set( findall( h, '-property', 'FontName' ), 'FontName', 'Helvetica' );
set( findall( h, '-property', 'FontWeight' ), 'FontWeight', 'Bold' );
set( h, 'Position', [100 100 600 600] );
set( h, 'PaperPositionMode', 'auto' );

print(h,'t0_differences','-dpng');

%% Remove events when the time difference is greater than 10 seconds

kill_idx = (abs(time_diff) > 2);
inl(kill_idx) = []; % remove those events from the INL catalog

eqt = eqt(eqt_ID); % first pull out just the EQT events we want
eqt(kill_idx) = []; % remove those events from the EQT catalog

%% Now let's start matching the P waves

counter = 0; % this is going to counter all matching events

for ii = 1 : numel(inl)
    
    inl_c = inl( ii ); % get the INL info for this event
    eqt_c = eqt( ii ); % get the EQT info for this event
    
    fprintf('Processing event %d\n',ii);
    fprintf('INL time %s; EQT time %s\n',...
        datestr(inl_c.otime, 'yyyy-mm-dd HH:MM:SS.FFF'),...
        datestr(eqt_c.otime, 'yyyy-mm-dd HH:MM:SS.FFF') )
    dt = eqt_c.otime - inl_c.otime;
    fprintf('Origin time difference: %0.2f [s]\n', time_diff(ii) );
    
    inl_stats = {inl_c.P.stat}; % the INL p-wave stations
    eqt_stats = {eqt_c.P.stat}; % the EQT p-wave stations
    
    % find which stations are in both station lists
    [c1,ia,ib] = intersect(inl_stats, eqt_stats); 
    
    % loop over the stations that match
    for jj = 1 : numel(c1)
        counter = counter + 1; % add a P-wave pick
        p_diff(counter) = ( inl_c.P(ia(jj)).time - eqt_c.P(ib(jj)).time ) * 24*3600; % make sure to convert to seconds
    end
    
end

%%

hist_x = linspace(-0.5, 0.5, 101);
kill_pdx = (abs(p_diff)>2); %remove any p differences that were greater than 2 seconds in difference 
p_diff(kill_pdx) =[]; 
mean_p = mean(p_diff);
std_p = std(p_diff);
mae_p = mae(p_diff);

% Make a histogram of the P-pick differences
h = figure('Color','w');
h1 = histogram(p_diff, hist_x ); grid on; grid minor; axis('tight');
title('Manual minus EQT time');
text(hist_x(1)-hist_x(1)*0.1, max(h1.Values)-max(h1.Values)*0.1,sprintf('\\mu = %0.2f [s]',mean_p));
text(hist_x(1)-hist_x(1)*0.1, max(h1.Values)-max(h1.Values)*0.15,sprintf('\\sigma = %0.2f [s]',std_p));
text(hist_x(1)-hist_x(1)*0.1, max(h1.Values)-max(h1.Values)*0.2,sprintf('MAE = %0.2f [s]',mae_p));
xlabel('P-pick time difference [s]'); ylabel('No. Picks');
set( findall( h, '-property', 'FontSize' ), 'FontSize', 18 );
set( findall( h, '-property', 'FontName' ), 'FontName', 'Helvetica' );
set( findall( h, '-property', 'FontWeight' ), 'FontWeight', 'Bold' );
set( h, 'Position', [100 100 600 600] );
set( h, 'PaperPositionMode', 'auto' );

print(h,'p_differences','-dpng');


%% Now let's start matching the S waves

counter = 0; % this is going to counter all matching events

for ii = 1 : numel(inl)
    
    inl_c = inl( ii ); % get the INL info for this event
    eqt_c = eqt( ii ); % get the EQT info for this event
    
    fprintf('Processing event %d\n',ii);
    fprintf('INL time %s; EQT time %s\n',...
        datestr(inl_c.otime, 'yyyy-mm-dd HH:MM:SS.FFF'),...
        datestr(eqt_c.otime, 'yyyy-mm-dd HH:MM:SS.FFF') )
    dt = eqt_c.otime - inl_c.otime;
    fprintf('Origin time difference: %0.2f [s]\n', time_diff(ii) );
    
    inl_stats = {inl_c.S.stat}; % the INL p-wave stations
    eqt_stats = {eqt_c.S.stat}; % the EQT p-wave stations
    
    % find which stations are in both station lists
    [c1,ia,ib] = intersect(inl_stats, eqt_stats); 
    
    % loop over the stations that match
    for jj = 1 : numel(c1)
        counter = counter + 1; % add a P-wave pick
        s_diff(counter) = ( inl_c.S(ia(jj)).time - eqt_c.S(ib(jj)).time ) * 24*3600; % make sure to convert to seconds
    end
    
end

%%

mean_s = mean(s_diff);
std_s = std(s_diff);
mae_s = mae(s_diff);
kill_sdx = (abs(s_diff)>2);
s_diff(kill_sdx)=[]; %remove any s pick differences greater than 2 seconds 
% Make a histogram of the P-pick differences
h = figure('Color','w');
h1 = histogram(s_diff, hist_x ); grid on; grid minor; axis('tight');
title('Manual minus EQT time');
text(hist_x(1)-hist_x(1)*0.1, max(h1.Values)-max(h1.Values)*0.1,sprintf('\\mu = %0.4f [s]',mean_s));
text(hist_x(1)-hist_x(1)*0.1, max(h1.Values)-max(h1.Values)*0.15,sprintf('\\sigma = %0.2f [s]',std_s));
text(hist_x(1)-hist_x(1)*0.1, max(h1.Values)-max(h1.Values)*0.2,sprintf('MAE = %0.2f [s]',mae_s));
xlabel('S-pick time difference [s]'); ylabel('No. Picks');
set( findall( h, '-property', 'FontSize' ), 'FontSize', 18 );
set( findall( h, '-property', 'FontName' ), 'FontName', 'Helvetica' );
set( findall( h, '-property', 'FontWeight' ), 'FontWeight', 'Bold' );
set( h, 'Position', [100 100 600 600] );
set( h, 'PaperPositionMode', 'auto' );

print(h,'s_differences','-dpng');
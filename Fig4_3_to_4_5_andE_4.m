clear all; 
% Read in catalogs 
EQT_Challis=readtable('EQ_Catalogs/EQT_Challis_Gradient.csv'); % gradient model for challis 
EQT_Sulphur=readtable('EQ_Catalogs/EQT_Sulphur_sulvel.csv'); % sulphur with local velocity model 
EQT_Stanley=readtable('EQ_Catalogs/EQT_Stanley_Gradient.csv'); % gradient velocity model
% For stanley 
% EQT_Stanley_AK.Herr<=2.5&EQT_Stanley_AK.Verr<=4.2)
A_stan_idx=find(strcmp('A', EQT_Stanley.quality));
B_stan_idx=find(strcmp('B',EQT_Stanley.quality));
C_stan_idx=find(strcmp('C',EQT_Stanley.quality));
D_stan_idx=find(strcmp('D',EQT_Stanley.quality));
% A_B_stan=[EQT_Stanley.quality(A_stan_idx) & EQT_Stanley.quality(B_stan_idx)]; % I CHANGED THIS TO REFLECT ONLY THE PARAMETERS
T1_stan_idx =find(EQT_Stanley.Herr<=5&EQT_Stanley.Verr<=5);
% For Challis 
EQT_Challis.lon=EQT_Challis.lon*-1;
A_chal_idx=find(strcmp('A',EQT_Challis.quality));
B_chal_idx=find(strcmp('B',EQT_Challis.quality));
C_chal_idx=find(strcmp('C',EQT_Challis.quality));
D_chal_idx=find(strcmp('D',EQT_Challis.quality));
T1_chal_idx=find(EQT_Challis.Herr<=5&EQT_Challis.Verr<=5);
% T2_chal_idx=find(EQT_Challis.Herr<1&EQT_Challis.Verr<1&EQT_Challis.Herr>=1&EQT_Challis.Verr>=2);

% For sulphur
A_sul_idx=find(strcmp('A',EQT_Sulphur.quality));
B_sul_idx=find(strcmp('B',EQT_Sulphur.quality));
C_sul_idx=find(strcmp('C',EQT_Sulphur.quality));
D_sul_idx = find(strcmp('D', EQT_Sulphur.quality)); 
T1_sul_idx=find(EQT_Sulphur.Herr<=5&EQT_Sulphur.Verr<=5);

%% STanley Cross Correlations 

figure(1)
subplot(2, 3, 1)
plot(EQT_Stanley.rmse(D_stan_idx(:)),EQT_Stanley.Herr(D_stan_idx(:)), 'mo');hold on;  %D stans 
plot(EQT_Stanley.rmse(C_stan_idx(:)),EQT_Stanley.Herr(C_stan_idx(:)), 'go'); %c stans
plot(EQT_Stanley.rmse(B_stan_idx(:)),EQT_Stanley.Herr(B_stan_idx(:)), 'bo'); %b stans 
plot(EQT_Stanley.rmse(A_stan_idx(:)),EQT_Stanley.Herr(A_stan_idx(:)), 'ro'); %A stans 

xlim([0 5]); 
ylim([0 55]);
xlabel('RMSE (sec)')
ylabel('Horizontal Error (km)')
title('Stanley EQT  Gradient Model ')
legend('D events','C events', 'B events', 'A events')
hold off; 


subplot(2,3,2)
plot(EQT_Stanley.rmse(D_stan_idx(:)),EQT_Stanley.Verr(D_stan_idx(:)), 'mo'); hold on; %D stans 
plot(EQT_Stanley.rmse(C_stan_idx(:)),EQT_Stanley.Verr(C_stan_idx(:)), 'go'); %c stans
plot(EQT_Stanley.rmse(B_stan_idx(:)),EQT_Stanley.Verr(B_stan_idx(:)), 'bo'); %b stans 
plot(EQT_Stanley.rmse(A_stan_idx(:)),EQT_Stanley.Verr(A_stan_idx(:)), 'ro'); %A events 

xlim([0 10]); 
ylim([0 105]);
xlabel('RMSE (sec)')
ylabel('Vertical Error (km)')
title('Stanley EQT Gradient Model ')
legend('D events','C events', 'B events', 'A events')
hold off; 

subplot(2,3,3)
plot(EQT_Stanley.Herr(D_stan_idx(:)),EQT_Stanley.Verr(D_stan_idx(:)), 'mo'); hold on; %D stans 
plot(EQT_Stanley.Herr(C_stan_idx(:)),EQT_Stanley.Verr(C_stan_idx(:)), 'go'); %c stans
plot(EQT_Stanley.Herr(B_stan_idx(:)),EQT_Stanley.Verr(B_stan_idx(:)), 'bo'); %b stans 
plot(EQT_Stanley.Herr(A_stan_idx(:)),EQT_Stanley.Verr(A_stan_idx(:)), 'ro'); %A events 
xlim([0 10]); 
ylim([0 105]);
legend('D events','C events', 'B events', 'A events')
xlabel('Horizontal Error (km)')
ylabel('Vertical Error (km)')
title('Stanley EQT Gradient Model ')


subplot(2, 3, 4)
plot(EQT_Stanley.rmse(D_stan_idx(:)),EQT_Stanley.Herr(D_stan_idx(:)), 'mo');hold on;  %D stans 
plot(EQT_Stanley.rmse(C_stan_idx(:)),EQT_Stanley.Herr(C_stan_idx(:)), 'go'); %c stans
plot(EQT_Stanley.rmse(B_stan_idx(:)),EQT_Stanley.Herr(B_stan_idx(:)), 'bo'); %b stans 
plot(EQT_Stanley.rmse(A_stan_idx(:)),EQT_Stanley.Herr(A_stan_idx(:)), 'ro'); %A stans 

xlim([0 3]); 
ylim([0 5]);
xlabel('RMSE (sec)')
ylabel('Horizontal Error (km)')
title('Stanley EQT  Gradient Model ')
legend('D events','C events', 'B events', 'A events')
hold off; 


subplot(2,3,5)
plot(EQT_Stanley.rmse(D_stan_idx(:)),EQT_Stanley.Verr(D_stan_idx(:)), 'mo'); hold on; %D stans 
plot(EQT_Stanley.rmse(C_stan_idx(:)),EQT_Stanley.Verr(C_stan_idx(:)), 'go'); %c stans
plot(EQT_Stanley.rmse(B_stan_idx(:)),EQT_Stanley.Verr(B_stan_idx(:)), 'bo'); %b stans 
plot(EQT_Stanley.rmse(A_stan_idx(:)),EQT_Stanley.Verr(A_stan_idx(:)), 'ro'); %A events 

xlim([0 3]); 
ylim([0 5]);
xlabel('RMSE (sec)')
ylabel('Vertical Error (km)')
title('Stanley EQT Gradient Model ')
legend('D events','C events', 'B events', 'A events')
hold off; 

subplot(2,3,6)
plot(EQT_Stanley.Herr(D_stan_idx(:)),EQT_Stanley.Verr(D_stan_idx(:)), 'mo'); hold on; %D stans 
plot(EQT_Stanley.Herr(C_stan_idx(:)),EQT_Stanley.Verr(C_stan_idx(:)), 'go'); %c stans
plot(EQT_Stanley.Herr(B_stan_idx(:)),EQT_Stanley.Verr(B_stan_idx(:)), 'bo'); %b stans 
plot(EQT_Stanley.Herr(A_stan_idx(:)),EQT_Stanley.Verr(A_stan_idx(:)), 'ro'); %A events 
xlim([0 5]); 
ylim([0 5]);
legend('D events','C events', 'B events', 'A events')
xlabel('Horizontal Error (km)')
ylabel('Vertical Error (km)')
title('Stanley EQT Gradient Model ')


%% For Challis 
figure(2)
subplot(2, 3, 1)
plot(EQT_Challis.rmse(D_chal_idx(:)),EQT_Challis.Herr(D_chal_idx(:)), 'mo');hold on;  %D chals 
plot(EQT_Challis.rmse(C_chal_idx(:)),EQT_Challis.Herr(C_chal_idx(:)), 'go'); %c chals
plot(EQT_Challis.rmse(B_chal_idx(:)),EQT_Challis.Herr(B_chal_idx(:)), 'bo'); %b chals 
plot(EQT_Challis.rmse(A_chal_idx(:)),EQT_Challis.Herr(A_chal_idx(:)), 'ro'); %A chals 

xlim([0 5]); 
ylim([0 55]);
xlabel('RMSE (sec)')
ylabel('Horizontal Error (km)')
title('Challis EQT  Gradient Model ')
legend('D events','C events', 'B events', 'A events')
hold off; 


subplot(2,3,2)
plot(EQT_Challis.rmse(D_chal_idx(:)),EQT_Challis.Verr(D_chal_idx(:)), 'mo'); hold on; %D events 
plot(EQT_Challis.rmse(C_chal_idx(:)),EQT_Challis.Verr(C_chal_idx(:)), 'go'); %c chals
plot(EQT_Challis.rmse(B_chal_idx(:)),EQT_Challis.Verr(B_chal_idx(:)), 'bo'); %b chals 
plot(EQT_Challis.rmse(A_chal_idx(:)),EQT_Challis.Verr(A_chal_idx(:)), 'ro'); %A chals 

xlim([0 10]); 
ylim([0 105]);
xlabel('RMSE (sec)')
ylabel('Vertical Error (km)')
title('Challis EQT Gradient Model ')
legend('D events','C events', 'B events', 'A events')
hold off; 

subplot(2,3,3)

plot(EQT_Challis.Herr(D_chal_idx(:)),EQT_Challis.Verr(D_chal_idx(:)), 'mo'); hold on; %D stans 
plot(EQT_Challis.Herr(C_chal_idx(:)),EQT_Challis.Verr(C_chal_idx(:)), 'go'); %c stans
plot(EQT_Challis.Herr(B_chal_idx(:)),EQT_Challis.Verr(B_chal_idx(:)), 'bo'); %b stans 
plot(EQT_Challis.Herr(A_chal_idx(:)),EQT_Challis.Verr(A_chal_idx(:)), 'ro'); %A events 
xlim([0 10]); 
ylim([0 105]);
legend('D events','C events', 'B events', 'A events')
xlabel('Horizontal Error (km)')
ylabel('Vertical Error (km)')
title('Challis EQT Gradient Model ')


subplot(2, 3, 4)
plot(EQT_Challis.rmse(D_chal_idx(:)),EQT_Challis.Herr(D_chal_idx(:)), 'mo');hold on;  %D chals 
plot(EQT_Challis.rmse(C_chal_idx(:)),EQT_Challis.Herr(C_chal_idx(:)), 'go'); %c chals
plot(EQT_Challis.rmse(B_chal_idx(:)),EQT_Challis.Herr(B_chal_idx(:)), 'bo'); %b chals 
plot(EQT_Challis.rmse(A_chal_idx(:)),EQT_Challis.Herr(A_chal_idx(:)), 'ro'); %A chals 

xlim([0 1.5]); 
ylim([0 5]);
xlabel('RMSE (sec)')
ylabel('Horizontal Error (km)')
title('Challis EQT  Gradient Model ')
legend('D events','C events', 'B events', 'A events')
hold off; 


subplot(2,3,5)
plot(EQT_Challis.rmse(D_chal_idx(:)),EQT_Challis.Verr(D_chal_idx(:)), 'mo'); hold on; %D events 
plot(EQT_Challis.rmse(C_chal_idx(:)),EQT_Challis.Verr(C_chal_idx(:)), 'go'); %c chals
plot(EQT_Challis.rmse(B_chal_idx(:)),EQT_Challis.Verr(B_chal_idx(:)), 'bo'); %b chals 
plot(EQT_Challis.rmse(A_chal_idx(:)),EQT_Challis.Verr(A_chal_idx(:)), 'ro'); %A chals 

xlim([0 1.5]); 
ylim([0 5]);
xlabel('RMSE (sec)')
ylabel('Vertical Error (km)')
title('Challis EQT Gradient Model ')
legend('D events','C events', 'B events', 'A events')
hold off; 


subplot(2,3,6)
plot(EQT_Challis.Herr(D_chal_idx(:)),EQT_Challis.Verr(D_chal_idx(:)), 'mo'); hold on; %D stans 
plot(EQT_Challis.Herr(C_chal_idx(:)),EQT_Challis.Verr(C_chal_idx(:)), 'go'); %c stans
plot(EQT_Challis.Herr(B_chal_idx(:)),EQT_Challis.Verr(B_chal_idx(:)), 'bo'); %b stans 
plot(EQT_Challis.Herr(A_chal_idx(:)),EQT_Challis.Verr(A_chal_idx(:)), 'ro'); %A events 
xlim([0 5]); 
ylim([0 5]);
legend('D events','C events', 'B events', 'A events')
xlabel('Horizontal Error (km)')
ylabel('Vertical Error (km)')
title('Challis EQT Gradient Model ')


%% For sulphur
figure(3)
subplot(2, 3, 1)
plot(EQT_Sulphur.rmse(D_sul_idx(:)),EQT_Sulphur.Herr(D_sul_idx(:)), 'mo');hold on;  %D suls 
plot(EQT_Sulphur.rmse(C_sul_idx(:)),EQT_Sulphur.Herr(C_sul_idx(:)), 'go'); %c suls
plot(EQT_Sulphur.rmse(B_sul_idx(:)),EQT_Sulphur.Herr(B_sul_idx(:)), 'bo'); %b suls 
plot(EQT_Sulphur.rmse(A_sul_idx(:)),EQT_Sulphur.Herr(A_sul_idx(:)), 'ro'); %A suls 

xlim([0 5]); 
ylim([0 55]);
xlabel('RMSE (sec)')
ylabel('Horizontal Error (km)')
title('Sulphur EQT  Gradient Model ')
legend('D suls','C suls', 'B suls', 'A suls')
hold off; 


subplot(2,3,2)
plot(EQT_Sulphur.rmse(D_sul_idx(:)),EQT_Sulphur.Verr(D_sul_idx(:)), 'mo'); hold on; %D suls 
plot(EQT_Sulphur.rmse(C_sul_idx(:)),EQT_Sulphur.Verr(C_sul_idx(:)), 'go'); %c suls
plot(EQT_Sulphur.rmse(B_sul_idx(:)),EQT_Sulphur.Verr(B_sul_idx(:)), 'bo'); %b suls 
plot(EQT_Sulphur.rmse(A_sul_idx(:)),EQT_Sulphur.Verr(A_sul_idx(:)), 'ro'); %A suls 

xlim([0 10]); 
ylim([0 105]);
xlabel('RMSE (sec)')
ylabel('Vertical Error (km)')
title('Sulphur EQT Gradient Model ')
% legend('D suls','C events', 'B events', 'A events')
hold off; 

subplot(2,3,3)
plot(EQT_Sulphur.Herr(D_sul_idx(:)),EQT_Sulphur.Verr(D_sul_idx(:)), 'mo'); hold on; %D stans 
plot(EQT_Sulphur.Herr(C_sul_idx(:)),EQT_Sulphur.Verr(C_sul_idx(:)), 'go'); %c stans
plot(EQT_Sulphur.Herr(B_sul_idx(:)),EQT_Sulphur.Verr(B_sul_idx(:)), 'bo'); %b stans 
plot(EQT_Sulphur.Herr(A_sul_idx(:)),EQT_Sulphur.Verr(A_sul_idx(:)), 'ro'); %A events 
xlim([0 10]); 
ylim([0 105]);
legend('D events','C events', 'B events', 'A events')
xlabel('Horizontal Error (km)')
ylabel('Vertical Error (km)')
title('Sulphur EQT Gradient Model ')


subplot(2, 3, 4)
plot(EQT_Sulphur.rmse(D_sul_idx(:)),EQT_Sulphur.Herr(D_sul_idx(:)), 'mo');hold on;  %D chals 
plot(EQT_Sulphur.rmse(C_sul_idx(:)),EQT_Sulphur.Herr(C_sul_idx(:)), 'go'); %c chals
plot(EQT_Sulphur.rmse(B_sul_idx(:)),EQT_Sulphur.Herr(B_sul_idx(:)), 'bo'); %b chals 
plot(EQT_Sulphur.rmse(A_sul_idx(:)),EQT_Sulphur.Herr(A_sul_idx(:)), 'ro'); %A chals 

xlim([0 1]); 
ylim([0 5]);
xlabel('RMSE (sec)')
ylabel('Horizontal Error (km)')
title('Sulphur EQT  Gradient Model ')
legend('D events','C events', 'B events', 'A events')
hold off; 


subplot(2,3,5)
plot(EQT_Sulphur.rmse(D_sul_idx(:)),EQT_Sulphur.Verr(D_sul_idx(:)), 'mo'); hold on; %D events 
plot(EQT_Sulphur.rmse(C_sul_idx(:)),EQT_Sulphur.Verr(C_sul_idx(:)), 'go'); %c chals
plot(EQT_Sulphur.rmse(B_sul_idx(:)),EQT_Sulphur.Verr(B_sul_idx(:)), 'bo'); %b chals 
plot(EQT_Sulphur.rmse(A_sul_idx(:)),EQT_Sulphur.Verr(A_sul_idx(:)), 'ro'); %A chals 

xlim([0 1]); 
ylim([0 5]);
xlabel('RMSE (sec)')
ylabel('Vertical Error (km)')
title('Sulphur EQT Gradient Model ')
legend('D events','C events', 'B events', 'A events')
hold off; 


subplot(2,3,6)
plot(EQT_Sulphur.Herr(D_sul_idx(:)),EQT_Sulphur.Verr(D_sul_idx(:)), 'mo'); hold on; %D stans 
plot(EQT_Sulphur.Herr(C_sul_idx(:)),EQT_Sulphur.Verr(C_sul_idx(:)), 'go'); %c stans
plot(EQT_Sulphur.Herr(B_sul_idx(:)),EQT_Sulphur.Verr(B_sul_idx(:)), 'bo'); %b stans 
plot(EQT_Sulphur.Herr(A_sul_idx(:)),EQT_Sulphur.Verr(A_sul_idx(:)), 'ro'); %A events 
xlim([0 5]); 
ylim([0 5]);
legend('D events','C events', 'B events', 'A events')
xlabel('Horizontal Error (km)')
ylabel('Vertical Error (km)')
title('Sulphur EQT Gradient Model ')


%% Create Histograms 
figure(4); hold on
set(gcf,'DefaultTextFontSize',12); hold on
set(gcf, 'DefaultAxesFontSize',12); 
subplot(3,2,1)
N0 = EQT_Stanley.Verr<=5 & EQT_Stanley.Herr <=5;
N= EQT_Stanley.Verr(N0);
X_min=0;
X_max=5;
B=51;                           % # of bins between X_min and X_max
dB = 0.1 ;% bin size
BE=(X_min-0.05):dB:(X_max-0.05);    % bin edges; first and last edges are at X_min-dB and X_max+dB
H=histcounts(N,BE); % histogram of X
                    
Bc=(BE(2:end)+BE(1:(end-1)))/2; % bin centroids
% Cumulative mass function
CMF=cumsum(H);
yyaxis left
bar(Bc, H);
ylabel('Number of Events')
yyaxis right
hax=gca; 
hax.XGrid = 'on'; 
ylim([0,1])
hax.YTick = 0:.2:1;
yl = arrayfun(@(y)yline(hax, y,'LineStyle',hax.GridLineStyle,'Color',hax.GridColor,'Alpha', hax.GridAlpha), hax.YTick);
% set(hax,'Xgrid','on')
% set(hax,'Ygrid','on')
set(hax, 'XLimSpec', 'Tight')
ylim([0 1])
h = line(Bc, CMF/length(N));
h.LineWidth = 2;
xlabel('Vertical Error km')

m_stv =mean(N);
st_stv = std(N);
m = mode(N);
MAD = mad(N); 
ND = length(N);
str = sprintf('Mean = %.1f', m_stv);
str2 = sprintf('Std = %.1f', st_stv);
str3 = sprintf('N = %d', ND);
str4 = {'Stanley'};
str5 = sprintf('MAD = %.1f', MAD);
str6 = sprintf('Mode = %.1f', m);
text(4,0.5,str)
text(4,0.4, str2)
text(4,0.3, str3)
text(0.01, .9, str4)
text(4, 0.2, str5)
text(4, 0.1, str6)
subplot(3,2,2)

N = EQT_Stanley.Herr(N0); 
H=histcounts(N,BE);             % histogram of X
CMF=cumsum(H);

yyaxis left
bar(Bc, H);
yyaxis right
hax=gca; 
hax.XGrid = 'on'; 
ylim([0,1])
hax.YTick = 0:.2:1;
yl = arrayfun(@(y)yline(hax, y,'LineStyle',hax.GridLineStyle,'Color',hax.GridColor,'Alpha', hax.GridAlpha), hax.YTick);
 
% grid on 
% set(hax,'Xgrid','on')
% set(hax,'Ygrid','on')
set(hax, 'XLimSpec', 'Tight')
h = line(Bc, CMF/length(N));
h.LineWidth = 2;
ylim([0 1])
%  xlim([0 5])

xlabel('Horizontal Error km')
ylabel('Event Percentage')



m_stv =mean(N);
st_stv = std(N);
m = mode(N);
MAD = mad(N); 
ND = length(N);
str = sprintf('Mean = %.1f', m_stv);
str2 = sprintf('Std = %.1f', st_stv);
str3 = sprintf('N = %d', ND);
st4 = {'Stanley'};
str5 = sprintf('MAD = %.1f', MAD);
str6 = sprintf('Mode = %.1f', m);
text(4,0.5,str)
text(4,0.4, str2)
text(4,0.3, str3)
text(-0.05, .9, str4)
text(4, 0.2, str5)
text(4, 0.1, str6)
%Challis 
subplot(3,2,3)
N0 = EQT_Challis.Verr<=5 & EQT_Challis.Herr <=5;
N = EQT_Challis.Verr(N0);
% X_max=max(N);
% B=49;                           % # of bins between X_min and X_max
% dB=(X_max-X_min)/B;             % bin size
% BE=(X_min-dB):dB:(X_max+dB);    % bin edges; first and last edges are at X_min-dB and X_max+dB
H=histcounts(N,BE);             % histogram of X
                    % mass density function
% Bc=(BE(2:end)+BE(1:(end-1)))/2; % bin centroids
% Cumulative mass function
CMF=cumsum(H);
yyaxis left
bar(Bc, H);
ylabel('Number of Events')
yyaxis right
hax=gca; 
hax.XGrid = 'on'; 
ylim([0,1])
hax.YTick = 0:.2:1;
yl = arrayfun(@(y)yline(hax, y,'LineStyle',hax.GridLineStyle,'Color',hax.GridColor,'Alpha', hax.GridAlpha), hax.YTick);

set(hax, 'XLimSpec', 'Tight')
ylim([0 1])
h = line(Bc, CMF/length(N));
h.LineWidth = 2;
% xlim([0 5])
 xlabel('Vertical Error km')
 

m_stv = mean(N);
st_stv = std(N);
m = mode(N);
MAD = mad(N); 
ND = length(N);
str5 = sprintf('MAD = %.1f', MAD);
str6 = sprintf('Mode = %.1f', m);
str = sprintf('Mean = %.1f', m_stv);
str2 = sprintf('Std = %.1f', st_stv);
str3 = sprintf('N = %d', ND);
str4 = {'Challis'};
text(4,0.5,str)
text(4,0.4, str2)
text(4,0.3, str3)
text(0.1, .9, str4)
text(4, 0.2, str5)
text(4, 0.1, str6)
subplot(3,2,4)
N = EQT_Challis.Herr(N0);
 
% X_min=min(N);
% X_max=max(N);
% B=49;                           % # of bins between X_min and X_max
% dB=(X_max-X_min)/B;             % bin size
% BE=(X_min-dB):dB:(X_max+dB);    % bin edges; first and last edges are at X_min-dB and X_max+dB
H=histcounts(N,BE);             % histogram of X
%                     % mass density function
% Bc=(BE(2:end)+BE(1:(end-1)))/2; % bin centroids
% Cumulative mass function
CMF=cumsum(H);

yyaxis left
bar(Bc, H);
yyaxis right
hax=gca; 
hax.XGrid = 'on'; 
ylim([0,1])
hax.YTick = 0:.2:1;
yl = arrayfun(@(y)yline(hax, y,'LineStyle',hax.GridLineStyle,'Color',hax.GridColor,'Alpha', hax.GridAlpha), hax.YTick);

set(hax, 'XLimSpec', 'Tight')
ylim([0 1])
h = line(Bc, CMF/length(N));
h.LineWidth = 2;



xlabel('Horizontal Error km')
ylabel('Event Percentage')


m_stv =mean(N);
st_stv = std(N);
m = mode(N);
MAD = mad(N); 
ND = length(N);
str = sprintf('Mean = %.1f', m_stv);
str2 = sprintf('Std = %.1f', st_stv);
str3 = sprintf('N = %d', ND);
str4 ={'Challis'};
str5 = sprintf('MAD = %.1f', MAD);
str6 = sprintf('Mode = %.1f', m);
text(4,0.5,str)
text(4,0.4, str2)
text(4,0.3, str3)
text(0.1, .9, str4)
text(4, 0.2, str5)
text(4, 0.1, str6)
% Sulphur 
subplot(3,2,5)
N0 = EQT_Sulphur.Verr<=5 & EQT_Sulphur.Herr <=5;
N = EQT_Sulphur.Verr(N0);

% X_min=min(N);
% X_max=max(N);
% B=49;                           % # of bins between X_min and X_max
% dB=(X_max-X_min)/B;             % bin size
% BE=(X_min-dB):dB:(X_max+dB);    % bin edges; first and last edges are at X_min-dB and X_max+dB
H=histcounts(N,BE);             % histogram of X
                    % mass density function
% Bc=(BE(2:end)+BE(1:(end-1)))/2; % bin centroids
% Cumulative mass function
CMF=cumsum(H);
yyaxis left
bar(Bc, H);
ylabel('Number of Events')
yyaxis right
hax=gca; 
hax.XGrid = 'on'; 
ylim([0,1])
hax.YTick = 0:.2:1;
yl = arrayfun(@(y)yline(hax, y,'LineStyle',hax.GridLineStyle,'Color',hax.GridColor,'Alpha', hax.GridAlpha), hax.YTick);

set(hax, 'XLimSpec', 'Tight')
ylim([0 1])
h = line(Bc, CMF/length(N));
h.LineWidth = 2;


 xlabel('Vertical Error km')


m_stv =mean(N);
st_stv = std(N);
m = mode(N);
MAD = mad(N); 
ND = length(N);
str = sprintf('Mean = %.1f', m_stv);
str2 = sprintf('Std = %.1f', st_stv);
str3 = sprintf('N = %d', ND);
str4 = {'Sulphur'};
str5 = sprintf('MAD = %.1f', MAD);
str6 = sprintf('Mode = %.1f', m);
text(4,0.5,str)
text(4,0.4, str2)
text(4,0.3, str3)
text(0.01, .9, str4)
text(4, 0.2, str5)
text(4, 0.1, str6)
subplot(3,2,6)
N = EQT_Sulphur.Herr(N0);

% X_min=min(N);
% X_max=max(N);
% B=49;                           % # of bins between X_min and X_max
% dB=(X_max-X_min)/B;             % bin size
% BE=(X_min-dB):dB:(X_max+dB);    % bin edges; first and last edges are at X_min-dB and X_max+dB
H=histcounts(N,BE);             % histogram of X
                    % mass density function
% Bc = linspace(0.1,4.9, 50);
% Bc=(BE(2:end)+BE(1:(end-1)))/2; % bin centroids
% Cumulative mass function
CMF=cumsum(H);

yyaxis left
bar(Bc, H);
yyaxis right
hax=gca; 
hax.XGrid = 'on'; 
ylim([0,1])
hax.YTick = 0:.2:1;
yl = arrayfun(@(y)yline(hax, y,'LineStyle',hax.GridLineStyle,'Color',hax.GridColor,'Alpha', hax.GridAlpha), hax.YTick);

set(hax, 'XLimSpec', 'Tight')
h = line(Bc, CMF/length(N));
h.LineWidth = 2;
ylim([0 1])
xlim([0 5])

xlabel('Horizontal Error km')
ylabel('Event Percentage')



m_stv =mean(N);
st_stv = std(N);
m = mode(N);
MAD = mad(N); 
ND = length(N);
str = sprintf('Mean = %.1f', m_stv);
str2 = sprintf('Std = %.1f', st_stv);
str3 = sprintf('N = %d', ND);
str4 = {'Sulphur'};
str5 = sprintf('MAD = %.1f', MAD);
str6 = sprintf('Mode = %.1f', m);
text(4,0.5,str)
text(4,0.4, str2)
text(4,0.3, str3)
text(0.01, .9, str4)
text(4, 0.2, str5)
text(4, 0.1, str6)
ax=gca;




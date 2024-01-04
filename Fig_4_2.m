clear; close all; clc;

vpvs_sv = 1.78; %sawtooth model vp/vs ratio
vpvs_cv = 1.74; %challis (Schemeta) model vp/vs ratio
vpvs_suv = 1.73; %Sulphur (Brumbaugh) model vp/vs ratio


%Read in the different text files needed to run script 
suv = fopen('Velocity_Profiles/sulphur_velocity.txt','r');
sv = fopen('Velocity_Profiles/sawtooth_p_model.txt', 'r');
cv = fopen('Velocity_Profiles/challis_velocity.txt', 'r');
fix = fopen('Velocity_Profiles/ak_135.txt','r');
fid = fopen('Velocity_Profiles/ak_gradient.txt','r');

%Create textscan format for reading the decimals 
sul = textscan(suv,'%f%f','Delimiter',',');%sulphur
AK = textscan(fix, '%f%f%f', 'Delimiter',','); %ak135
G = textscan(fid, '%f%f%f', 'Delimiter',','); %ak135chal = textscan(cv, '%f%f%f', 'Delimiter', ',');%challis
saw = textscan(sv, '%f%f%f', 'Delimiter', ',');
fclose(suv);
fclose(fix);
fclose(cv);
fclose(sv);
fclose(fid);


% read in sawstooth model P and S vel
% Vp_sv = saw{1};
% depth_sv = saw{2};

%read in sulphur model P velocity and depths
Vp_suv = sul{1};
depth_suv = sul{2};

%read in challis model P velocity and depths
% % Vp_cv = chal{1};
% depth_cv = chal{2};

%%read in gradient model P and S vel
Vp_g = G{2};
Vs_g = G{3};

%read in ak135f model P and S vel
Vp_ak = AK{2};
Vs_ak = AK{3};

%Read in depths for gradient and ak models
depth_g = G{1};
depth_ak = AK{1};


% Vs_sv = Vp_sv ./ vpvs_sv; % compute the Vs velocity from the Vp velocity for stanley
Vs_suv = Vp_suv./vpvs_suv;
% Vs_cv = Vp_cv./vpvs_cv;

%Plot below 
h = figure;
% plot(Vp_sv, depth_sv, 'k','LineWidth',3); hold on; axis('ij');
% plot(Vs_sv, depth_sv, 'k--','LineWidth',3);hold on; axis('ij');
plot(Vp_g, depth_g, 'r','LineWidth',3);hold on; axis('ij');
plot(Vs_g, depth_g, 'r--','LineWidth',3);hold on; axis('ij');
plot(Vp_suv, depth_suv, 'g','LineWidth',3); hold on; axis('ij');
plot(Vs_suv, depth_suv, 'g--','LineWidth',3);hold on; axis('ij');
plot(Vp_ak, depth_ak, 'b','LineWidth',3);hold on; axis('ij');
plot(Vs_ak, depth_ak, 'b--','LineWidth',3);hold on; axis('ij');
% plot(Vp_cv, depth_cv, 'b','LineWidth',3);hold on; axis('ij');

legend({'P Gradient','S Gradient', 'P Brumbaugh', 'S Brumbaugh', 'P AK-135', 'S AK-135'}); ylabel('Depth [km]'); xlabel('Velocity [km/s]');
axis('tight'); xlim([2, 8.5]); ylim([-3.1,50]); grid on;
title('AK-135 Model vs Local Velocity Models')
set( findall( h, '-property', 'FontSize' ), 'FontSize', 18 );
set( findall( h, '-property', 'FontName' ), 'FontName', 'Helvetica' );
set( findall( h, '-property', 'FontWeight' ), 'FontWeight', 'Bold' );

set( h, 'Position', [100 100 500 700] );
set( h, 'PaperPositionMode', 'auto' );

print(h,'stanley_velocity_model','-dpng','-r300');


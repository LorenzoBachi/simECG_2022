%% ECG SIMULATOR MAIN SCRIPT: USER-DEFINED PARAMETERS, OUTPUT VISUALIZATION

% VERSION: 1.0.1 - July 2023, updated function descriptions

% For the scientific materials of this code, please see
%   L. Bachi, H. Halvaei, C. Pérez, A. Martín-Yebra, A. Petrėnas,
%   A. Sološenko, L. Johnson, V. Marozas, J. P. Martínez, E. Pueyo,
%   M. Stridh, P. Laguna, L. Sörnmo, "ECG Modeling for Simulation of
%   Arrhythmias in Time-Varying Conditions", IEEE Transactions on
%   Biomedical Engineering, 2023, https://doi.org/10.1109/TBME.2023.3288701

% This is a contribution to the previous work by
%   A. Petrėnas, V. Marozas, A. Sološenko, R. Kubilius, J. Skibarkienė,
%   J. Oster and L. Sörnmo, "Electrocardiogram modeling duringparoxysmal
%   atrial fibrillation: application to the detection of brief episodes",
%   Physiological Measurement, 2017,
%   https://doi.org/10.1088/1361-6579/aa9153

% This ECG simulator version was developed by
% Lorenzo Bachi, PhD student at Scuola Superiore Sant'Anna, Pisa,
% l.bachi@santannapisa.it
% Cristina Perez, PhD student at Universidad de Zaragoza, Spain,
% cperez@unizar.es
% Hesam Halvaei, PhD student at Lund University, Lund, Sweden
% hesam.halvaei@bme.lth.se
% together with the help of all other authors.

% For bug reporting, please email l.bachi@santannapisa.it.

% This software is licensed under GNU General Public License version 3:
% https://www.gnu.org/licenses/gpl-3.0.html This software is free and open
% source. However, please note that that if you distribute or publicly
% release software that is based on or incorporates GPL-licensed software,
% you must also distribute it under the GPL-3.0 license.

% The ECG simulator requires data files to run, which can be downloaded at
% the link found in the readme file.

% ECG is simulated @1 KHz.

clear; clc;

%% General Parameters
sigLength = 5*60;   %desired ECG length in seconds;
onlyRR = 0;         % 1 - only RR intervals are generated, 0 - multilead ECG is generated
realRRon = 0;       % 1 - real RR series are used, 0 - synthetic
realVAon = 0;       % 1 - real ventricular activity is used, 0 - synthetic
realAAon = 0;       % 1 - real atrial activity is used, 0 - synthetic
simECGdata.fs = 1000; %sampling frequency

%% Atrial Fibrillation
B_af = 0;           % AF burden. 0 - the entire signal is SR, 1 - the entire signal is AF
d_af = 100;         % Median episode length > in beats <
% default from the MIt-BIH Arrhythmia and Atrial Fibrillation Database: 85

%% Atrial Tachycardia
B_at = 0;     % AT burden
load('ATDist.mat');
%ATDist = sqrt(ATDist); ATDist = ATDist.^(1/4); %ATDist = ones(1,50); %ATDist = zeros(1,50); ATDist(1) = 1;
% ATDist represents the desired distribution of the atrial tachycardia episodes,
% in relative weights. If the sum of the weights exceeds 1, they are
% automatically rescaled in order to have a sum of 1. The first element of
% the distribution represents the probability weight of isolated APBs.
% The second element of the distribution relates to coupltes, and so on.
% The maximum length of ATDist represents the maximum desidred length of an
% atrial tachycardia episode (default: 50)
at_x = 1:50;
apb_p = [0.3,0.3,0.3,0]; %[0.3,0.3,0.3,0.1]; % probability of the four APB classes

%% Bigeminy, Trigeminy
B_bt = 0;           % Bigeminy and trigeminy burden
p_bt = [0.5, 0.5];  % differential probability of bigeminy vs trigeminy
% default: [0.72, 1-0.72]
d_bt = 8;           % Median episode length (in beats) for bigeminy and trigeminy
% default from the MIt-BIH Arrhythmia Database: 8

%% Ventricular Premature Beats
B_vpb = 0;          % VPB burden
vpb_p = [0.5,0.5,0];%[0.475,0.475,0.05]; % probability of the three VPB classes (SR only)
multiform_vpbs = 0; % if this setting is different from 0, different shapes of VPBs may be used in the same record

%% Noise Parameters
% Type of noise, a vector with the codes of all types of desired noise:
% 0 - No noise added (noise RMS = 0 mV)
% 1 - Motion artefacts
% 2 - Electrode movement
% 3 - Baseline wander
% 4 - Mixture of noises
% 5 - Bw + ma because em has residual ECG;
% 6 - Simulated muscular noise %CPerez 07/2022
% 7 - Real exercise stress test noise (from R. Bailón)
% 8 - Motion artifacts
noiseType = [3,           6,        8];
% Noise level in millivolts, a vector with each RMS level according to the
% selected noise sources
noiseRMS =  [0.15,      0.01,      2]; 
%*Note: if noiseType == 6 or 7, it can be defined a SNR value (in dB) instead of mV value.

%Motion artifacts parameters
simECGdata.MA_Prob = 0.05;  %the probability of success, i.e., spikes 
%Thumb-ECG parameter (taking into account in simulated muscular noise)
simECGdata.MA_Flag = 0;     % 0 - Holter recording  1 - Thumb-ECG

%% Exercise stress test parameters %CPerez 03/2022
%  WARNING : cannot select real atrial activity and synthetic ventricular activity
simECGdata.ESTflag = 0;                         % 1- Exercise Stress Test flag, 0 - other cases
if simECGdata.ESTflag
    simECGdata.Basal = randi([1,3],1)*60;       %Basal area before Exercise Stress Test starts, in seconds. %Cris 04/2022
    simECGdata.Exercise = randi([7,12],1)*60;   % Duration of Exercise in Exercise Stress Test in seconds. %Cris 04/2022
    simECGdata.Recovery = randi([3,5],1)*60;    %Duration of Recovery Stress Test in seconds. %Cris 04/2022
    simECGdata.Basal2 = randi([1,2],1)*60;      %Basal area before Exercise Stress Test ends. %Cris 04/2022
    
    simECGdata.Duration = simECGdata.Basal +simECGdata.Exercise + simECGdata.Recovery + simECGdata.Basal2; %Duration of Exercise Stress Test in seconds. %Cris 04/2022
    simECGdata.peak = simECGdata.Basal +simECGdata.Exercise;    % peak of Exercise Stress Test in seconds. %Cris 04/2022
    sigLength =  simECGdata.Duration;
    simECGdata.RRini = 60/randi([55,65],1);     % onset HR express according to RR in seconds. %Cris 04/2022
    simECGdata.RRpeak = 60/randi([170,185],1);  % HR express according to RR in seconds at exercise peak. %Cris 04/2022
    simECGdata.RRend = simECGdata.RRini*1.2;    % final HR express according to RR in seconds. %Cris 04/2022
    
    simECGdata.Frini = simECG_random_number(0.25, 0.35);    % onset Respiratory frequency in Hz. %Cris 04/2022
    simECGdata.Frpeak = simECG_random_number(0.65, 0.75);   % Respiratory frequency at exercise peak in Hz. %Cris 04/2022
    simECGdata.Frend = simECG_random_number(0.25, 0.35);    % final Respiratory frequency in Hz. %Cris 04/2022
end

%% Amplitude scaling
% WARNING: This feature is still in development. Use with caution! This
% amplitude factor only works for simulated ventricular and atrial
% activity.
ecg_amp = rand*10; %ECG signal amplitude scaling factor, a real number. Default is rand*10, minimum is zero

%% ECG Generator
arrhythmiaParameters.B_af = B_af;
arrhythmiaParameters.d_af = d_af;
arrhythmiaParameters.B_at = B_at;
arrhythmiaParameters.at_dist = ATDist;
arrhythmiaParameters.at_x = at_x;
arrhythmiaParameters.apb_p = apb_p;
arrhythmiaParameters.B_vpb = B_vpb;
arrhythmiaParameters.vpb_p = vpb_p;
arrhythmiaParameters.B_bt = B_bt;
arrhythmiaParameters.p_bt = p_bt./(sum(p_bt)+eps);
arrhythmiaParameters.d_bt = d_bt;
simECGdata.ecg_amp = ecg_amp;
simECGdata.multiform_vpbs = multiform_vpbs;

for nr = 1:1
    
    [simECGdata, initialParameters, annotations] = simECG_generator(sigLength,realRRon, realVAon, realAAon, noiseType, noiseRMS, onlyRR, arrhythmiaParameters, simECGdata);
    disp(['Simulation ' num2str(nr) '  completed.']);
    
end

%% Result Data and Parameters
% Initial parameters
fibFreqz = initialParameters.fibFreqz;       % Dominant fibrillatory frequency
% Data
rr = simECGdata.rr;                          % RR intervals in miliseconds
multileadECG = simECGdata.multileadECG;      % 15 lead ECG
%multileadVA = simECGdata.multileadVA;       % 15 lead ventricular activity
%multileadAA = simECGdata.multileadAA;       % 15 lead atrial activity
%multileadNoise = simECGdata.multileadNoise; % 15 lead physiological noise
%multileadECG = multileadVA + multileadAA + multileadNoise;
QRSindex = simECGdata.QRSindex;              % QRS index. Shows sample number when R peak occurs
targets_beats = simECGdata.targets_beats;    % Marks sinus rhythm, AF, premature atrial and premature ventricular beats
ecgLength = simECGdata.ecgLength;            % ECG length in samples
state_history = simECGdata.state_history;

%% Plots
lead = 2;

if ishandle(1), clf(1); end
figure(1);
ax(1)=subplot(211);
plot(cumsum(rr)./1000,rr, 'o--');
xlabel('Time [s]'); ylabel('RR [ms]');
xlim([0,sigLength]);
%axis([135, 145, 300, 1300]);

switch lead
    case 1
        
        line = 'I';
    case 2
        line = 'II';
    case 3
        line = 'III';
    case 4
        line = 'aVR';
    case 5
        line = 'aVL';
    case 6
        line = 'aVF';
    case {7,8,9,10,11,12}
        line = ['V',num2str(lead-6)];
end
ax(2)=subplot(212);
t=(0:1:length(multileadECG(7,:))-1)./1000;
plot(t,multileadECG(lead,:)');
hold on
plot(cumsum(rr)./1000,multileadECG(lead,cumsum(rr))','ro');
xlabel('Time [s]'); ylabel(['Lead ', line,' [mV]']);
xlim([0,sigLength]);
ylim([ -2.5, 2.5 ]);
%axis([135, 145, -0.2, 1]);

linkaxes(ax,'x');

set(gcf,'units','centimeters','position',[0,10,50,15]);
%set(gcf, 'Position', get(0, 'Screensize'));


if ishandle(2), clf(2); end
figure(2);
stem(cumsum(rr)./1000,targets_beats, 'r');
xlabel('Beat number');
ylabel('Type of beat');
ylim([0.75 4.25]);
yticks([1,2,3,4]); yticklabels({'N','AF','APB','V'});
xlim([0,sigLength]);
set(gcf, 'Position', get(0, 'Screensize'));

%% 12 lead plot
spacing = 3;
if ishandle(3), close(3); end
figure(3);
ax2(1)=subplot(1,2,1);
hold on;
for l=1:6
    plot(t,multileadECG(l,:)+(6-l)*spacing*ones(1,length(t)),'k');
end
axis([t(1),t(end),-spacing,6*spacing]);
yticks([0,spacing,2*spacing,3*spacing,4*spacing,5*spacing]);
yticklabels({'aVF','aVL','aVR','III','II','I'});
box on;

ax2(2)=subplot(1,2,2);
hold on;
for l=1:6
    plot(t,multileadECG(l+6,:)+(6-l)*spacing*ones(1,length(t)),'k');
end
axis([t(1),t(end),-spacing,6*spacing]);
yticks([0,spacing,2*spacing,3*spacing,4*spacing,5*spacing]);
yticklabels({'V6','V5','V4','V3','V2','V1'});
box on;

linkaxes(ax2,'x');
set(gcf, 'Position', get(0, 'Screensize'));





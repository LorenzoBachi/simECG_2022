%% ECG SIMULATOR MAIN SCRIPT: USER-DEFINED PARAMETERS, OUTPUT VISUALIZATION

% VERSION: 2022B
% A contribution to the previous works by
%     A. Petrėnas et al, 2017: Electrocardiogram modeling duringparoxysmal
% atrial fibrillation: application to the detection of brief episodes
%     R. Alcaraz et al, 2019: Reference database and performance evaluation
%of methods for extraction of atrial fibrillatory waves in the ECG

% Lorenzo Bachi, PhD student at Scuola Superiore Sant'Anna, Pisa,
% l.bachi@santannapisa.it, under the supervision of Leif Sörnmo.

% Cristina Perez, PhD student at Universidad de Zaragoza, Spain,
% cperez@unizar.es, under the supervision of Leif Sörnmo.

% ECG is simulated @1 KHz.

%Last update: CPerez 18/05/2022

clc, clear all
addpath('D:\Doctorado\2.Senales_Simuladas\Simulator_Lund\PAFsim_v2022\')
%% Initial parameters
%--> General Parameters
sigLength = 2*60;   %desired ECG length in seconds;
onlyRR = 0;         % 1 - only RR intervals are generated, 0 - multilead ECG is generated

medEpis = 100;      % Median episode length ---> in seconds <----
stayInAF = 1-log(2)/(medEpis);	% Probability to stay in AF state
AFburden = 0;     % AF burden. 0 - the entire signal is SR, 1 - the entire signal is AF

realRRon = 0;       % 1 - real RR series are used, 0 - synthetic
realVAon = 0;       % 1 - real ventricular activity is used, 0 - synthetic
realAAon = 0;       % 1 - real atrial activity is used, 0 - synthetic

APBph = 0;         % Number of median APBs per hour
ARRate = 0.85;       % Rate of atrial runs, a number between 0 and 1. 0: no runs, 1: every APB begins a run.
% ARDist represents the desired distribution of the atrial runs, in
% relative weights. If the sum of the weights exceeds 1, they are
% automatically rescaled in order to have a sum of 1. The first element of
% the distribution represents the probability weight of a run of TWO APBs.
% The maximum length of ARDist represents the maximum desidred length of an
% atrial run minus one (ARDist(1) refers to couplets).
ARDist = ones(1,49); %1 + (1:49)./20; %ones(1,49); %1./((1:49).^2);
VPBph = 0;         % Number of median VPBs per hour


%--> Noise Parameters
noiseType = 0;      % Type of noise. 4 - a mixture of all noises; 5 - bw + ma because em has residual ECG; % Type 6 - real EST
noiseRMS = 0.25;    % Noise level in milivolts. 

%--> Exercise stress test parameters %CPerez 03/2022
ecgParameters.ESTflag = 0;     % 1- Exercise Stress Test flag, 0 - other cases
if ecgParameters.ESTflag
    ecgParameters.Basal = 3*60;      %Basal area before Exercise Stress Test starts, in seconds. %Cris 04/2022
    ecgParameters.Exercise = 8*60;    % Duration of Exercise in Exercise Stress Test in seconds. %Cris 04/2022
    ecgParameters.Recovery = 3.3*60; %Duration of Recovery Stress Test in seconds. %Cris 04/2022
    ecgParameters.Basal2 = 2*60; %Basal area before Exercise Stress Test ends. %Cris 04/2022
    
    ecgParameters.Duration = ecgParameters.Basal +ecgParameters.Exercise + ecgParameters.Recovery + ecgParameters.Basal2; %Duration of Exercise Stress Test in seconds. %Cris 04/2022
    ecgParameters.peak = ecgParameters.Basal +ecgParameters.Exercise;    % peak of Exercise Stress Test in seconds. %Cris 04/2022
    
    ecgParameters.RRini = 60/60;    % onset HR express according to RR in seconds. %Cris 04/2022
    ecgParameters.RRpeak = 60/180;    % HR express according to RR in seconds at exercise peak. %Cris 04/2022
    ecgParameters.RRend = 60/70;    % final HR express according to RR in seconds. %Cris 04/2022
    
    ecgParameters.Frini = 0.3;    % onset Respiratory frequency in Hz. %Cris 04/2022
    ecgParameters.Frpeak = 0.7;    % Respiratory frequency at exercise peak in Hz. %Cris 04/2022
    ecgParameters.Frend = 0.3;    % final Respiratory frequency in Hz. %Cris 04/2022
end

% Note: cannot select real atrial activity and synthetic ventricular activity

%% ECG generator 
%resuldir = 'D:\investigacion\QTadaptation\matlab\simulated_ECGs_QTadapt\';
clc

arrhythmiaParameters.AFburden = AFburden;
arrhythmiaParameters.stayInAF = stayInAF;
arrhythmiaParameters.APBph = APBph;
arrhythmiaParameters.ARRate = ARRate;
arrhythmiaParameters.ARDist = ARDist;
arrhythmiaParameters.VPBph = VPBph;


for nr = 1:1
[simECGdata, initialParameters, annotations, ecgParameters] = simECG_generator(sigLength,realRRon, realVAon, realAAon, noiseType, noiseRMS, onlyRR, arrhythmiaParameters, ecgParameters);
fprintf('finish!\n')

%cd(resuldir)
fa=1000;
ecgLength=simECGdata.ecgLength;
Duration=ecgLength./fa;
sinal=simECGdata.multileadECG([1,2,7:12],:)';

% ecgnr = 'sim_02_synt';
% ecgnr = 'sim_Example_nTV';
% save(['.\Simulated_results\' ecgnr '.mat'], 'simECGdata','sinal','fa','Duration', 'ecgnr','sinal_clean','ecgParameters'); fprintf('saved!\n')


% ecgnr = 'sim_02_n1_ECG-LR';
% ecgnr = 'sim_65_ECG-LR';
% save(['.\DelayQT_sim\Results\' ecgnr '.mat'], 'simECGdata','sinal','fa','Duration', 'ecgnr', 'sinal_clean'); fprintf('saved!\n')


% save(strcat(ecgnr,'.mat'), 'simPAFdata','sinal','fa','Duration', 'ecgnr') 
% save(['.\Results\sim_' num2str(ns) '.mat'], 'simECGdata','-v7.3')
%cd ..
end


%% Returned data and initial parameters
% Initial parameters
fibFreqz = initialParameters.fibFreqz;      % Dominant fibrillatory frequency
% Data
rr = simECGdata.rr;                         % RR intervals in miliseconds
multileadECG = simECGdata.multileadECG;     % 15 lead ECG
multileadVA = simECGdata.multileadVA;       % 15 lead ventricular activity
multileadAA = simECGdata.multileadAA;       % 15 lead atrial activity
multileadNoise = simECGdata.multileadNoise; % 15 lead physiological noise
QRSindex = simECGdata.QRSindex;             % QRS index. Shows sample number when R peak occurs
targets_beats = simECGdata.targets_beats;   % Marks sinus rhythm, AF, premature atrial and premature ventricular beats
pafBoundaries = simECGdata.pafBoundaries;   % Start and the end of each PAF episode
pafEpisLength = simECGdata.pafEpisLength;   % Length (in beats) of each PAF episode
ecgLength = simECGdata.ecgLength;           % ECG length in samples

%% Plots
if ishandle(1), close(1); end
figure(1)
ax(1)=subplot(211);
plot(cumsum(rr)./1000,rr, 'o--');
xlabel('Time [s]'); ylabel('RR [ms]');
% xlim([0,sigLength]);
%axis([135, 145, 300, 1300]);

l=2;
switch l
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
        line = ['V',num2str(l-6)];
end
ax(2)=subplot(212);
t=(0:1:length(simECGdata.multileadECG(7,:))-1)./1000;
plot(t,simECGdata.multileadECG(10,:)');
hold on
plot(cumsum(rr)./1000,simECGdata.multileadECG(l,cumsum(rr))','ro');
xlabel('Time [s]'); ylabel(['Lead ', line,' [mV]']);
xlim([0,sigLength]);
%axis([135, 145, -0.2, 1]);

linkaxes(ax,'x');

set(gcf, 'Position', get(0, 'Screensize'));

if ishandle(2), close(2); end
figure(2)
hold on;
plot(rr);
plot(targets_beats, 'r');
hold off;
xlabel('Beat number');
ylabel('Type of beat');
ylim([0.75 3.25]);
yticks([1,2,3]); yticklabels({'N','AF','APB'});
set(gcf, 'Position', get(0, 'Screensize'));


%% LEGACY CODE
% 
% %% PAF detection
% tarECG = multileadECG(7,:); % Target lead. Must be lead V1
% refECG = multileadECG(1,:); % Reference lead. Lead I, or V6.
% if ishandle(3), close(3); end
% figure(3);
% hold on;
% plot(tarECG+3);
% plot(refECG, 'r');
% hold off;
% legend('Target lead', 'Reference lead');
% 
% %% ECG preprocessing
% % Baseline removal < 0.5 Hz
% b = 0.997781*[1 -2 1];
% a = [1 -1.9955571 0.99556697];
% % Noise removal > 40 Hz
% h  = fdesign.lowpass('N,F3dB', 5, 40, 1000);
% Hd = design(h, 'butter');
% tarECG = filtfilt(b,a, tarECG);
% tarECG = filtfilt(Hd.sosMatrix,Hd.ScaleValues,tarECG);
% refECG = filtfilt(b,a, refECG);
% refECG = filtfilt(Hd.sosMatrix,Hd.ScaleValues,refECG);

% f waves extraction using ESN
% reservoirSize = 100;      
% QRSdetectionOn = 1;       % QRS detection performed inside the function func_fWaveExtractionESN
% rIndex = [];
% [estAA, estVA, rIndex, rrInt, tarECG, refECG, refECGqrs] = func_fWaveExtractionESN(reservoirSize, tarECG, refECG, QRSdetectionOn, rIndex);
% rr = rrInt/1000;  % RR interval length in seconds
% 
% PAF detection
% [R, ~] = func_LowComplexityDetector(rr, 0.1, 0.725);
% N = func_noise_index(rIndex,estAA);    
% [F, P] = func_P_f_waves_index(rIndex,estAA);   %both indexes are updated compared to occult paper
% [O, Othr] = func_occult_AF_detector(R, P, F, N); % fuzzy based AF detector
%  
% figure(1)
% subplot(2,1,1)
% hold on
% plot(tarECG)
% plot(estAA+mean(tarECG)+0.3, 'r')
% stairs(QRSindex/4, 2*targets_SR_AF-1, 'b')
% hold off
% xlabel('Sample number')
% ylabel('Amplitude')
% title('ECG & extracted f-waves')
% 
% subplot(2,1,2)
% hold on
% plot(O, 'Linewidth', 2, 'Color', [0 1 0])
% % plot(R, 'r')
% line([0 length(O)], [0.6 0.6], 'Color', [1 0 0])
% plot(targets_SR_AF, 'b')
% hold off
% ylim([0 1])
% legend('O')
% xlabel('RR interval number')
% title('Detector output')




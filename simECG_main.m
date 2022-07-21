%% ECG SIMULATOR MAIN SCRIPT: USER-DEFINED PARAMETERS, OUTPUT VISUALIZATION

% VERSION: 2022C
% A contribution to the previous works by
%     A. Petrėnas et al, 2017: Electrocardiogram modeling duringparoxysmal
% atrial fibrillation: application to the detection of brief episodes
%     R. Alcaraz et al, 2019: Reference database and performance evaluation
%of methods for extraction of atrial fibrillatory waves in the ECG

% Lorenzo Bachi, PhD student at Scuola Superiore Sant'Anna, Pisa,
% l.bachi@santannapisa.it, under the supervision of Leif Sörnmo.

% Cristina Perez, PhD student at Universidad de Zaragoza, Spain,
% cperez@unizar.es, under the supervision of Leif Sörnmo.

% Hesam Halvaei, PhD student at Lund University, Lund, Sweden
% hesam.halvaei@bme.lth.se


% ECG is simulated @1 KHz.

clear; clc;

%% Initial parameters
%--> General Parameters
sigLength = 3*60;   %desired ECG length in seconds;
onlyRR = 0;         % 1 - only RR intervals are generated, 0 - multilead ECG is generated
realRRon = 0;       % 1 - real RR series are used, 0 - synthetic
realVAon = 0;       % 1 - real ventricular activity is used, 0 - synthetic
realAAon = 0;       % 1 - real atrial activity is used, 0 - synthetic
ecgParameters.fs = 1000; %sampling frequency

%--> Atrial fibrillation
medEpis = 100;      % Median episode length > in beats <
stayInAF = 1-log(2)/(medEpis);	% Probability to stay in AF state
AFburden = 0;     % AF burden. 0 - the entire signal is SR, 1 - the entire signal is AF

%--> Atrial tachycardia
APBph = 10;         % Number of APBs per hour
load('ATDist.mat'); %comment for custom probability distribution
% ATDist represents the desired distribution of the atrial tachycardia episodes,
% in relative weights. If the sum of the weights exceeds 1, they are
% automatically rescaled in order to have a sum of 1. The first element of
% the distribution represents the probability weight of isolated APBs.
% The second element of the distribution relates to coupltes, and so on.
% The maximum length of ATDist represents the maximum desidred length of an
% atrial tachycardia episode (default: 50)
ATDist = sqrt(ATDist);%ones(1,50);

% Ventricular premature beats
VPBph = 0;         % Number of VPBs per hour

% Bigeminy, trigeminy
BT_r = 0; % rate of bigeminy and trigeminy
BT_p = [0, 0]; % differential probability of bigeminy vs trigeminy
%Setting both probabilities to zero deactivates the BT state
BT_medEpis = 0;    % Median episode length (in beats) for bigeminy and trigeminy

%--> Noise Parameters
noiseType = 6;      % Type of noise. 
% 0 - no noise;
% 1 - motion artifact;
% 2 - electrode movement
% 3 - baseline wander;
% 4 - a mixture of all noises;
% 5 - bw + ma because em has residual ECG;
% 6 - Simulated Muscular Noise;

noiseRMS = 0.01;    % Noise level in millivolts. 

%--> Exercise stress test parameters %CPerez 03/2022
ecgParameters.ESTflag = 1;     % 1- Exercise Stress Test flag, 0 - other cases
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
clc
arrhythmiaParameters.AFburden = AFburden;
arrhythmiaParameters.stayInAF = stayInAF;
arrhythmiaParameters.APBph = APBph;
arrhythmiaParameters.ATDist = ATDist;
arrhythmiaParameters.BT_r = BT_r;
arrhythmiaParameters.BT_p = BT_p./(sum(BT_p)+eps);
arrhythmiaParameters.BT_medEpis = BT_medEpis;
if arrhythmiaParameters.APBph >0
    arrhythmiaParameters.ATDist = arrhythmiaParameters.ATDist./ sum(arrhythmiaParameters.ATDist);
end
arrhythmiaParameters.VPBph = VPBph;


for nr = 1:1

[simECGdata, initialParameters, annotations, ecgParameters] = simECG_generator(sigLength,realRRon, realVAon, realAAon, noiseType, noiseRMS, onlyRR, arrhythmiaParameters, ecgParameters);
fprintf('Simulation completed.\n')

end


%% Returned data and initial parameters
% Initial parameters
fibFreqz = initialParameters.fibFreqz;      % Dominant fibrillatory frequency
% Data
rr = simECGdata.rr;                         % RR intervals in miliseconds
multileadECG = simECGdata.multileadECG;     % 15 lead ECG
%multileadVA = simECGdata.multileadVA;       % 15 lead ventricular activity
%multileadAA = simECGdata.multileadAA;       % 15 lead atrial activity
%multileadNoise = simECGdata.multileadNoise; % 15 lead physiological noise
QRSindex = simECGdata.QRSindex;             % QRS index. Shows sample number when R peak occurs
targets_beats = simECGdata.targets_beats;   % Marks sinus rhythm, AF, premature atrial and premature ventricular beats
pafBoundaries = simECGdata.pafBoundaries;   % Start and the end of each PAF episode
pafEpisLength = simECGdata.pafEpisLength;   % Length (in beats) of each PAF episode
ecgLength = simECGdata.ecgLength;           % ECG length in samples
state_history = simECGdata.state_history;

%% Plots
if ishandle(1), close(1); end
figure(1);
ax(1)=subplot(211);
plot(cumsum(rr)./1000,rr, 'o--');
xlabel('Time [s]'); ylabel('RR [ms]');
% xlim([0,sigLength]);
%axis([135, 145, 300, 1300]);

l=11;
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
t=(0:1:length(multileadECG(7,:))-1)./1000;
plot(t,multileadECG(l,:)');
hold on
plot(cumsum(rr)./1000,multileadECG(l,cumsum(rr))','ro');
xlabel('Time [s]'); ylabel(['Lead ', line,' [mV]']);
xlim([0,sigLength]);
%axis([135, 145, -0.2, 1]);

linkaxes(ax,'x');

set(gcf, 'Position', get(0, 'Screensize'));

if ishandle(2), close(2); end
figure(2);
stem(targets_beats, 'r');
xlabel('Beat number');
ylabel('Type of beat');
ylim([0.75 4.25]);
yticks([1,2,3,4]); yticklabels({'N','AF','APB','V'});
set(gcf, 'Position', get(0, 'Screensize'));

% %% Provisional code for the paper figures
% clear simECGdata;
% spacing = 2;
% if ishandle(3), close(3); end
% figure(3);
% hold on;
% limits = [20, 30]; %seconds
% idx = round(limits(1)*1000) +1: round(limits(2)*1000);
% tf = t( idx );
% plot(tf-limits(1),multileadECG(7,idx)','k');
% plot(tf-limits(1),multileadECG(2,idx)'+spacing*ones(1,length(idx)),'k');
% set(gcf,'units','centimeters','position',[15,15,18,10]);
% set(gca,'FontName','Times','Fontsize',10);
% box on;
% set(gca,'ytick',[-0.5,0,0.5,1.5,2,2.5]); set(gca,'yticklabel',{'-0.5','0','0.5','-0.5','0','0.5'});

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




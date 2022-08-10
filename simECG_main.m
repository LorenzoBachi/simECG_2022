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
sigLength = 2*60;   %desired ECG length in seconds;
onlyRR = 0;         % 1 - only RR intervals are generated, 0 - multilead ECG is generated
realRRon = 0;       % 1 - real RR series are used, 0 - synthetic
realVAon = 0;       % 1 - real ventricular activity is used, 0 - synthetic
realAAon = 0;       % 1 - real atrial activity is used, 0 - synthetic
ecgParameters.fs = 1000; %sampling frequency

%--> Atrial fibrillation
medEpis = 0;      % Median episode length > in beats <
AFburden = 0;     % AF burden. 0 - the entire signal is SR, 1 - the entire signal is AF

%--> Atrial tachycardia
APBph = 0;         % Number of APBs per hour
load('ATDist.mat'); %comment for custom probability distribution
% ATDist represents the desired distribution of the atrial tachycardia episodes,
% in relative weights. If the sum of the weights exceeds 1, they are
% automatically rescaled in order to have a sum of 1. The first element of
% the distribution represents the probability weight of isolated APBs.
% The second element of the distribution relates to coupltes, and so on.
% The maximum length of ATDist represents the maximum desidred length of an
% atrial tachycardia episode (default: 50)
%ATDist = %sqrt(ATDist); %ones(1,50);

% Ventricular premature beats
VPBph = 0;         % Number of VPBs per hour

% Bigeminy, trigeminy
BT_r = 0; % rate of bigeminy and trigeminy
BT_p = [0.5, 0.5]; % differential probability of bigeminy vs trigeminy
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
    ecgParameters.Basal = randi([2,4],1)*60;      %Basal area before Exercise Stress Test starts, in seconds. %Cris 04/2022
    ecgParameters.Exercise = randi([7,10],1)*60;    % Duration of Exercise in Exercise Stress Test in seconds. %Cris 04/2022
    ecgParameters.Recovery = randi([3,5],1)*60; %Duration of Recovery Stress Test in seconds. %Cris 04/2022
    ecgParameters.Basal2 = randi([1,2],1)*60; %Basal area before Exercise Stress Test ends. %Cris 04/2022
    
    ecgParameters.Duration = ecgParameters.Basal +ecgParameters.Exercise + ecgParameters.Recovery + ecgParameters.Basal2; %Duration of Exercise Stress Test in seconds. %Cris 04/2022
    ecgParameters.peak = ecgParameters.Basal +ecgParameters.Exercise;    % peak of Exercise Stress Test in seconds. %Cris 04/2022
    sigLength =  ecgParameters.Duration;
    ecgParameters.RRini = 60/randi([55,70],1);    % onset HR express according to RR in seconds. %Cris 04/2022
    ecgParameters.RRpeak = 60/randi([160,180],1);    % HR express according to RR in seconds at exercise peak. %Cris 04/2022
    ecgParameters.RRend = ecgParameters.RRini*1.2;    % final HR express according to RR in seconds. %Cris 04/2022
    
    ecgParameters.Frini = (0.35 - 0.25)*rand(1) + 0.25;    % onset Respiratory frequency in Hz. %Cris 04/2022
    ecgParameters.Frpeak = (0.75 - 0.65)*rand(1) + 0.65;    % Respiratory frequency at exercise peak in Hz. %Cris 04/2022
    ecgParameters.Frend = (0.35 - 0.25)*rand(1) + 0.25;    % final Respiratory frequency in Hz. %Cris 04/2022
end

% Note: cannot select real atrial activity and synthetic ventricular activity

%% ECG generator
clc
if medEpis < 1, medEpis = 1; end
if BT_medEpis < 1, medEpis = 1; end
arrhythmiaParameters.AFburden = AFburden;
stayInAF = 1-log(2)/(medEpis);	% Probability to stay in AF state
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
xlim([0,sigLength]);
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

%% Provisional code for the paper figures
% clear simECGdata;
% spacing = 2;
% if ishandle(4), close(4); end
% figure(4);
% hold on;
% limits = [288, 294]; %seconds
% idx = round(limits(1)*1000) +1: round(limits(2)*1000 +1);
% tf = t( idx );
% lh1 = plot(tf-limits(1),multileadECG(7,idx)','k');
% lh2 = plot(tf-limits(1),multileadECG(2,idx)'+spacing*ones(1,length(idx)),'k');
% axis([tf(1)-limits(1),tf(end)-limits(1),-1.5,3.5]);
% set(gcf,'units','centimeters','position',[15,15,18,6]);
% set(gca,'FontName','Times','Fontsize',10);
% box on;
% gridecg2(0.25,0.25);
% %set(gca,'ytick',[-0.5,0,0.5,1.5,2,2.5]); set(gca,'yticklabel',{'-0.5','0','0.5','-0.5','0','0.5'});
% set(gca,'ytick',[0,2]); set(gca,'yticklabel',{'V1','II'});
% set(gca,'TickLength',[0 0])
% delete(lh1); delete(lh2);
% lh1 = plot(tf-limits(1),multileadECG(7,idx)','k');
% lh2 = plot(tf-limits(1),multileadECG(2,idx)'+spacing*ones(1,length(idx)),'k');
% paper_example.multileadECG = multileadECG;
% paper_example.limits = limits;
%save('.\ECG simulator\2022 07 Paper examples\paper_example5.mat','paper_example');





function [rr, annTime, annBeats, annRhythm, targets_beats, pafBoundaries, pafEpisodeLength, ecgParameters, state_history] = simECG_global_rr_intervals(sigLength, fibFreqz, realRRon, arrhythmiaParameters, ecgParameters)
% [] = simECG_global_rr_intervals() generates the global RR series, which
% may include sinus rhythm, atrial fibrillation or other arrhythmias.

% The dominant length of PAF episodes, and the percentage of the total time AF is
% present (AF burden) is determined by a first-order two-state Markov chain.
%
% An atrioventricular node model in which the ventricles are assumed to be
% activated by the arriving atrial impulses according to a Poisson process
% is used to generate realistic RR series during an episode of AF

% Ventricular rhythm during SR is simulated using RR interval generator
% proposed by in which both the impact of parasympathetic stimulation
% (respiratory sinus arrhythmia) and baroreflex regulation (Mayer waves) is
% modeled by a bimodal power spectrum

% In case of real data the entire MIT-BIH Normal Sinus Rhythm database,
% consisting of 18 long-term ECG recordings, was taken as a basis for the set
% of 18 RR interval series of sinus rhythm.
%
% The Long Term Atrial Fibrillation database was used to compose a set
% of AF rhythm. In total, 69 RR interval series were extracted out of
% 84 long-term ECG recordings. The remaining records were excluded since
% they did not meet the criteria of ?5000 RR intervals of AF.

% With respect to the previous contributions to the simulator, in this
% version the RR series is generated from four separate pacing
% processes: the sinus node pacing, the atrial ectopic pacing, the
% ventricular ectopic pacing and the atrial fibrillation pacing. The four
% stimulation processes are combined in a single RR series. The pacing
% processes may behave in different manners depending on the dominating
% rhythm.

% Additionally, the simulator outputs beat and rhythm annotations in the
% same format of the MIT-BIH Arrhythmia Database.

% The desired ECG signal length in number of RR intervals 
% rrLength = round(sigLength/0.2);
rrLength = 5*sigLength;

% signal length in ms
sigLengthMs = sigLength * 1000;
ecgParameters.Fr = [];

% fetching simulation parameters
AFburden = arrhythmiaParameters.AFburden;
stayInAF = arrhythmiaParameters.stayInAF;
APBph = arrhythmiaParameters.APBph;
ATDist = arrhythmiaParameters.ATDist;
VPBph = arrhythmiaParameters.VPBph;
BT_p = arrhythmiaParameters.BT_p;
BT_medEpis = arrhythmiaParameters.BT_medEpis;
% beat codes for annotation production
beatCodes = ['N','A','V','+']; % please note: if new beat types are added, '+' should always be last
rhythmCodes = {'(N','(AFIB','(SVTA','(B','(T'};

% Choice of rhythm to generate
if ecgParameters.ESTflag == 1 %Lorenzo 06/2022
    rhythmType = 3; % SR and EST
    if AFburden ~= 0
        AFburden = 0;
        disp('Warning: excercise stress test activated, AF burden set to zero.');
    end
elseif AFburden == 0
    rhythmType = 0; % SR
elseif AFburden == 1
    rhythmType = 1; % AF
elseif AFburden > 0 && AFburden < 1
    rhythmType = 2; % PAF
else
    error('AF burden must be a value between 0 and 1')
end

disp('Generating RR intervals ...');

%% RR generation
rr_af = [];
switch rhythmType  % 0 - sinus rhythm, 1 - AF, 2 - PAF

    case 0 % The entire rhythm is SR
        if realRRon == 1 % Use real RR series
            rr_sr = simECG_get_real_RR_intervals(0, rrLength); % sinus rhythm rr series (real)
            hrMean = (1./diff(rr_sr))*60 ;
        else % Generate sinus activity with the McSharry spectral model
            [rr_sr,hrMean,ecgParameters] = simECG_generate_sinus_rhythm(rrLength,ecgParameters); % sinus rhythm rr series
%             rr_sr = rr_sr(1:rrLength);
        end
        
    case 1 % The entire rhythm is AF
        if realRRon == 1 % Use real RR series
            rr_af = simECG_get_real_RR_intervals(1, rrLength); %atrial fibrillation rr series
        else % Generate atrial fibrillation activity with the AV model of V. Corino
            rr_af = simECG_generate_AF_intervals(fibFreqz,rrLength); %atrial fibrillation rr series
            rr_af = rr_af(1:rrLength);
        end
        
    case 2 % PAF
        
        % Generate sinus node pacing activity
        srLength = rrLength;
        if realRRon == 1 % Use real RR series
            rr_sr = simECG_get_real_RR_intervals(0, srLength); % sinus rhythm rr series (real)
            hrMean = (1./diff(rr_sr))*60 ;
        else % Use simulated RR series
            [rr_sr,hrMean,ecgParameters] = simECG_generate_sinus_rhythm(srLength,ecgParameters); % sinus rhythm rr series
            sigLength = ceil(ecgParameters.Duration);
%             rr_sr = rr_sr(1:srLength);
        end
        
        % Generate atrial fibrillation pacing activity
        afLength = round(AFBurden*rrLength);
        if realRRon == 1 % Use real RR series
            rr_af = simECG_get_real_RR_intervals(1, afLength);  %atrial fibrillation rr series
        else
            rr_af = simECG_generate_AF_intervals(fibFreqz,afLength); %atrial fibrillation rr series
            rr_af = rr_af(1:afLength);
        end
        
        % Make sure that the average RR value during SR is larger than that in AF
        if mean(rr_sr) < mean(rr_af)
            rr_sr = rr_sr + (mean(rr_af) - mean(rr_sr));
        end
        
        
    case 3 % The entire rhythm is SR but for EST signal %CPerez 06/2021
        srLength = [];
        if realRRon == 1 % Use real RR series
            [rr_sr, ecgParameters.peak, fa ecgParameters.ecgnr] = simECG_get_real_RR_intervals(2,rrLength); % If opt == 1 - AF, If opt == 2 - SR for stress test
             
            L = 300;
            if  realRRon %CPerez 04/2022. To take into account the influence of the "QT memory" before starting the EST
                rep = ceil(L/rr_sr(1));
                rr_sr = [repmat(rr_sr(1),rep,1); rr_sr];
            end
             
            dHR = 60./rr_sr; %instantaneous HR
            timebeats = [0; cumsum(rr_sr(2:end))]; %In seconds
            [pos, dHRu] = resamp(dHR,timebeats.*fa,fa,4);%4Hz resample
            pos = pos./fa;
            [b,a] = butter(2,0.03); %Fc = 0.03Hz
            hrMeanu = filtfilt(b,a,dHRu);
            hrMean = interp1(pos, hrMeanu, timebeats);
            ecgParameters.Duration = sum(rr_sr); %in seconds
            
            %Recalculate the exercise peak
            posPeak = sort(find(round(rr_sr,3)==min(round(rr_sr,3)))); %peak position
            if length(posPeak)>1
                posPeak = posPeak(round(length(posPeak)/2));
            end
            ecgParameters.peak = sum(rr_sr(1:posPeak));
            
            
            ecgParameters.Frini = [];
            ecgParameters.Frpeak = [];
            ecgParameters.Frend = [];
            
            sigLength = ceil(ecgParameters.Duration);
            sigLengthMs = sigLength * 1000;
            
        else % Use synthetic RR series CPerez 04/2022
            [rr_sr, hrMean, ecgParameters] = simECG_generate_sinus_rhythm(rrLength, ecgParameters);
            sigLength = fix(ecgParameters.Duration); %in sec
            sigLengthMs = sigLength * 1000; %in msec
        end
        rhythm_states = ones(1,sigLength);
        rr_af = Inf;
end

%% Markov chain transition matrix
% State names
%  state 1: SR
%  state 2: AF
%  state 3: AT
%  state 4: BT
%  states 5, 6, 7: VPBs
sN = {'1';'2';'3';'4';'5';'6';'7'};

% Probabilities from previous version
goToNS = 1-stayInAF;
goToAF = (goToNS*AFburden)/(1-AFburden);

p_SR_VPB = VPBph / (60 * hrMean);
p_AF_VPB = VPBph / (60 * hrMean);
p_AT_VPB = 0;
p_AF_AF = stayInAF;
p_SR_AF = goToAF;
p_AF_SR = 1 - p_AF_VPB - p_AF_AF;
p_SR_AT = APBph / (60 * hrMean);
p_AT_SR = 1;
p_SR_BT = (BT_p(1) | BT_p(2))*VPBph / (60 * hrMean);
p_BT_SR = log(2) / BT_medEpis;
p_SR_SR = 1 - p_SR_AF - p_SR_AT - p_SR_VPB - p_SR_BT;

transM = zeros(7);
transM(1,1) = p_SR_SR;
transM(1,2) = p_SR_AF;
transM(2,2) = p_AF_AF;
transM(2,1) = p_AF_SR;
transM(1,3) = p_SR_AT;
transM(3,1) = p_AT_SR;
transM(1,4) = p_SR_BT;
transM(4,1) = p_BT_SR;
transM(4,4) = 1 - p_BT_SR;
transM(1,5) = p_SR_VPB;
transM(5,1) = 1;
transM(2,6) = p_AF_VPB;
transM(6,2) = 1;
transM(3,7) = p_AT_VPB;
transM(7,3) = 1;

%trans = [stayInNS,1-stayInNS;1-stayInAF, stayInAF];
%[~,rhythm_states] = hmmgenerate(sigLength,trans,[0; 1]);

% turn SR rr intervals into progressive time stamps
sp = round(cumsum(rr_sr)*1000);
% starting (previous) state
ps = 1;
% time counter
t = 0;
% rr, target beats and state history vectors initialization
targets_beats = zeros(1,rrLength);
rr = zeros(1,rrLength);
state_history = zeros(1,rrLength);
% beat counter
k = 0;
% annotation vectors initialization
annTime = zeros(1,2*rrLength);
annBeats = zeros(1,2*rrLength);
annRhythm = cell(1,2*rrLength);
kann = 0;
% rhythm check, used to annotate rhythm changes
rc = 0;
% rhythm counters
c_SR = 1;
c_AF = 1;
% flags used to determine wheter an AT or BT episode has begun
at_flag = 0;
bt_flag = 0;
% RR series generation loop
while t<=sigLengthMs
    % the "hmmgenerate" command assumes that the starting state of the HMM
    % is always the first one. However, in the simulator, we generate one
    % state at a time. This means that the initial state needs to be
    % swapped with state 1 after each state is drawn from the Markov
    % chain.
    states = (1:7)';
    temp = states(ps);
    states(ps) = 1;
    states(1) = temp;
    [~,ss] = hmmgenerate(1,transM(states,states),states,'Statenames',sN(states));
    s = str2double(ss);
    k = k + 1;
    state_history(k) = s;
    switch s
        case 1
            % SR
            % insert normal beat
            rr(k) = rr_sr(c_SR);
            % check if rhythm annotation should change
            if rc~=1
                % rhythm annotation - N
                kann=kann+1;
                annTime(kann) = t + round((rr(k)*1000)/2);
                annBeats(kann) = beatCodes(end);
                annRhythm{kann}= rhythmCodes(1);
                % rhythm check
                rc = 1;
            end
            % update time counter
            t = t + round(rr(k)*1000);
            % sp series is adjusted
            sp(c_SR:end) = sp(c_SR:end) + t-sp(c_SR);
            % increase SR counter
            c_SR = c_SR + 1;
            % beat annotation
            kann=kann+1;
            annTime(kann) = round(t);
            annBeats(kann) = beatCodes(1);
            annRhythm{kann} = [];
            % update targets array
            targets_beats(k) = 1;
        case 2
            % AF
            % insert af beat from Corino's AV model
            rr(k) = rr_af(c_AF);
            % check if rhythm annotation should change
            if rc~=2
                % rhythm annotation - AFIB
                kann=kann+1;
                annTime(kann) = t + round((rr(k)*1000)/2);
                annBeats(kann) = beatCodes(end);
                annRhythm{kann}= rhythmCodes(2);
                % rhythm check
                rc = 2;
            end
            % update time counter
            t = t + round(rr(k)*1000);
            % increase AF counter
            c_AF = c_AF + 1;
            % beat annotation
            kann=kann+1;
            annTime(kann) = round(t);
            annBeats(kann) = beatCodes(2);
            annRhythm{kann} = [];
            % update targets array
            targets_beats(k) = 2;
        case 3
            % AT
            if at_flag == 0
                % sampling the tachycardia episode length distribution
                x = find(rand<=cumsum(ATDist));
                % duration is the first nonzero index of rand<=cumDist
                d = x(1);
                if d == 1
                    % Isolated APB
                    apbtype = randi([1,4]);
                    switch apbtype
                        case 1 % APBs with sinus reset
                            % prematurity factor
                            beta_p = 0.55 + rand*0.4;
                            % rr interval of the apb
                            rrtemp = rr_sr(c_SR) * beta_p;
                            % No compensatory pause
                        case 2 % APBs with delayed sinus reset
                            % prematurity factor
                            beta_p = 0.55 + rand*0.4;
                            % rr interval of the apb
                            rrtemp = rr_sr(c_SR) * beta_p;
                            % delay factor
                            beta_f = 1.1 + rand*0.25;
                            % apb compensatory pause
                            rr_sr(c_SR) = rr_sr(c_SR) * beta_f;
                        case 3 % APBs with full compensatory pause
                            % prematurity factor
                            beta_p = 0.55 + rand*0.4;
                            % rr interval of the apb
                            rrtemp = rr_sr(c_SR) * beta_p;
                            % apb compensatory pause
                            rr_sr(c_SR) = (2*rr_sr(c_SR)) - rrtemp;
                        case 4 % Interpolated APBs
                            % prematurity factor
                            beta_p = 0.55 + rand*0.1;
                            % rr interval of the apb
                            rrtemp = rr_sr(c_SR) * beta_p;
                            % apb compensatory pause
                            rr_sr(c_SR) = rr_sr(c_SR) - rrtemp;
                    end
                else
                    % Tachycardia episode begins
                    bpm_at = ( rand*(2-1.1) + 1.1 ) * hrMean;
                    % excessive rates are discarded
                    while bpm_at > 200
                        bpm_at = ( rand*(2-1.1) + 1.1 ) * hrMean;
                    end
                    beta_at = hrMean / bpm_at;
                    % prematurity factor
                    beta_p = 0.55 + rand*0.4;
                    rrtemp = rr_sr(c_SR) * beta_p;
                    % delay factor
                    beta_f = rand*(1.5-0.7) + 0.7;
                    % apb compensatory pause
                    c_SR = c_SR + 1;
                    rr_sr(c_SR) = rr_sr(c_SR) * beta_f;
                    % tachycardia annotation
                    if rc~=3
                        % rhythm annotation - SVTA
                        kann=kann+1;
                        annTime(kann) = t + round(0.5*(rrtemp*1000));
                        annBeats(kann) = beatCodes(end);
                        annRhythm{kann}= rhythmCodes(3);
                        rc=3;
                    end
                    % updating the transition probabilities in the Markov
                    % chain
                    transM(3,7) = p_SR_VPB;
                    transM(3,3) = 1 - p_SR_VPB;
                    transM(3,1) = 0;
                    % the episode has begun
                    at_flag = d-1;
                end
            else
                % tachycardia episode continues
                rrtemp = 0;
                while rrtemp<0.3
                    % beat to beat variation
                    delta = (rand*0.1*2) - 0.1 ; % min(max( (1 + randn*0.1), 0.9 ), 1.1 ) ;
                    % current rr interval of the episode
                    f = beta_at;
                    rrtemp = f * rr_sr(c_SR-1) + delta ;
                end
                % is the tachycardia episode still ongoing?
                at_flag = at_flag - 1;
                if at_flag <= 0
                    % the tachycardia episode is over
                    % transition probabilities return to initial values
                    transM(3,7) = 0;
                    transM(3,3) = 0;
                    transM(3,1) = 1;
                end
            end
            % insert atrial beat
            rr(k) = rrtemp;
            % update targets array
            targets_beats(k) = 3;
            % update time counter
            t = t + round(rr(k)*1000);
            % sinus rhythm rr series is delayed
            sp(c_SR:end) = sp(c_SR:end) + t + round( (rrtemp + rr_sr(c_SR))* 1000 ) - sp(c_SR);
            % beat annotation
            kann=kann+1;
            annTime(kann,1) = round(t);
            annBeats(kann,1) = beatCodes(2);
            annRhythm{kann,1} = [];
        case 4
            % BT
            if ps ~= 4
                % decision between trigeminy and bigeminy
                x = rand;
                if x>BT_p(1)
                    % trigeminy
                    BT_type = 3;
                    ann = rhythmCodes(5);
                else
                    % bigeminy
                    BT_type = 2;
                    ann = rhythmCodes(4);
                end
                bt_flag = BT_type;
            end
            if bt_flag>1
                % insert normal beat
                rr(k) = rr_sr(c_SR);
                % check if rhythm annotation should change
                if rc~=4
                    % rhythm annotation - B or T
                    kann=kann+1;
                    annTime(kann) = t + round((rr(k)*1000)/2);
                    annBeats(kann) = beatCodes(end);
                    annRhythm{kann}= ann;
                    % rhythm check
                    rc = 4;
                    % BT minimum duration
                    BT_md =2;
                end
                % update time counter
                t = t + round(rr(k)*1000);
                % sp series is adjusted
                sp(c_SR:end) = sp(c_SR:end) + t-sp(c_SR);
                % update SR counter
                c_SR = c_SR + 1;
                % beat annotation
                kann=kann+1;
                annTime(kann) = round(t);
                annBeats(kann) = beatCodes(1);
                annRhythm{kann} = [];
                % update targets array
                targets_beats(k) = 1;
                % beat insterted
                bt_flag = bt_flag - 1;
                % the transition probabilites of the Markov chain change
                transM(4,1) = 0;
                transM(4,4) = 1;
            else
                % ventricular beat
                % prematurity factor
                beta_p = 0.75 + rand*0.1;
                % rr interval of the apb
                rr(k) = rr_sr(c_SR) * beta_p;
                % delay factor
                beta_f = 1.1 + rand*0.2;
                % apb compensatory pause
                rr_sr(c_SR) = rr_sr(c_SR) * beta_f;
                % update time counter
                t = t + round(rr(k)*1000);
                % sp series is adjusted
                sp(c_SR:end) = sp(c_SR:end) + t-sp(c_SR);
                % beat annotation
                kann=kann+1;
                annTime(kann,1) = round(t);
                annBeats(kann,1) = beatCodes(3);
                annRhythm{kann,1} = [];
                % update targets array
                targets_beats(k) = 4;
                % update the BT minimum duration variable
                BT_md = BT_md - 1;
                % check whether the minimum duration has expired
                if BT_md<=0
                    % minimum duration of the BT episode has expired
                    transM(4,1) = p_BT_SR;
                    transM(4,4) = 1 - p_BT_SR;
                end
                % the next beat should be normal if Bt remains active
                bt_flag = BT_type;
            end
        case 5
            % Isolated VPB
            vpbtype = randi([1,3]);
            switch vpbtype
                case 1 % VPBs with full compensatory pause
                    % prematurity factor
                    beta_p = 0.55 + rand*0.4;
                    % rr interval of the vpb
                    rrtemp = rr_sr(c_SR) * beta_p;
                    % vpb compensatory pause
                    rr_sr(c_SR) = (2*rr_sr(c_SR)) - rrtemp;
                case 2 % VPBs with non compensatory pause
                    % prematurity factor
                    beta_p = 0.55 + rand*0.4;
                    % rr interval of the apb
                    rrtemp = rr_sr(c_SR) * beta_p;
                case 3 % Interpolated VPBs
                    % prematurity factor
                    beta_p = 0.55 + rand*0.1;
                    % rr interval of the apb
                    rrtemp = rr_sr(c_SR) * beta_p;
                    % apb compensatory pause
                    rr_sr(c_SR) = rr_sr(c_SR) - rrtemp;
            end
            % insert ventricular beat
            rr(k) = rrtemp;
            % update time counter
            t = t + round(rr(k)*1000);
            % sp series is adjusted
            sp(c_SR:end) = sp(c_SR:end) + t-sp(c_SR);
            % beat annotation
            kann=kann+1;
            annTime(kann,1) = round(t);
            annBeats(kann,1) = beatCodes(3);
            annRhythm{kann,1} = [];
            % update targets array
            targets_beats(k) = 4;
        case 6
            % Isolated VPB during AF
            % prematurity factor
            beta_p = 0.55 + rand*0.4;
            % VPBs should not happen too close to other beats in AF
            rrtemp = rr_af(c_AF) * beta_p;
            while rrtemp < 250
                beta_p = 0.55 + rand*0.4;
                rrtemp = rr_af(c_AF) * beta_p;
            end
            % insert ventricular beat
            rr(k) = rrtemp;
            % update time counter
            t = t + round(rr(k)*1000);
            % afp series is adjusted
            afp(c_AF:end) = afp(c_AF:end) + t-afp(c_AF);
            % beat annotation
            kann=kann+1;
            annTime(kann,1) = round(t);
            annBeats(kann,1) = beatCodes(3);
            annRhythm{kann,1} = [];
            % update targets array
            targets_beats(k) = 4;
        case 7
            % isolated VPB during AT
            % prematurity factor
            beta_p = 0.55 + rand*0.4;
            % VPBs should not happen too close to other beats in AT
            rrtemp = rr_sr(c_SR-1) * beta_p;
            while rrtemp < 0.250
                beta_p = 0.55 + rand*0.4;
                rrtemp = rr_sr(c_SR-1) * beta_p;
            end
            % insert ventricular beat
            rr(k) = rrtemp;
            % update time counter
            t = t + round(rr(k)*1000);
            % sp series is adjusted
            sp(c_SR:end) = sp(c_SR:end) + t-sp(c_SR);
            % beat annotation
            kann=kann+1;
            annTime(kann,1) = round(t);
            annBeats(kann,1) = beatCodes(3);
            annRhythm{kann,1} = [];
            % update targets array
            targets_beats(k) = 4;
    end
    ps = s;
end
targets_beats = targets_beats(1:k);
rr = rr(1:k);
state_history = state_history(1:k);
annTime = annTime(1:kann);
annBeats = annBeats(1:kann);
annRhythm = annRhythm{1:kann};

%% Number and boundaries of PAF episodes
k = 1;
pafEpisodeLength = [];
for p = 1:length(state_history)-1
    if state_history(p) == 2
        if state_history(p+1) == 2
            k = k + 1;
            if p == (length(state_history)-1)
                pafEpisodeLength = [pafEpisodeLength k];
            end
        else
            pafEpisodeLength = [pafEpisodeLength k];
            k = 1;
        end
    end
end
pafBoundaries=0;
% % Find boundaries of each PAF episode
% if prod(state_history == 1)==1 % The entire signal is SR
%     pafBoundaries(1,1) = 0;
%     pafBoundaries(1,2) = 0;
% elseif prod(state_history == 2)==1 % The entire signal is AF
%     pafBoundaries(1,1) = 1;
%     pafBoundaries(1,2) = length(rr);
% else % The signal with PAF
%     diffTar = diff(state_history);
%     j = 1;
%     k = 1;
%     flag = 1;
%     for i = 1:length(diffTar)
%         if diffTar(i) == 1
%             pafBoundaries(j,1) = i + 1;
%             j = j + 1;
%             flag = 1;
%         end
%         if diffTar(i) == -1
%             pafBoundaries(k,2) = i;
%             k = k + 1;
%             flag = 2;
%         end
%         
%         if i == length(diffTar)
%             if flag == 1
%                 pafBoundaries(k,2) = i+1;
%             end
%         end
%     end
% end

end
function rr = simECG_generate_AF_intervals(lamba,nRR)
% [] = simECG_generate_AF_intervals() returns AF RR series modelled according
% to the paper by Corino et al. An atrioventricular node model for analysis
% of the ventricular response during atrial fibrillation. IEEE Transactions
% on Biomedical Engineering. 2011, 58(12), 3386-3395.
% 
% Licensed under GNU General Public License version 3:
% https://www.gnu.org/licenses/gpl-3.0.html 

if nRR < 70       
    time = 1;
else
    time = nRR/70;
end

prob_alpha = 0.6;
P = [0.15 0.25];
difftau = 0.2;
slope_beta = 10;

% tau1 non è costante ma dipendente dall'rr
% P(1)=slope
% P(2)=intercetta
tau1_in=polyval(P,1);
tau1=polyval(P,0:.001:3);
P2=[P(1) P(2)+difftau];
tau2=polyval(P2,0:.001:3);
tau2_in=polyval(P2,1);

Vt=-40;
Vr=-90;
dVdt=0;

aa=exprnd(1/lamba,1000000,1);% tolto 27 aprile: +50/1000; % minimo 50msec
t_a=cumsum(aa);
aa=aa(t_a<time*60); % 30 min const
atr_t=0;
% vtr_t=0;
avj_tmA=0;
vtr_tmA=0;
phase='phase4';
Vm=Vr;
Ts=0.001; 
t=0;
% avj_t=0;
nextAA=aa(1);
% vtr_nA=0;
i_R=0;
vtr_t=0;
vtr_nA=0;
avj_t=0;avj_nA=0;
R=zeros(length(aa),1);

nAA=0;
prob=zeros(length(aa),1);
prob1=round(length(aa)*prob_alpha); % prob_alpha è la prob che tau sia =tau1
prob(1:prob1)=1;
prob=prob(randperm(length(prob)));
prob_vera=zeros(size(prob));
tautau=zeros(size(aa));

tempo_x=0:.001:3;

beta=-slope_beta*tempo_x+1;
% % % più veloce beta(find(beta<0))=0;
beta(beta<0)=0;

% % beta=zeros(length(x));
% prob_beta=0 se passa, =1 se viene bloccato
prob2=round(length(aa)*beta); 
Nbeta=length(find(prob2>0));
prob_beta=zeros(length(aa),length(Nbeta));
% tic
for ii=1:Nbeta %length(x)
    prob_beta(1:prob2(ii),ii)=1;
    prob_beta(:,ii)=prob_beta(randperm(length(aa)),ii);
end
% toc   
blocked_beta=zeros(length(aa),1);
blocked_beta_t=zeros(length(tempo_x),1);
non_blocked_beta_t=zeros(length(tempo_x),1);
tempo_tot=zeros(length(aa),1);

% true=zeros(length(aa),1);
while nAA<length(aa)
    
%%%%%%%%%%%% function UpdateAtTs
%%%%%%%%%%%% At every msec this is what happens: Vm is linearly increased
%%%%%%%%%%%% by dVdt*Ts if in phase4; or the avj_t is increased

    t=t+Ts; % update time
    %i=i+1;disp([i vtr_nA])
    atr_t=atr_t+Ts;
    vtr_t=vtr_t+Ts;
    
    avj_tmA=avj_tmA+Ts;
    vtr_tmA=vtr_tmA+Ts;
    
    % update timers, do i need them??????? tmA for AVJ and vtr

    if strcmp(phase,'phase4') % otherwise it's in refractory period

        Vm=Vm+dVdt*Ts; 
        
    else avj_t=avj_t+Ts;
    end
    
%%%%%%%%%%%% function AnteHitAVJ
%%%%%%%%%%%% an AA arrives at the AVJ: Vm is increased of deltaV if in phase4 or tau (=RP) is prolonged

    if atr_t>=nextAA
        
%         t1=atr_t;
        atr_t=0;
        nAA=nAA+1; % increments AF counter
        if nAA<length(aa)
            nextAA=aa(nAA+1); % next AA
        end

        if strcmp(phase,'phase4') % otherwise it's in refractory period
%         phase=='phase4' 
            
            beats=find(R>0);
            if beats>0
                    %R(1)>0 %length(R)>=1
                tempo=t-R(beats(end))-tau;
                [m, pos]=min(abs(tempo_x-tempo));
                if pos<=Nbeta
                    if prob_beta(nAA,pos)==0
                        deltaV=(Vt-Vr)+1; % sempre soprasoglia!!
                        Vm=Vm+deltaV;
                        non_blocked_beta_t(pos)=non_blocked_beta_t(pos)+1;
                    else blocked_beta(nAA)=1;
                        blocked_beta_t(pos)=blocked_beta_t(pos)+1;
    %                     blocked_beta=blocked_beta+1;
                        tempo_tot(nAA)=tempo;
                    end 
                else % passano tutti!!
                    deltaV=(Vt-Vr)+1; % sempre soprasoglia!!
                    Vm=Vm+deltaV;
                    non_blocked_beta_t(pos)=non_blocked_beta_t(pos)+1;
                end
            else 
                deltaV=(Vt-Vr)+1; % sempre soprasoglia!!
                Vm=Vm+deltaV;
            end            
        end       
    end  
%%%%%%%%%%%% function VtrSense
%%%%%%%%%%%% there is a wave in the ventricle (vtr_nA>0) and the ventricle is not in refractory period (vtr_tmA>=AntDly)
 
    if vtr_nA>0 %& vtr_tmA>=AntDly

        i_R=i_R+1;
        R(i_R)=t;
        vtr_tmA=0;
        vtr_t=0;
        vtr_nA=0;
    end
      
%%%%%%%%%%%% function ActivateAVJ + StartAVJref if in phase4:
%%%%%%%%%%%% when Vm>=Vt AVJ is activated (avj_nA=avj_nA+1;) and then the
%%%%%%%%%%%% RP starts in AVJ

%%%%%%%%%%%% OR StartAVJph4 if in phase0
%%%%%%%%%%%% when avj_t>tau, i.e. the RP of AVJ finishes, phase4 again

    if strcmp(phase,'phase4') % otherwise it's in refractory period
        if Vm>=Vt       
%           %% ActivateAVJ 
            avj_nA=avj_nA+1; 
            vtr_nA=vtr_nA+1; % cohen non ha ritardo AV      
            phase='phase0';
            avj_t=0;
            
            if prob(nAA)==0 % va attraverso tau2
                
                beats=find(R>0);
                if beats>0
                    new_RR=t-R(beats(end));
                    [m xx]=min(abs(tempo_x-new_RR));
                    tau=tau2(xx);
                    prob_vera(nAA)=2;
                else
                    tau=tau2_in;
                    prob_vera(nAA)=2;
                end               
            else % va attraverso tau1
                beats=find(R>0);
                if beats>0
                    new_RR=t-R(beats(end));
                    [m, xx]=min(abs(tempo_x-new_RR));
                    tau=tau1(xx);
                    prob_vera(nAA)=1;
                else
                    tau=tau1_in;
                    prob_vera(nAA)=1;
                end
            end
            tautau(nAA)=tau;
        end 
        
    elseif exist('tau1','var')==1                  
        if avj_t>tau
            %% StartAVJph4
            phase='phase4';
            Vm=Vr;
        end

    end
   
end
R=R(R>0);
rr=diff(R);

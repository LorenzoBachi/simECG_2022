function multilead_P_waves = simECG_generate_multilead_P_waves(num)
%
% multilead_P_waves = simECG_generate_multilead_P_waves() returns multilead (15
% lead) simulated P waves. P waves are simulated using Hermite functions as
% a basis.
%
% num - number of P waves to be generated.
%
% Generated leads:
% multilead_P_waves(1,:) - I        multilead_P_waves(7,:) - V1           
% multilead_P_waves(2,:) - II       multilead_P_waves(8,:) - V2        
% multilead_P_waves(3,:) - III      multilead_P_waves(9,:) - V3        
% multilead_P_waves(4,:) - aVR      multilead_P_waves(10,:) - V4 
% multilead_P_waves(5,:) - aVL      multilead_P_waves(11,:) - V5 
% multilead_P_waves(6,:) - aVF      multilead_P_waves(12,:) - V6 
%
% multilead_P_waves(13,:) - X
% multilead_P_waves(14,:) - Y
% multilead_P_waves(15,:) - Z
%
multilead_P_waves = zeros(15,num,120);
p_frank = zeros(3,120);

% Guillem transform matrix  
D = [-0.266  0.027  0.065  0.131 0.203 0.220  0.370 -0.154;
      0.088 -0.088  0.003  0.042 0.047 0.067 -0.131  0.717;
     -0.319 -0.198 -0.167 -0.099 0.009 0.060  0.184  0.114];
% Compute inversion of Guillem's transformation matrix for the purpose to
% generate 12 leads
T=pinv(D'*D)*D';  

% % Time array (generates 150 samples)
% t = -3.8:0.05:3.7-0.03;

% Time array (generates 120 samples) %Cris 03/2022 
t = -3:0.05:3-0.01;

% Get parameters for P wave in Frank lead X
kx0 = simECG_random_number(0.04, 0.09);
kx1 = simECG_random_number(-0.02, 0.02);
kx2 = simECG_random_number(0, 0.02);
bx0 = simECG_random_number(0.8, 0.9); %Cris 03/2022 
bx1 = simECG_random_number(0.5, 0.9); %Cris 03/2022 
bx2 = simECG_random_number(0.5, 0.9); %Cris 03/2022 

% Get parameters for P wave in Frank lead Y
ky0 = simECG_random_number(0.05, 0.125); %Cris 03/2022 
ky1 = simECG_random_number(-0.03, 0.03);
ky2 = simECG_random_number(0, 0.03); 
by0 = simECG_random_number(0.8, 0.9); %Cris 03/2022 
by1 = simECG_random_number(0.5, 0.9); %Cris 03/2022 
by2 = simECG_random_number(0.5, 0.9);%Cris 03/2022 

% Get parameters for P wave in Frank lead Z
kz0 = simECG_random_number(-0.02, 0.02);  
kz1 = simECG_random_number(-0.02, -0.05);  
kz2 = simECG_random_number(-0.02, 0); 
bz0 = simECG_random_number(0.5, 0.9); %Cris 03/2022 
bz1 = simECG_random_number(0.7, 0.9); %Cris 03/2022 
bz2 = simECG_random_number(0.6, 0.9); %Cris 03/2022 

% Generate dominant frequency of P variability (between 1 and 3 Hz)
fk0 = simECG_random_number(1, 3);
fk1 = simECG_random_number(1, 3);
fk2 = simECG_random_number(1, 3);

fb0 = simECG_random_number(1, 3);
fb1 = simECG_random_number(1, 3);
fb2 = simECG_random_number(1, 3);

% Generate phase for each variability related frequency
if rand(1,1)> 0.5
    phaseX = 0;
else 
    phaseX = pi;
end

if rand(1,1)> 0.5
    phaseY = 0;
else 
    phaseY = pi;
end

if rand(1,1)> 0.5
    phaseZ = 0;
else 
    phaseZ = pi;
end

% Generate required number of P waves in Frank leads X, Y, Z
for n = 1:num
  
kx0m = simECG_variable_P_morphology(fk0, phaseX, kx0, 0.04, 0.09, n);
kx1m = simECG_variable_P_morphology(fk1, phaseX, kx1, -0.02, 0.02, n);
kx2m = simECG_variable_P_morphology(fk2, phaseX, kx2, 0, 0.02, n);

bx0m = simECG_variable_P_morphology(fb0, phaseX, bx0, 0.8, 1, n);
bx1m = simECG_variable_P_morphology(fb1, phaseX, bx1, 0.5, 1, n);
bx2m = simECG_variable_P_morphology(fb2, phaseX, bx2, 0.5, 1, n);
    
[p_frank(1,:), ~, ~, ~] =  simECG_generate_single_P_wave(t,kx0m,kx1m,kx2m,bx0m,bx1m,bx2m);
 
ky0m = simECG_variable_P_morphology(fk0, phaseY, ky0, 0.05, 0.15, n);
ky1m = simECG_variable_P_morphology(fk1, phaseY, ky1, -0.03, 0.03, n);
ky2m = simECG_variable_P_morphology(fk2, phaseY, ky2, 0, 0.03, n);

by0m = simECG_variable_P_morphology(fb0, phaseY, by0, 0.8, 1, n);
by1m = simECG_variable_P_morphology(fb1, phaseY, by1, 0.5, 1, n);
by2m = simECG_variable_P_morphology(fb2, phaseY, by2, 0.5, 1, n);
    
[p_frank(2,:), ~, ~, ~] =  simECG_generate_single_P_wave(t,ky0m,ky1m,ky2m,by0m,by1m,by2m);
 
kz0m = simECG_variable_P_morphology(fk0, phaseZ, kz0, -0.02, 0.02, n);
kz1m = simECG_variable_P_morphology(fk1, phaseZ, kz1, 0.02, 0.05, n);
kz2m = simECG_variable_P_morphology(fk2, phaseZ, kz2, -0.02, 0, n);

bz0m = simECG_variable_P_morphology(fb0, phaseZ, bz0, 0.5, 1, n);
bz1m = simECG_variable_P_morphology(fb1, phaseZ, bz1, 0.7, 1, n);
bz2m = simECG_variable_P_morphology(fb2, phaseZ, bz2, 0.6, 1, n);
    
[p_frank(3,:), ~, ~, ~] =  simECG_generate_single_P_wave(t,kz0m,kz1m,kz2m,bz0m,bz1m,bz2m);

% p_frank(1,:) = p_frank(1,:)*1.5;
% p_frank(2,:) = p_frank(2,:)*1.2;
% p_frank(3,:) = p_frank(3,:)*1.5;

p_sim = T*p_frank; 
 
 multilead_P_waves(1:8,n,:) = p_sim(1:8,:); % V1, V2, V3, V4, V5, V6, I, II, 
 multilead_P_waves(9,n,:) = p_sim(8,:) - p_sim(7,:);     % III
 multilead_P_waves(10,n,:) = -(p_sim(7,:) + p_sim(8,:))/2; % aVR
 multilead_P_waves(11,n,:) =  p_sim(7,:) - p_sim(8,:)/2;   % aVL
 multilead_P_waves(12,n,:) =  p_sim(8,:) - p_sim(7,:)/2;   % aVF
    
 multilead_P_waves(13,n,:) = p_frank(1,:); % X
 multilead_P_waves(14,n,:) = p_frank(2,:); % Y
 multilead_P_waves(15,n,:) = p_frank(3,:); % Z

end

p_waves_temp = multilead_P_waves;
multilead_P_waves(1:6,:,:) = multilead_P_waves(7:12,:,:);
multilead_P_waves(7:12,:,:) = p_waves_temp(1:6,:,:);

end
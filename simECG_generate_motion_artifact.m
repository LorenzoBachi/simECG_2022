function [simuMA] = simECG_generate_motion_artifact(ecgLength, simECGdata, noiseRMS)
% [] = simECG_generate_motion_artifact() returns a 15-lead simulated motion
% artifact.
% 
% Input arguments:
% ecgLength - length of signal, in samples.
% simECGdata - struct of ECG simulation parameters defined in the main
% script.
% noiseRMS - desired noise level amplitude.
% 
% Output arguments:
% simuMA - simulated 15-lead motion artifact.
% 
% Author: Hesam Halvaei, Lund University
% 
% Licensed under GNU General Public License version 3:
% https://www.gnu.org/licenses/gpl-3.0.html

sigmay = 0.03*noiseRMS; 

% if simECGdata.MA_Flag %Thumb-ECG
%     L = 1;
% else %8 indipendent leads or 1-lead
%     L = 8;
% end
L = 8;
for Li = 1:L
    bernogauss = simECG_Bernoulli_Gaussian(ecgLength, simECGdata.MA_Prob, sigmay);
    bernogauss_conv(Li,:) = simECG_Bernoulli_Gaussian_convolution(ecgLength, bernogauss, simECGdata.MA_Flag)';
end

simuMA = bernogauss_conv; 

% if ~simECGdata.MA_Flag % Transform to the 15 leads
    simuMA_8 = simuMA; %leadcalc(simuMA,'stan');% V1,V2,V3,V4,V5,V6,I,II
    simuMA_12 = leadcalc(simuMA_8,'extr');% V1,V2,V3,V4,V5,V6,aVL,I,-aVR,II,aVF,III
    simuMA_VCG = leadcalc(simuMA_8,'synt');
    
    simuMA_15 = vertcat(simuMA_8(7,:),simuMA_8(8,:),simuMA_12(12,:),...
        -simuMA_12(9,:),-simuMA_12(7,:),simuMA_12(11,:),simuMA_8(1:6,:),simuMA_VCG);
    
    simuMA = simuMA_15;
% end

end


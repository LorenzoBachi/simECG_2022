function [simNoise] = simECG_Bernoulli_Gaussian_convolution(len, bernogauss, MA_Flag);
% func_BernoulliGaussianConvolution 
% Convolution of the a realizaiton of bernoulli gaussian process and filtering
% INPUTS
% bernogauss: a bernoulli gaussian realization
% MA_Flag: 0 Holter /1 ThumbECG

% OUTPUT 
% simNoise: simulated noise
%--------
% Hesam Halvaei, Lund University
%--------

ar_coeffs = [1.0000   -0.2758    0.3098   -0.1025    0.2495]; % Coeffients for MA
% OR INSTEAD 
% load("DATA_AR_MN_Dictionary.mat");
%ar_coeffs = AR_MN(randi([1, size(AR_MN, 1)]), :, randi([1, size(AR_MN, 3)]));
bernogauss_conv = zeros(size(bernogauss));
N = 400;
inds = find(bernogauss ~= 0); % finds the non-zero indices from the process
if ~isempty(inds)
for q = 1:numel(inds)
    h = zeros(N, 1);
    a1 = unifrnd(.7,.9);
    a2 = unifrnd(.7,.9);
    K = randi([125, 375]);
    n = 1:K-1; 
    h(n) = a1.^(-n);
    n = K:N;
    h(n) = a1.^(-(K)) * (a2.^(n-K));
    h = h/(max(abs(h)));
    spike = zeros(size(bernogauss));
    spike(inds(q)) = bernogauss(inds(q));
    spike_conv = conv(spike, h, 'same'); 
    bernogauss_conv = bernogauss_conv + spike_conv;
end

ar_input = decimate(bernogauss_conv, 5);
ar_model = idpoly(ar_coeffs, 'integratenoise', MA_Flag);
simNoise_decimated= sim(ar_model, ar_input);
if MA_Flag
    [b, a] = butter(4, [.677]/100, 'high');
    simNoise_decimated = filtfilt(b, a, simNoise_decimated);
end
simNoise = interp(simNoise_decimated, 5);
simNoise = simNoise(1:len);


    
end
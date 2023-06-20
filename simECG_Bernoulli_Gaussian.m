function [bernogauss] = simECG_Bernoulli_Gaussian(n, prob, sigma)
% Realization of a Bernoulli-Gaussian process
%INPUT:
% n: length of signal
% prob: probability of success in Bernoulli Gaussian Process
% sigma Noise level
% standard deviation of guassian noise
%OUTPUT
% bernogauss: Bernoulli-Gaussian Process
%--------
% Hesam Halvaei, Lund University
%--------
N = ones(n, 1);
p = N.*prob;
berno = binornd(N, p);
gauss= normrnd(0, sigma, [n, 1]);
bernogauss = gauss.*berno;
end

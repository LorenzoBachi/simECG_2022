function [bernogauss] = simECG_Bernoulli_Gaussian(n, prob, sigma)
% Realization of a Bernoulli-Gaussian process
N = ones(n, 1);
p = N.*prob;
berno = binornd(N, p);
gauss= normrnd(0, sigma, [n, 1]);
bernogauss = gauss.*berno;
end

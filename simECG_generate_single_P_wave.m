function [P, fi0, fi1, fi2] = simECG_generate_single_P_wave(t, k0,k1,k2,b0,b1,b2)
%
% [] = simECG_generate_single_P_wave() returns a single P wave.
%
% Hermite functions
    fi0 = k0*(1/(sqrt(b0*sqrt(pi))))*exp((-t.^2)/(2*b0.^2));
    fi1 = -1*k1*((sqrt(2))/(sqrt(b1*(sqrt(pi)))))*(t/b1).*exp((-t.^2)/(2*(b1.^2)));
    fi2 = k2*(1/(sqrt(2*b2*sqrt(pi))))*((2*((t.^2)/(b2.^2)))-1).*exp((-t.^2)/(2*(b2.^2)));
    % Resulting P wave
    P = fi0 + fi1 + fi2;
end
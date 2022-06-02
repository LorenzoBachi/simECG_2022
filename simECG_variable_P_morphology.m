function par = simECG_variable_P_morphology(f, phase, par, rangeLow, rangeHigh, n)
    % 
    % par = simECG_variable_P_morphology returns parameter values that are used for
    % altering simulated P wave morphology in time.
    %
    Td = 0.05;
    par = par + (2*par*(rand(1,1)-0.5)/10) + ((rangeHigh-rangeLow)*(1+sin(2*pi*f*Td*n+phase))/2+rangeLow)/10;
end
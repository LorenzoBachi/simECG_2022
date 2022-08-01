function [vpb] = simECG_generate_XYZ_ventricular_VA(fitCoeff,J)
%
% [] = simECG_generate_ventricular_XYZ_VA() returns simulated QRST for
% ventricular beats, generated from the Hermite and logistic function
% coefficients resulting from fitting the VPB from Alzaraz database. The
% fitting was performed according to Sörnmo et al. 1981 and Bock et al.,
% 2021. Generated VPBs are 540 samples long at 500 Hz, in the Frank lead
% space.

%get size of fitCoeff
dims = size(fitCoeff);
%number of VPBs in patch
N=dims(1);
%preallocate result matrix
vpb = zeros(3,N,540);
%numer of Hermite functions used for QRS complex and T wave
N_QRS = 6;
N_T = 4;

for k=1:N
    for l=1:3
        %global time vector
        t_all = 1:700;
        %QRS time vector
        t_QRS = 1:J(k);
        %T wave time vector
        t_T = (J(k)+1):700;
        %preallocate QRS portion of the beat
        y_QRS_hat = zeros(1,length(t_all));
        %fetch QRS complex Hermite amplitude coefficients
        c_QRS = squeeze(fitCoeff(k,l,1:N_QRS));
        %fetch QRS complex Hermite function width and displacement
        %coefficients
        temp = squeeze(fitCoeff(k,l,N_QRS+(1:2)));
        sigma_QRS = temp(1); tau_QRS = temp(2);
        %fetch QRS complex logistic function amplitude, witdh and
        %displacement coefficients
        temp = squeeze(fitCoeff(k,l,N_QRS+(3:5)));
        %rebuild QRS complex
        c_s_QRS = temp(1); sigma_s_QRS = temp(2); tau_s_QRS = temp(3);
        for n=1:N_QRS
            hf = simECG_hermite_function(n-1,sigma_QRS,t_all-tau_QRS);
            y_QRS_hat = y_QRS_hat + c_QRS(n)*hf;
        end
        %preallocate T wave portion of the T wave
        y_T_hat = zeros(1,length(t_all));
        %fetch T wave Hermite amplitude coefficients
        c_T = squeeze(fitCoeff(k,l,N_QRS+5+(1:N_T)));
        %fetch T wave Hermite function width and displacement coefficients
        temp = squeeze(fitCoeff(k,l,N_QRS+N_T+(6:7)));
        sigma_T = temp(1); tau_T = temp(2);
        %fetch T wave logistic function amplitude, witdh and displacement
        %coefficients
        temp = squeeze(fitCoeff(k,l,N_QRS+N_T+(8:10)));
        c_s_T = temp(1); sigma_s_T = temp(2); tau_s_T = temp(3);
        %rebuild T wave
        for n=1:N_T
            hf = simECG_hermite_function(n-1,sigma_T,t_all-tau_T);
            y_T_hat = y_T_hat + c_T(n)*hf;
        end
        %combine QRS complex and T wave Hermite functions
        y_all_hat = y_QRS_hat + y_T_hat;
        %check whether logistic function is needed
        if (c_s_T~=0)&&(c_s_QRS~=0)
            %add logistic function
            e_QRS_hat = c_s_QRS * simECG_logistic_function( -(t_QRS - tau_s_QRS )/(sigma_s_QRS) );
            e_T_hat = c_s_T * simECG_logistic_function( ( t_T - tau_s_T )/(sigma_s_T) );
        else
            %no logistic function needed
            e_QRS_hat = zeros(1,length(t_QRS));
            e_T_hat = zeros(1,length(t_T));
        end
        %combine logistic functions
        e_all_hat = [e_QRS_hat, e_T_hat];
        %combine Hermite and logistic functions
        x_hat = y_all_hat + e_all_hat;
        %the last padding section is discarded
        vpb(l,k,:) = x_hat(1:540);
    end
end

end

function [f] = simECG_hermite_function(n,sigma,t)
%SIMECG_HERMITE_FUNCTION computes the n-th Hermite function. 

f = 1/sqrt(sigma*(2^n)*factorial(n)*sqrt(pi)) * ( exp((-(t.^2))/(2*(sigma^2))) .* simECG_hermite_polynomial(n,t/sigma) ) ;

end

function [p] = simECG_hermite_polynomial(n,x)
%SIMECG_HERMITE_POLYNOMIAL returns the n-th phycisit's Hermite polynomial calculated
%over x. Only first eight polynomials are considered.

switch n
    case 0
        p = ones(size(x));
    case 1
        p = 2*x;
    case 2
        p = 4*x.^2 - 2;
    case 3
        p = 8*x.^3 - 12*x;
    case 4
        p = 16*x.^4 - 48*x.^2 + 12;
    case 5
        p = 32*x.^5 - 160*x.^3 + 120*x;
    case 6
        p = 64*x.^6 - 480*x.^4 + 720*x.^2 - 120;
    case 7
        p = 128*x.^7 - 1344*x.^5 + 3360*x.^3 - 1680*x;
    otherwise
        p=zeros(size(x));
end

end

function [f] = simECG_logistic_function(t)
%SIMECG_LOGISTIC FUNCTION computes the logistic function.

f = 1 ./ (1 + exp(t));

end




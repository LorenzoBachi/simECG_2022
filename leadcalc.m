function y = leadcalc(x,type)

%Function
% leadcalc
% 
%Purpose
% various lead calculations
%
%Synopsis
% y = leadcalc(x,type)
% y = leadcalc(x)
%
%Description
% The function calculates extremity/limb leads for the standard 12-lead ECG
% or synthesized VCG leads based on either an (8xn) or (9xn) matrix
% containing V1,V2,V3,V4,V5,V6,I,II and, eventually, III.
%
% The desired lead configuration is specified by the character string 'type'.
% The string 'extr' produces the missing extremity/limb leads and the result-
% ing output matrix contains the leads V1,V2,V3,V4,V5,V6,aVL,I,-aVR,II,aVF,III
% (stored according to the panoramic display format).
%
% The string 'synt' results in a synthesized VCG (leads X, Y and Z) using
% the inverse Dower matrix for lead synthesis.
%
% The string 'syntpd' result in a synthesized VCG (leads X, Y and Z) using
% the inverse Dower matrix for pedriatic ECG's lead synthesis. cf. Edenbrandt et al.:
% VCGs Derived from 12-Lead ECGs; Pediatr Cardiol 15:21-26, 1994.
%
% The string 'hcuzst' result in a synthesized VCG (leads X, Y and Z) using
% the inverse Dower matrix when V2 is substituted by RV4 in stress test of
% the University Clinic Hospital of Zaragoza in stress tests.
%
% The string 'stan' produces the standard leads and the resulting output
% matrix contains the leads V1,V2,V3,V4,V5,V6,I,II and III.  The input
% matrix must contain the frank leads X,Y,Z.
%
% The string 'levk' results in a synthesized VCG (leads X, Y and Z) using
% Levkov's T1 matrix for lead synthesis, cf. Med. & Biol. Eng. & Comp., March
% 1987, 25, 155-164.
%
% The string 'kors' results in a synthesized VCG (leads X, Y and Z) using
% Kors Regression method, cf. Reconstruction of the Frank vectorcardiogram from
% standard electrocardiographic leads: diagnostic comparison of different methods
% European Heart Journal, 1990, 11, 1083-1092 
%
% The default string is 'extr'.
%
% Copyright (c), Leif Sornmo, Lund University, 98-01-26
% Last modification: Juan Bolea, march 2009
%----------------------------------------------------------------------------


%----------------------------------------------------------------------------
if size(x,1) > size(x,2)
	disp('WARNING: input ECG data is transposed')
	x = x';
end
no_leads = size(x,1);

if nargin == 1        
    type = 'extr';
else
    if ~strcmp(type,'stan')
        if no_leads < 8 || no_leads > 9
            error('only 8 or 9 leads accepted');
        end
    else
        if no_leads ~= 3
            error('only 3 leads accepted');
        end
    end
end            

%----------------------------------------------------------------------------


%----------------------------------------------------------------------------
if ~isempty(strfind(type,'extr'))
% extremity leads
    if strcmp(type,'extr')  % Label obsolete
        type = 'extr3';
    end
	y = zeros(12,size(x,2));
    y(1:6,:) = x(1:6,:);				% V1, V2, V3, V4, V5, V6
    if str2double(type(5)) == 1 %extr1 I lead is lost so lead 7 is II and lead 8 is III
        y(10,:) = x(7,:);				% II
        y(8,:) = x(7,:)-x(8,:);  % I = II - III
        y(12,:) = x(8,:); % III
    elseif str2double(type(5)) == 2 % extr2 II lead is lost so lead 7 is I and lead 8 is III
        y(8,:) = x(7,:);					% I
        y(10,:) = x(8,:) + x(7,:);
        y(12,:) = x(8,:); % III
    elseif str2double(type(5)) == 3 % extr3 III lead is lost so lead 7 is I and lead 8 is II
        y(8,:) = x(7,:);	
        y(10,:) = x(8,:); %RUTE% II
        if no_leads == 8
            y(12,:) = x(8,:) - x(7,:);		% III = II-I
        else
            y(12,:) = x(9,:);
        end
    end
	y(7,:) = x(7,:) - 0.5 * y(10,:);		% aVL = I-II*0.5	
	y(9,:) = (y(8,:)+y(10,:))/2;			% -aVR = (I+II) / 2 	
	y(11,:) = y(10,:) - 0.5 * y(8,:);	% aVF = II-I/2
	
%----------------------------------------------------------------------------
	

%----------------------------------------------------------------------------
% synthesised VCG with inverse Dower matrix
elseif strcmp(type,'synt')
	t = [-.172 -.074 .122 .231 .239 .194 .156  -.010;...
		  .057 -.019 -.106 -.022 .041 .048 -.227 .887;...
	     -.229 -.310 -.246 -.063 .055 .108 .022 .102]; 
	y = t * x(1:8,:);

%----------------------------------------------------------------------------


%--------------------------------------------------------------------------
% synthesised VCG with inverse Dower matrix for pediatric ECG's
% Based on the image surface data of Frank, the coefficients of lead V3 in
% the Dower matrix (0.882; 0.098; -1.277) were replaced with corresponding
% coefficients for lead V_4R (-0.537; 0.096; -0.272). 
elseif strcmp(type,'syntpd')
    t = [-.122 .009 -.128 .275 .251 .185 .160 -.013;...
        .019 -.087 .073 -.065 .025 .051 -.235 .891;...
        -.278 -.439 -.072 -.189 -.016 .084 -.023 .128];
    y = t * x(1:8,:);

%----------------------------------------------------------------------------


%--------------------------------------------------------------------------
% synthesised VCG with inverse Dower matrix when V2 is substituted by RV4
% (HCUZ stress test recordings)
% Raquel Bailon Luesma, Feb2003
elseif strcmp(type,'hcuzst')

	t = [-.182 -.137 .088 .207 .226 .189 .147 -.002;...
		  .032 .066 -.121 -.019 .053 .063 -.217 .882;...
          -.371 -.132 -.415 -.126 .066 .160 .038 .109]; 
	y = t * x(1:8,:);

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% stardard leads with Dower transformation matrix
elseif strcmp(type,'stan')
    
    D = [-0.515  0.157   -0.917;
        0.044  0.164   -1.387;
        0.882  0.098   -1.277;
        1.213  0.127   -0.601;
        1.125  0.127   -0.086;
        0.831  0.076    0.230;
        0.632 -0.235    0.059;
        0.235  1.066   -0.132];
    y = D * x;

%--------------------------------------------------------------------------


%----------------------------------------------------------------------------
% synthesised VCG with Levkov's T1 matrix
elseif strcmp(type,'levk')
	t = [0.20 -0.56 -0.11 -0.02  0.09  0.11  0.18 0.27; ...
		-0.01 -0.91 -0.01 -0.01 -0.02 -0.02 -0.04 0; ...
		 0.17  0.31 -0.34 -0.20 -0.14 -0.06 -0.25 0.35];

	RA = -x(8,:); %RA = -II;
	if size(x,1) == 9
		LA = -x(9,:);  %LA = -III;
	else
		LA = x(8,:) - x(7,:);  %LA = -III;
	end
	C = x(1:6,:) + kron(ones(6,1),(RA+LA) ./ 3);
	y = t * [LA; RA; C];
    
elseif strcmp(type,'kors')
	t = [-0.13 0.05 -0.01 0.14 0.06 0.54 0.38 -0.07; ...
		0.06 -0.02 -0.05 0.06 -0.17 0.13 -0.07 0.93; ...
		 -0.43 -0.06 -0.14 -0.20 -0.11 0.31 0.11 -0.23];

	RA = -x(8,:); %RA = -II;
	if size(x,1) == 9
		LA = -x(9,:);  %LA = -III;
	else
		LA = x(8,:) - x(7,:);  %LA = -III;
	end
	C = x(1:6,:) + kron(ones(6,1),(RA+LA) ./ 3);
	y = t * [LA; RA; C];
    
else
	error('unknown type')
end
%----------------------------------------------------------------------------

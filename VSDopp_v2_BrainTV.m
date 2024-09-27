function velocityEstimate = VSDopp_v2_BrainTV(PRF,avgIQ_data,nlambda,GATE)
%% FCN: _[] = VSDopp_v2_BrainTV()
%
% *PURPOSE:*    This routine calculate brain tissue displacement based on
%               the SDopp autocorrelation velocity estimator (See Hoeks et al 1994)
%               - One time lag is T
%               - Subsample volume of 3 (1 burst covers 6 GATE (3 sample volumes)
%               - Equation 7a is the in-phase quadrature component
%               - Equation 7b is the shifted quadrature component
%               - Equation 9 is the coefficients for the complex autocorrelation function
%               - PL is the Package Length or no. of pulses over which the mean
%                 velocity is estimated
%               - ND is the number of IQ sample volumes to estimate velocity over 
%               - T is the period
%               - lambda is the wavelength of the ultrasound
%               
%
% *INPUT:*      PRF - Pulse Repetition Frequency
%               avgIQ_data - IQ data averaged over a user defined number of
%               samples
%               nlambda - lambda/4 * 1000 (Nyquist m)
%               GATE - no. of gates selected by the user to analyse (34)
%
% *OUTPUT:*     velocityEstimate - velocity in microns/s estimated from SDopp
%               autocorrelator velocity estimator
%                         
%
% *AUTHOR(S):*
%
% Algorithm Implementation by Dr. C. Banahan 24/06/2016
% Updated by BEng Student Jakub Markiewicz 1/04/2024
%--------------------------------------------------------------------------
%%

% T = 1 / PRF; test
T = PRF; %?

% For the first time lag we average over 4 IQ sample volumes and over 2
% pulse lengths. Subsequent times are averaged over 3 IQ sample volumes and
% over 1 pulse length

% Select package length and no. of subsamples
PL = 2;
ND = 4;

% Set up Sums 
% These are used for looping (=Sum symbol parameters in algorithm)
b = 2 - PL:1:0; %b = [0]
a = 2 - ND:1:0; %a = [-2,-1,0]

% Initialize matrix for storing phase value from Algorithm output
% Size must fit pair of IQ signals
phaseIQ = zeros(length(avgIQ_data)-2,2);

% Initialize matrix for storing phase values.
% Size must fit pair of IQ signals and allow for looping (== Sum in algorithm)
% For 4 subsamples, there are 3 columns.
phaseIQ_Matrix = zeros(length(avgIQ_data)-2,length(a));

% Initialize matrix for storing the final phase values (Loops/Sum finished)
% Size must fit pair of IQ signals at all depths/gates (no. of gates=34)
% For all depths, there are 34 columns in the final matrix
phaseDifference_Matrix = zeros(length(avgIQ_data)-2,GATE);

% Algorithm must execute on all IQ pairs
% Iterate over first gate via index by pulse no. (n) and depth (m) of IQ signal
for n = 1:length(avgIQ_data)-2 
    m = 1; 

    % Must average over ND subsamples
    % Introduce column index variable and iterate until all elements have a value
    for k = 1:length(a)
        
        % Initiliaze sum for averaging
        phaseIQ_Matrix(n,k) = 0;

        % Introduce row index variable and iterate until all elements have a
        % value
        for j = 1:length(b)

               % As j increases in each iteration, compute phase values
               % using SDopp/Autocorrelation algorithm

               if n==1 && b(j) == 0
                   % Algorithm operation
                   phaseIQ(j,k) = avgIQ_data(n+2-b(j),...
                       m-a(k)).*conj(avgIQ_data(n+1-b(j),m-a(k)));
                   
                   % Update phaseIQ_Matrix with result of algorithm output
                   phaseIQ_Matrix(n,k) = phaseIQ_Matrix(n,k)+ phaseIQ(j,k);
              
               elseif n==1 && b(j) ==1
                   % Algorithm operation
                   phaseIQ(j,k) = avgIQ_data(n-1+b(j),...
                       m-a(k)).*conj(avgIQ_data(n-1 + b(j)+1,m-a(k)));
                  
                   % Update phaseIQ_Matrix with result of algorithm output
                   phaseIQ_Matrix(n,k) = phaseIQ_Matrix(n,k) + phaseIQ(j,k);
               
               else
                   % Algorithm operation
                   phaseIQ(j,k) = avgIQ_data(n-b(j)+1,...
                       m-a(k)).*conj(avgIQ_data(n-b(j),m-a(k)));
                   
                   % Update phaseIQ_Matrix with result of algorithm output
                   phaseIQ_Matrix(n,k) = phaseIQ_Matrix(n,k)+ phaseIQ(j,k);
               end
        end
        % Once all iterations are complete (=Sum is complete), add all
        % phase results to final IQ Data x Gate matrix
        phaseDifference_Matrix(n,m) = phaseDifference_Matrix(n,m) + phaseIQ_Matrix(n,k); 
    end
end

% Select package length and no. of subsamples for subsequent time lags
PL = 1; 
ND = 3;

% Change indices so time 0 corresponds to index 1
% Prevents index errors
if PL == 1
    b = [1,0];
else
    b = 2 - PL:1:0;
end

% If there are no IQ subsamples, a = 0 and SDopp reduces to Autocorrelation
if ND == 1 
    a = 0; 
else
    a = 2 - ND:1:0;
end    

% Must iterate for all gates/depths of IQ signal
% Set up parameters based on subsample approach
if ND > 2
    l=GATE-1;
else
    l=GATE;
end

% Operation same as before, but now averaging over all subsamples at all
% depths
for m = 2:l 
    for n=1:length(avgIQ_data)-2 
        
        %Column index variable as before
        for k=1:length(a)
            phaseIQ_Matrix(n,k) = 0;
           
            %Row index variable as before
            for j=1:length(b)
               
                % using SDopp/Autocorrelation algorithm

                if n==1 && b(j) == 0
                   % Algorithm operation as before
                   phaseIQ(j,k) = avgIQ_data(n+2-b(j),...
                       m-a(k)).*conj(avgIQ_data(n+1-b(j),m-a(k)));
                   
                   phaseIQ_Matrix(n,k) = phaseIQ_Matrix(n,k)+ phaseIQ(j,k); 
              
                elseif n==1 && b(j) ==1
                   % Algorithm operation as before
                   phaseIQ(j,k) = avgIQ_data(n-1+b(j),...
                       m-a(k)).*conj(avgIQ_data(n-1 + b(j)+1,m-a(k)));
                   
                   phaseIQ_Matrix(n,k) = phaseIQ_Matrix(n,k) + phaseIQ(j,k);
               
                else
                   % Algorithm operation as before
                   phaseIQ(j,k) = avgIQ_data(n-b(j)+1,...
                       m-a(k)).*conj(avgIQ_data(n-b(j),m-a(k)));
                   
                   phaseIQ_Matrix(n,k) = phaseIQ_Matrix(n,k)+ phaseIQ(j,k);
               end
            end
            % Final matrix now holds all phase information for all
            % subsamples and depths of IQ signal
            phaseDifference_Matrix(n,m) = phaseDifference_Matrix(n,m) + phaseIQ_Matrix(n,k);
            
        end

    end
end

% Using Doppler equation, calculate velocity
% Use angle() for arg()
velocityEstimate = (nlambda/(2*pi*T)).*unwrap(angle(phaseDifference_Matrix)); 

%V_SDopp = (nlambda/(pi)).*unwrap(angle(phaseDifference_Matrix));  %TEST



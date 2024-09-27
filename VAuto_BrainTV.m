function autocorrelationVelocity = VAuto_BrainTV(PRF, avgIQ_data, nlambda, GATE)
%% FCN: [] = VAuto_BrainTV()
%
% PURPOSE: Verify the implementation of the SDopp algorithm by introducing
%          Autocorrelation and comparing results
%
%          Method will remain the same to have comparable results:
%          - One time lag is T
%          - Subsample volume of 3
%          - PL is the Package Length or no. of pulses over which the mean
%            velocity is estimated
%          - T is the period
%          - lambda is the wavelength of ultrasound
%
%
% INPUT:   PRF - Pulse repetition frequency
%          avgIQ_data - IQ data averged over a user defined number of
%          samples
%          nlambda - lambda/4 * 1000 (Nyquist m)
%          GATE - no. of gates selected by the user to analyse (34)
%
% OUTPUT:  autocorrelationVelocity - velocity in micron/s estimated from 
%          autocorrelation velocity estimator
%
% AUTHOR(S):
%
% Ideas borrowed from SDopp implementation by Dr. C. Banahan
% Implementation by BEng Student Jakub Markiewicz 1/04/2024
%--------------------------------------------------------------------------

%For first time lag average over 2 pulse lengths. Subsequent times are
%averaged over 1 pulse length. (Same as SDopp implementation).

T = PRF; %?
PL = 2;

% Initialize velocity matrix
autocorrelationVelocity = zeros(length(avgIQ_data)-2, GATE);

% Iterate over first gate
for m = 1
    % Iterate over IQ data points
    for n = 1:length(avgIQ_data)-2
        % Initialize sum for averaging
        autocorrelationSum = 0;
        
        % Average autocorrelation over PL pulse lengths
        for j = 0:PL-1
            % Calculate autocorrelation
            autocorrelation = avgIQ_data(n+2-j,m).*conj(avgIQ_data(n+1-j,m));
            
            % Add to sum
            autocorrelationSum = autocorrelationSum + autocorrelation;
        end
        
        % Compute phase difference
        phaseDifference = angle(autocorrelationSum);
        
        % Apply Kasai autocorrelation velocity estimator formula
        autocorrelationVelocity(n,m) = (nlambda/(2*pi*T)) * phaseDifference;
    end
end

PL = 1;

% Iterate over all other gates
for m = 2:GATE
    % Iterate over IQ data points
    for n = 1:length(avgIQ_data)-2
        % Initialize sum for averaging
        autocorrelationSum = 0;
        
        % Average autocorrelation over PL pulse lengths
        for j = 0:PL-1
            % Calculate autocorrelation
            autocorrelation = avgIQ_data(n+2-j,m).*conj(avgIQ_data(n+1-j,m));
            
            % Add to sum
            autocorrelationSum = autocorrelationSum + autocorrelation;
        end
        
        % Compute phase difference
        phaseDifference = angle(autocorrelationSum);
        
        % Apply Kasai autocorrelation velocity estimator formula
        autocorrelationVelocity(n,m) = (nlambda/(2*pi*T)) * phaseDifference;
    end
end

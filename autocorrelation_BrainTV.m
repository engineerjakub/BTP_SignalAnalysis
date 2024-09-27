function velocityEstimate = VSDopp_v2_BrainTV(T,avgIQ_data,nlambda,mm)
%% FCN: _[] = VSDopp_BrainTV()
%
% *PURPOSE:*    This routine calculate brain tissue displacement based on
%               the SDopp autocorrelation velocity estimator (See Hoeks et al 1994)
%               - One time lag is T
%               - Subsample volume of 3 (1 burst covers 6 mm (3 sample volumes)
%               - Equation 7a is the in-phase quadrature component
%               - Equation 7b is the shifted quadrature component
%               - Equation 9 is the coefficients for the complex autocorrelation function
%               - PL is the Packet Length or no. of pulses over which the mean
%                 velocity is estimated
%               - ND is the number of IQ sample volumes to estimate velocity over 
%               - T is the period
%               - lambda is the wavelength of the ultrasound
%               
%
% *INPUT:*      T - period
%               avgIQ_data - IQ data averaged over a user defined number of
%               samples
%               nlambda - lambda/4 * 1000 (Nyquist m)
%               mm - no. of gates selected by the user to analyse
%
% *OUTPUT:*     V_SDopp- displacement in microns estimated from SDopp
%               autocorrelator velocity estimator
%                         
%
% *AUTHOR(S):*
%
%               Dr. C. Banahan updated 24/06/2016
%               Jakub Markiewicz updated 20/04/2024
%--------------------------------------------------------------------------
%%

%For the first time lag we average over 4 IQ sample volumes and over 2
%pulse lengths. Subsequent times are averaged over 3 IQ sample volumes and
%over 1 pulse length.
%%

%Set PL and ND parameters
PL = 2;
ND = 4;

%Set up sums 
b = 2 - PL:1:0; %b = [0]
a = 2 - ND:1:0; %a = [-2,-1,0]

%Initiate zeros matrix array for storing values.
%Must have the size of average data -> samples * channels.
%There are 34 samples of IQ signals for each channel.
phaseIQ_Matrix = zeros(length(avgIQ_data)-2,length(a));
disp(phaseIQ_Matrix)

%Matrix with length of avgIQ_data cell rows and 2 columns.
phaseIQ = zeros(length(avgIQ_data)-2,2);

%Matrix with length of avgIQ_data cell rows and mm (no. of gates=34) columns.
phaseDifference_Matrix = zeros(length(avgIQ_data)-2,mm);

%for loop is on until n = length(avgIQ_data)-2
%e.g. if length is 4, then n = 1, then n = 2. (4-2 = 2)
for n = 1:length(avgIQ_data)-2 %pulse number
    m = 1; %depth

    %k = 1, then k = 2, k = 3. Loops three times.
    for k = 1:length(a)
        %Index n,k in matrix VD and set it to 0.
        phaseIQ_Matrix(n,k) = 0;

        %Only loops once, length of b is 1.
        for j = 1:length(b)

            %This if statement is for the first iteration.
               if n==1 && b(j) == 0

                   %Index j,k in matrix v and set it to output of the
                   %operation below.
                   phaseIQ(j,k) = avgIQ_data(n-b(j)+2,...
                       m-a(k)).*conj(avgIQ_data(n-b(j)+1,m-a(k)));
                   
                   phaseIQ_Matrix(n,k) = phaseIQ_Matrix(n,k)+ phaseIQ(j,k); 
               
               elseif n==1 && b(j)==1
                   phaseIQ(j,k) = avgIQ_data(n-b(j)-1,...
                       m-a(k)).*conj(avgIQ_data(n+b(j),m-a(k)));
                   
                   phaseIQ_Matrix(n,k) = phaseIQ_Matrix(n,k) + phaseIQ(j,k);
               
               else
                   phaseIQ(j,k) = avgIQ_data(n-b(j)+1,...
                       m-a(k)).*conj(avgIQ_data(n-b(j),m-a(k)));
                   
                   phaseIQ_Matrix(n,k) = phaseIQ_Matrix(n,k)+ phaseIQ(j,k);
               end
        end
        phaseDifference_Matrix(n,m) = phaseDifference_Matrix(n,m) + phaseIQ_Matrix(n,k); 
    end
end
%%
PL = 1; 
ND = 3;
%need to change indices so time 0 corresponds to index 1
if PL==1
    b = [1,0];
else
    b = 2 - PL:1:0;
end
if ND == 1 %or 2 a=0
    a = 0; %not averaging over different IQ sample volumes (reduces to autocorrelator estimator)
else
    a = 2 - ND:1:0;
end    
    
if ND > 2
    l=mm-1;
else
    l=mm;
end
%%
for m = 2:l 
    for n=1:length(avgIQ_data)-2 %need to start at second pulse to average over pulse 1
        for k=1:length(a)
            phaseIQ_Matrix(n,k) = 0;
           for j=1:length(b)
               if n==1 && b(j) == 0
                   phaseIQ(j,k) = avgIQ_data(n-b(j)+2,...
                       m-a(k)).*conj(avgIQ_data(n-b(j)+1,m-a(k)));
                   phaseIQ_Matrix(n,k) = phaseIQ_Matrix(n,k)+ phaseIQ(j,k); 
               elseif n==1 && b(j) ==1
                   phaseIQ(j,k) = avgIQ_data(n-b(j)-1,...
                       m-a(k)).*conj(avgIQ_data(n+b(j),m-a(k)));
                   phaseIQ_Matrix(n,k) = phaseIQ_Matrix(n,k) + phaseIQ(j,k);
               else
                   phaseIQ(j,k) = avgIQ_data(n-b(j)+1,...
                       m-a(k)).*conj(avgIQ_data(n-b(j),m-a(k)));
                   phaseIQ_Matrix(n,k) = phaseIQ_Matrix(n,k)+ phaseIQ(j,k);
               end
           end
          phaseDifference_Matrix(n,m) = phaseDifference_Matrix(n,m) + phaseIQ_Matrix(n,k);     
            
        end

    end
end
%%

velocityEstimate = (nlambda/(pi*T)).*unwrap(angle(V_dop));  % argument 

%V_SDopp = (nlambda/(pi)).*unwrap(angle(V_dop));  %TEST



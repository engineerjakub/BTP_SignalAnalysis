function V_SDopp = VSDopp_v2_BrainTV(T,av_data,nlambda,mm)
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
%               av_data - IQ data averaged over a user defined number of
%               samples
%               nlambda - lambda/4 * 1000 (Nyquist m)
%               mm - no. of gates selected by the user to analyse
%
% *OUTPUT:*     D_SDopp- displacement in microns estimated from SDopp
%               autocorrelator velocity estimator
%                         
%
% *AUTHOR(S):*
%
%               Dr. C. Banahan updated 24/06/2016
%
%--------------------------------------------------------------------------
%%

%For the first time lag we average over 4 IQ sample volumes and over two
%pulse lengths. Subsequent times are averaged over 3 IQ sample volumes and
%over 1 pulse length.
%%
PL = 2;
ND = 4;
b = 2 - PL:1:0;
a = 2 - ND:1:0;

VD = zeros(length(av_data)-2,length(a));
v = zeros(length(av_data)-2,2);
V_dop = zeros(length(av_data)-2,mm);

for n = 1:length(av_data)-2
    m = 1;
    for k = 1:length(a)
        VD(n,k) = 0;
        for j=1:length(b)
               if n==1 && b(j) == 0
                   v(j,k) = av_data(n+2-b(j),...
                       m-a(k)).*conj(av_data(n+1-b(j),m-a(k)));
                   VD(n,k) = VD(n,k)+ v(j,k); 
               elseif n==1 && b(j) ==1
                   v(j,k) = av_data(n-1+b(j),...
                       m-a(k)).*conj(av_data(n-1 + b(j)+1,m-a(k)));
                   VD(n,k) = VD(n,k) + v(j,k);
               else
                   v(j,k) = av_data(n-b(j)+1,...
                       m-a(k)).*conj(av_data(n-b(j),m-a(k)));
                   VD(n,k) = VD(n,k)+ v(j,k);
               end
        end
        V_dop(n,m) = V_dop(n,m) + VD(n,k); 
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
    for n=1:length(av_data)-2 %need to start at second pulse to average over pulse 1
        for k=1:length(a)
            VD(n,k) = 0;
           for j=1:length(b)
               if n==1 && b(j) == 0
                   v(j,k) = av_data(n+2-b(j),...
                       m-a(k)).*conj(av_data(n+1-b(j),m-a(k)));
                   VD(n,k) = VD(n,k)+ v(j,k); 
               elseif n==1 && b(j) ==1
                   v(j,k) = av_data(n-1+b(j),...
                       m-a(k)).*conj(av_data(n-1 + b(j)+1,m-a(k)));
                   VD(n,k) = VD(n,k) + v(j,k);
               else
                   v(j,k) = av_data(n-b(j)+1,...
                       m-a(k)).*conj(av_data(n-b(j),m-a(k)));
                   VD(n,k) = VD(n,k)+ v(j,k);
               end
           end
          V_dop(n,m) = V_dop(n,m) + VD(n,k);     
            
        end

    end
end
%%

V_SDopp = (nlambda/(pi*T)).*unwrap(angle(V_dop));  % argument 

%V_SDopp = (nlambda/(pi)).*unwrap(angle(V_dop));  %TEST



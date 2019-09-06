function [signalbp] = iirhighpass2( fs, signal, highpassfreq, highpassstop )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%Set IIR Filter Parameters
highpass.HiStop                 = highpassstop          ;    % Hz 
highpass.HiPass                 = highpassfreq          ;    % Hz
highpass.Rp                     = 3                     ;    % dB passband ripple
highpass.Rs                     = 30                    ;    % dB attenuation in stopband
highpass.ID                     = 'highpass'            ;    % Identifier 
highpass.perform                = true                  ;    % toggle to apply or not apply 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HIGH PASS FILTERING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

% zero out the signal at its onset
fprintf(1,'Zeroing signal start-value to prevent filter instabilities at signal onset...\n');
signal = subtractfirstsamplevalue( signal );
    
% IIR bandpass definition using Butterworth filter
fprintf(1, 'Defining passband ripple and stopband attentuation...\n');
hp.Wp                          = [ highpass.HiPass ]/(fs/2); 
hp.Ws                          = [ highpass.HiStop ]/(fs/2);
hp.Rp                          = highpass.Rp; 
hp.Rs                          = highpass.Rs;

% calculate the minimum filter order
fprintf(1, 'Calculating minimum filter order...\n');
[hp.n,hp.Wn]                   = buttord( hp.Wp, hp.Ws, hp.Rp, hp.Rs);
 hp.n                          = hp.n + rem(hp.n,2);

% caclulate the filter coefficients in Zero-Pole-Gain design
%         [bp.z, bp.p, bp.k]             = butter( bp.n, bp.Wn, 'bandpass' );
fprintf(1, 'Calculating filter coefficients in Zero-Pole-Gain design...\n');
[hp.z, hp.p, hp.k]             = butter( hp.n, hp.Wn, 'high' );
[hp.sos, hp.g]                 = zp2sos( hp.z, hp.p, hp.k );
 hp.h                          = dfilt.df2sos( hp.sos, hp.g );

% do band-pass filtering...       
fprintf(1, 'Created Butterworth IIR highpass filter of order %d (n = %d filter coefficients)\n',hp.n,hp.n*2);
fprintf(1, 'High-pass filtering signal at %1.1f Hz\n',highpass.HiPass);
warning('off', 'signal:filtfilt:ParseSOS');

% check if a,b or zero-pole-gain filter design
if isfield(hp,'b'), 
    % filter signal using a,b coeffiienct 
    signalbp = single(filtfilt(hp.b,1,double(signal)));
else
    % filter signal using zero-pole-gain coefficients
    signalbp = single(filtfilt(hp.sos,hp.g,double(signal)));
end
    

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        


end

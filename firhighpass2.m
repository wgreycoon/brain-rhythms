function [signalbp] = firhighpass2(fs, signal, highpassfreq, highpassstop)

%Set FIR Filter Parameters
hipfiltobj.stopband             = highpassstop; %to nearest 1/10th Hz (target frequency for chosen level of attenuation)
hipfiltobj.passband             = highpassfreq;   %in Hz (target frequency above which no attenuation is desired)
hipfiltobj.passbandripple       = 0.1;            %normalized frequency, (*pi rad/sample)
hipfiltobj.stopbandattenuation  = 65;             %in percent, desired level of attenuation at stopband frequency
hipfiltobj.designmethod         = 'kaiserwin';    %for FIR hipass, could be this, or 'equiripple'


%% Build high pass and low pass filter objects to combine for band passing...

%design high pass filter
hpFilt = designfilt('highpassfir', ...
                    'SampleRate',          fs, ...
                    'StopbandFrequency',   hipfiltobj.stopband, ...
                    'PassbandFrequency',   hipfiltobj.passband, ...
                    'PassbandRipple',      hipfiltobj.passbandripple*fs/2, ...
                    'StopbandAttenuation', hipfiltobj.stopbandattenuation, ...
                    'DesignMethod',        hipfiltobj.designmethod);       
                
% hp = fvtool(hpFilt);                   
% addfilter(hp, lpFilt);


%% Zero data

    fprintf('Zeroing signal start-value to prevent filter instabilities at signal onset...')
[ signalzeroed ] = subtractfirstsamplevalue( signal );
    fprintf('done.\n')
    

%% High pass filter data with zero-phase-lag

% then high pass the low-passed data...
fprintf('High passing data with FIR filter at %2.1fHz...',highpassfreq)
signalbp = repmat(signalzeroed,1);
for ch_idx = 1:size(signal, 2)
    signalbp(:,ch_idx) = filter( hpFilt, signalzeroed(:,ch_idx) );
    signalbp(:,ch_idx) = filter( hpFilt, subtractfirstsamplevalue(flipud(signalbp(:,ch_idx))) );
    signalbp(:,ch_idx) = flipud( signalbp(:,ch_idx) );
end 
fprintf('done.\n')


    function [ dataout ] = subtractfirstsamplevalue( datain )

    %ensure expected dimensions
    if size(datain,2) > size(datain,1)
        datain = datain';
    end

    %subtract value of first sample from each time series (so all signals start at zero)
    d = repmat(datain,1);
    for idx = 1:size(datain,2)
        d(:,idx) = datain(:,idx) - datain(1,idx);
    end

    dataout = d; clear datain d

    end


end
































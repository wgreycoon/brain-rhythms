function [signalbp] = firbandpass2(fs, signal, highpassfreq, lowpassfreq, highstop, lowstop)
%% Coupling Palette 
% 
% -------------------------------------------------------------------------
% FIR band pass filtering
% ------------------------------------------------------------------------- 
%
% Dependencies
% ---------------------------
% functions:
%     subtractfirstsamplevalue.m
%
% toolbox:
%     MATLAB filter toolbox?
%
% ---------------------------
%
% Input 
% ---------------------------
%     fs           :    sampling frequency (in Hz)
%     signal       :    
%     highpassfreq :    
%     lowpassfreq  :    
%     highstop     :
%     lowstop      :    
%     
%
%     Example = couplingpalette( 100, ecog_data, 10, 15, 8, 17 );
%
% ---------------------------
%
% Output
% ---------------------------
%     signalbp     :    (time x channels) matrix 
%
% ---------------------------
%
% Author(s):   William Coon   wcoon@mgh.harvard.edu
%
% Last edited: March 17, 2019 1246hrs
%
% -------------------------------------------------------------------------
%% 

detrend.highpassfreq            = 0.1       ;    %in Hz, pass band frequency cutoff for IIR (butterworth) high pass filter 
detrend.filterorder             = 4          ;    %filter order

%Set FIR Filter Parameters
hipfiltobj.stopband             = highstop; %manual
hipfiltobj.passband             = highpassfreq;   %in Hz (target frequency above which no attenuation is desired)
hipfiltobj.passbandripple       = 0.1;            %normalized frequency, (*pi rad/sample)
hipfiltobj.stopbandattenuation  = 65;             %in percent, desired level of attenuation at stopband frequency
hipfiltobj.designmethod         = 'kaiserwin';    %for FIR hipass, could be this, or 'equiripple'

lowpfiltobj.stopband            = lowstop;  %manual
lowpfiltobj.passband            = lowpassfreq;    %in Hz (target frequency below which no attenuation is desired)
lowpfiltobj.passbandripple      = 0.1;            %normalized frequency, (*pi rad/sample), ripple in passband
lowpfiltobj.stopbandattenuation = 65;             %in percent, desired level of attenuation at stopband frequency
lowpfiltobj.designmethod        = 'kaiserwin';    %for FIR hipass, could be this, or 'equiripple'

if lowpfiltobj.stopband > fs/2
    lowpfiltobj.stopband = fs/2-1;
end


%% Build hipass and low pass filter objects to combine for band passing...

%design high pass filter
hpFilt = designfilt('highpassfir', ...
                    'SampleRate',          fs, ...
                    'StopbandFrequency',   hipfiltobj.stopband, ...
                    'PassbandFrequency',   hipfiltobj.passband, ...
                    'PassbandRipple',      hipfiltobj.passbandripple*fs/2, ...
                    'StopbandAttenuation', hipfiltobj.stopbandattenuation, ...
                    'DesignMethod',        hipfiltobj.designmethod);  

%design low pass filter
lpFilt = designfilt('lowpassfir', ...
                    'SampleRate',          fs, ...
                    'StopbandFrequency',   lowpfiltobj.stopband, ...
                    'PassbandFrequency',   lowpfiltobj.passband, ...
                    'PassbandRipple',      lowpfiltobj.passbandripple*fs/2, ...
                    'StopbandAttenuation', lowpfiltobj.stopbandattenuation, ...
                    'DesignMethod',        lowpfiltobj.designmethod);      
                
% hp = fvtool(hpFilt);                   
% addfilter(hp, lpFilt);


%% Zero data

    fprintf('Zeroing signal start-value to prevent filter instabilities at signal onset...')
[ signalzeroed ] = subtractfirstsamplevalue( signal );
    fprintf('done.\n')
    

%% De-trend data with IIR high pass 

%      fprintf('De-trending data with IIR high pass at %1.1fHz...',detrend.highpassfreq)
%  [ b, a ] = butter(detrend.filterorder, detrend.highpassfreq/(fs/2), 'high');
%  signaldt = filtfilt(b, a, double(signalzeroed)); signaldt = single(signaldt); 
%      fprintf('done.\n')
    signaldt = signalzeroed;


%% Filter data with zero-phase-lag

% first low pass the data...
fprintf('Low passing data with FIR filter at %dHz...',lowpassfreq)
signallp = repmat(signaldt,1);
for ch_idx = 1:size(signal, 2)
    signallp(:,ch_idx) = filter( lpFilt, signaldt(:,ch_idx) );
    signallp(:,ch_idx) = filter( lpFilt, subtractfirstsamplevalue(flipud(signallp(:,ch_idx))) );
    signallp(:,ch_idx) = flipud( signallp(:,ch_idx) );
end %clear signaldt
fprintf('done.\n')    

% then high pass the low-passed data...
fprintf('High passing data with FIR filter at %dHz...',highpassfreq)
signalbp = repmat(signallp,1);
for ch_idx = 1:size(signal, 2)
    signalbp(:,ch_idx) = filter( hpFilt, signallp(:,ch_idx) );
    signalbp(:,ch_idx) = filter( hpFilt, subtractfirstsamplevalue(flipud(signalbp(:,ch_idx))) );
    signalbp(:,ch_idx) = flipud( signalbp(:,ch_idx) );
end %clear signallp
fprintf('done.\n')

% fprintf('Completed alpha band extraction with FIR band pass procedure.\n')







% 
% function [hpFilt, lpFilt] = buildFIRfilters( fs, hipfiltobj, lowpfiltobj )
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %
% % This function is designed to build a low pass filter object and a high
% % pass filter object for "combining" to band pass a signak outside of this
% % function
% %
% %% Build hipass and low pass filters to combine for band passing...
% 
% %design high pass filter
% hpFilt = designfilt('highpassfir', ...
%                     'SampleRate',          fs, ...
%                     'StopbandFrequency',   hipfiltobj.stopband, ...
%                     'PassbandFrequency',   hipfiltobj.passband, ...
%                     'PassbandRipple',      hipfiltobj.passbandripple*fs/2, ...
%                     'StopbandAttenuation', hipfiltobj.stopbandattenuation, ...
%                     'DesignMethod',        hipfiltobj.designmethod);
% % hp = fvtool(hpFilt);     
% 
% %design low pass filter
% lpFilt = designfilt('lowpassfir', ...
%                     'SampleRate',          fs, ...
%                     'StopbandFrequency',   lowpfiltobj.stopband, ...
%                     'PassbandFrequency',   lowpfiltobj.passband, ...
%                     'PassbandRipple',      lowpfiltobj.passbandripple*fs/2, ...
%                     'StopbandAttenuation', lowpfiltobj.stopbandattenuation, ...
%                     'DesignMethod',        lowpfiltobj.designmethod);               
% lp = fvtool(lpFilt);  
% addfilter(hp, lpFilt);
% 
% end



%%

% % % ch = 21;
% % % close all; figure; plot(signalbp(12001:13200,ch))
% % % % [data, ~,~] = load_bcidat('/Users/wcoon/Dropbox/LabWork/Packages/FBO Analysis Package/data/core/AMC045/ECOGS001R03.dat');
% % % % hold on; plot(signallp(12001:14400,ch),'r')
% % % hold on; plot(signaldt((12001:13200),ch),'r')

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
































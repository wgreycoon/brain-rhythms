function [ convdata ] = wavbank3( sigdata, wavfreqs, fs, n, headerstr  )
%% Wavelet bank filtering 
% 
% -------------------------------------------------------------------------
% 'wavbank' and 'wavbank3' perform wavelet-based convolution (filtering) on
% an input time series.  Complex morlet wavelets are created, using the
% parameter 'n' (# of cycles) to control the time-frequency trade-off in
% specificity of time-frequency analysis (i.e., smaller 'n' gives higher
% temporal resolution at the expense of frequency resolution, while higher
% 'n' gives more stable frequency estimates at the expense of more
% smoothing over time (less temporal resolution). It returns complex-valued
% time series for each frequency requested, containing information about
% instantaneous PHASE and POWER.  You can use "abs(convdata)" to get the
% envelope of the signal, or "angle(convdata)" to get the instantaneous
% phase time series.
%
% these functions are based on code from Mike X. Cohen's delightful book on
% "Analyzing Neural Time Series Data".  I highly recommend it.
%
% 'wavbank3' performs the same signal processing as 'wavbank', but instead
% of reporting each convolution step, gives a progress bar inspired by code
% from:
% https://www.mathworks.com/matlabcentral/fileexchange/4863-ascii-progress-bar
% ------------------------------------------------------------------------- 
%
% Dependencies
% ---------------------------
% functions:
%     plotWavelet.m
%     constructWavelet.m
%
% ---------------------------
%
% Input 
% ---------------------------
%     fs           :    sampling frequency (in Hz)
%     sigdata      :    array of signal data (e.g., time x channels)
%     wavfreqs     :    vector of frequencies (in Hz) to return power for
%     n            :    scalar for wavelet construction (see above)
%     headerstr    :    string for progress bar feedback
%
%     Example = wavbank2( ecog_data, [1:30], 100, 7, 'Pz' );
%
% ---------------------------
%
% Output
% ---------------------------
%     convdata    :    (time x channels x frequency) matrix of complex
%                       values from wavelet convolution. Note that these
%                       are NOT NORMALIZED.
%
% ---------------------------
%
% Author(s):   William Coon   wcoon@mgh.harvard.edu
%
% Last edited: March 17, 2019 1231hrs
%
% -------------------------------------------------------------------------
%%

warning('No normalization, no 10*log10 conversion...');

% ensure data dimensions correct
dims = size(sigdata);
if dims(2) > dims(1)
    warning('Correcting dimensions...')
    sigdata = sigdata';
end

%pre-allocate
convdata = nan(size(sigdata,1),size(sigdata,2),length(wavfreqs));

%for status bar...
progress_bar_position = 0;
time_for_this_iteration = 0.01;
max = length(wavfreqs);    

    for wavfreq_idx = 1:length(wavfreqs)    
        
        thistic = tic;        
        
        %Construct wavelet
        wavcfg.srate                  = fs;                                         % get sampling rate
        wavcfg.halfwidth              = 1;                                          % in seconds, half-width of wavelet time series vector
        wavcfg.f                      = wavfreqs(wavfreq_idx);                      % frequency of wavelet in Hz
        wavcfg.n                      = n;%7;                                       % number of cycles in wavelet
        [o_wavelet, wavcfg]           = constructWavelet( wavcfg );

        %Plot Wavelet
        %plotWavelet( o_wavelet, wavcfg );

%         fprintf('Convolving signal data with %3.1f Hz complex wavelet...\n',wavcfg.f)

        for trial_idx = 1:size(sigdata,2)

            %Convolve wavelet with time series data
            % FFT parameters
            n_wavelet            = length(o_wavelet);
            n_data               = size(sigdata,1);
            n_convolution        = n_wavelet+n_data-1;
            half_of_wavelet_size = (length(o_wavelet)-1)/2;

            % FFT of wavelet and EEG data
            fft_wavelet = fft(o_wavelet,n_convolution);
            fft_data    = fft(sigdata(:,trial_idx),n_convolution,1); 

            % filtered signal from FFT (iFFT applied to product of multiplication in frequency domain)
            convolution_result_fft = ifft(repmat(fliplr(fft_wavelet)',1,size(sigdata(:,trial_idx),2)).*fft_data,n_convolution) * wavcfg.s;

            % cut off edges
            convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size,:);

            % get complex-valued time series containing phase and power information
            temppower = (reshape(convolution_result_fft,size(convolution_result_fft,1),size(convolution_result_fft,2))).^2;

            % store in struct
            convdata(:,trial_idx,wavfreq_idx)   = single(temppower);

        end
        
           progress_bar_position = progress_bar_position + 1 / max;
           clc;
           fprintf('|================ WAVELET BANK %s ====================|\n',headerstr);
           progress_string='|';       
           for counter = 1:floor(progress_bar_position * 100 / 2),
               progress_string = [progress_string, '#'];
           end
           disp(progress_string);
           disp(['|================== ',num2str(floor(progress_bar_position * 100)),'% completed ====================|']);
           % display progress per cent
           steps_remaining = max - wavfreq_idx;
           minutes = floor(time_for_this_iteration * steps_remaining / 60);
           seconds = rem(floor(time_for_this_iteration *  steps_remaining), 60);
           disp(' ');
           if (seconds > 9),
             disp(['            Estimated remaining time: ', num2str(minutes), ':', num2str(seconds)]);
           % show time indicators
           else
             disp(['            Estimated remaining time: ', num2str(minutes), ':0', num2str(seconds)]);
           end
           time_for_this_iteration = toc(thistic);

    end

end


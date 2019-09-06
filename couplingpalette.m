function [ couplingpalettevals ] = couplingpalette( phasesigdata, ampsigdata, phasewavfreqs, ampbpparams, fs, n )
%% Coupling Palette 
% 
% -------------------------------------------------------------------------
% coupling palette (see Kai Miller's work). this function takes input time
% series and analyzes Phase-Amplitude Coupling by binning power values by
% their associated corresponding phase values. 
%
% Complex morlet wavelets are created for PAC analysis, using the
% parameter 'n' (# of cycles) to control the time-frequency trade-off in
% specificity of time-frequency analysis (i.e., smaller 'n' gives higher
% temporal resolution at the expense of frequency resolution, while higher
% 'n' gives more stable frequency estimates at the expense of more
% smoothing over time (less temporal resolution). 
%
% amplitude time series are created with a bandpass (FIR) filter.
%
% ------------------------------------------------------------------------- 
%
% Dependencies
% ---------------------------
% functions:
%     constructWavelet.m
%     firbandpass2.m
%
% ---------------------------
%
% Input 
% ---------------------------
%     fs           :    sampling frequency (in Hz).
%     phasesigdata :    vector of phase-modulating channel time series.
%     ampsigdata   :    vector of power-modulated channel time series.
%     phasewavfreqs:    vector of frequencies of interest for phase.
%     ampbpparams  :    struct of parameters for band pass filtering. This
%                       will be the band of power modulated by phase.
%                             Example:
%                                 params.lowstop          = 8; 
%                                 params.lowpassfreq      = 10;     
%                                 params.highpassfreq     = 15;
%                                 params.highstop         = 17;
%     n            :    scalar for wavelet construction (see above)
%     
%
%     Example = couplingpalette( ecog_data, ecog_data, [1:30], [10 15 8 17], 7 );
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
% Last edited: September 6, 2019 1018hrs
%
% -------------------------------------------------------------------------
%% 
% warning('No normalization, no 10*log10 conversion...');

% ensure data dimensions correct
dims = size(phasesigdata);
if dims(2) > dims(1)
    warning('Correcting dimensions...')
    phasesigdata = phasesigdata';
end

% ensure data dimensions correct
dims = size(ampsigdata);
if dims(2) > dims(1)
    warning('Correcting dimensions...')
    ampsigdata = ampsigdata';
end

%% prepare signal to be analyzed for phase modulation

%get entire time series here...
bbandpwrall     = 10*log10(abs(hilbert(firbandpass2( fs, ampsigdata, ampbpparams(1), ampbpparams(2), ampbpparams(3), ampbpparams(4) ))).^2); %<--- SQUARED to get PWR

% %normalizing entire time series 
bbandpwrall     = (bbandpwrall - nanmean(bbandpwrall)) ./ nanstd(bbandpwrall);


%%

% set up for coupling pallettes...
num_bins             = 24;
bin_edges            = linspace(-pi, pi, num_bins + 1);
% bin_centers          = bin_edges(2:end)-(2*pi/num_bins/2);
couplingpalettevals  = nan(length(phasewavfreqs),num_bins);

fprintf('Convolving signal data from all channels with complex wavelets...\n')
for wavfreq_idx = 1:length(phasewavfreqs)            
    
    %NOW DO FOR ENTIRE TIME SERIES       
          
    %Construct wavelet
    wavcfg.srate                  = fs;                                         % get sampling rate
    wavcfg.halfwidth              = 1;                                          % in seconds, half-width of wavelet time series vector
    wavcfg.f                      = phasewavfreqs(wavfreq_idx);                      % frequency of wavelet in Hz
    wavcfg.n                      = n;                                          % number of cycles in wavelet
    [o_wavelet, wavcfg]           = constructWavelet( wavcfg );
        
    %Convolve wavelet with time series data
    % FFT parameters
    n_wavelet                     = length(o_wavelet);
    n_data                        = size(phasesigdata,1);
    n_convolution                 = n_wavelet+n_data-1;
    half_of_wavelet_size          = (length(o_wavelet)-1)/2;
    
    fprintf('%02.1f Hz wavelet...',wavcfg.f)  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ECoG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % FFT of wavelet and EEG data
        fft_wavelet               = fft(o_wavelet,n_convolution);
        fft_data                  = fft(phasesigdata,n_convolution,1);

        convolution_result_fft    = ifft(fliplr(fft_wavelet)'.*fft_data,n_convolution) * wavcfg.s;

        % cut off edges
        convolution_result_fft    = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size,:);

        % make LF phase & BB amplitude vectors of all trials in this channel
        phase_result_vector       = angle(convolution_result_fft);
        amp_result_vector         = bbandpwrall;

        % seperate the signal into evenly sized bins
        [~, which_bin]            = histc(phase_result_vector, bin_edges);

        % tally up bin values
        for idx = 1:num_bins

            flagBinMembers        = (which_bin == idx);
            bin_members           = amp_result_vector(flagBinMembers);
            couplingpalettevals(wavfreq_idx,idx)  = mean(bin_members);

        end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end ECoG %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('done.\n')
    
end
% %create vectors of 'Zmod' scores, complex numbers derived from the mean ripple band pwr associated with each phase bin
% zmodsbyphasebin  = couplingpalettevals .* exp(1i*repmat(bin_centers,size(couplingpalettevals,1),1));
% zmods2byphasebin = phaseavgdpwrvals2all .* exp(1i*repmat(bin_centers,size(phaseavgdpwrvals2all,1),1));
% 
% %calculate mean coupling vectors, 'Zmodmean')
% zmodmeans        = mean(zmodsbyphasebin,2);
% zmodmeans2       = mean(zmods2byphasebin,2);
% 
% %do statistics now
% for ix = 1:length(phasewavfreqs)
%     p_rippleband1(ix) = chi2cdf(sum( (couplingpalettevals(ix,:)-mean(couplingpalettevals(ix,:))).^2/mean(couplingpalettevals(ix,:)) ),num_bins-1,'upper'); %#ok
%     p_rippleband2(ix) = chi2cdf(sum((phaseavgdpwrvals2all(ix,:)-mean(phaseavgdpwrvals2all(ix,:))).^2/mean(phaseavgdpwrvals2all(ix,:))),num_bins-1,'upper'); %#ok
% end
% 
% %bonferonni correction
% p_rippleband1 = p_rippleband1*length(phasewavfreqs);
% p_rippleband2 = p_rippleband2*length(phasewavfreqs);
% p_rippleband1(p_rippleband1==0) = NaN;
% p_rippleband2(p_rippleband2==0) = NaN;

    function [wavelet_out, cfg] = constructWavelet( cfg )
        %% constructWavelet.m
        %
        % constructs a wavelet with parameters from struct 'cfg'
        % parameters...
        % srate = cfg.srate;        % sampling rate in Hz
        % f     = cfg.f;            % frequency of wavelet in Hz
        % n     = cfg.n;            % number of cycles in wavelet
        % time  = cfg.halfwidth;    % time, from -1 to 1 second in steps of 1/sampling-rate
        % s     = n/(2*pi*f);       % s = standard deviation of Gaussian, derive from 'n' number of cycles
        %
        % % % % %

            % parameters...
            srate = cfg.srate;                                  % sampling rate in Hz
            f     = cfg.f;                                      % frequency of wavelet in Hz
            N     = cfg.n;                                      % number of cycles in wavelet
            time  = -cfg.halfwidth:1/srate:cfg.halfwidth;       % time, from -t to t second in steps of 1/sampling-rate
            s     = N/(2*pi*f);                                 % s = standard deviation of Gaussian, derive from 'n' number of cycles

            % define complex Gaussian
            gausskernel = exp(-time.^2./(2*s^2));

            % define complex sine wave
            swave       = exp(2*pi*1i*f.*time);

            % define frequency-specific scaling factor
            A = sqrt(1/(s*sqrt(pi)));

            % and together they make a wavelet
            wavelet_out = gausskernel .* swave; 

        %     % scale wavelet 
        %     wavelet_out = A * wavelet_out;

            % add to struct
            cfg.time = time;
            cfg.s    = s;
            cfg.A    = A;
    end


end


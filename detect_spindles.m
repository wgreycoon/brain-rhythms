function [spindle_det, filtTraces] = detect_spindles( Data, ch_names, fs, bpparms, doWavelet, ampl_factor, minmaxdurs, sigproc, manualthreshs)
%%
% Well, guess there's some work to be done on the documentation side of
% things here....

% As of now, EVERYTHING works EXCEPT the wavelet-derived frequency measures
% (center frequency and frequency drift).  DO NOT TRUST THOSE NUMBERS I
% HAVEN'T DEBUGGED IT YET!!!!
%
% Everything else purrs like a kitten
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Will Coon, PhD, wcoon@mgh.harvard.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
warning('off','all'); warning
fprintf('============================================================\n')
fprintf('============================================================\n')
fprintf('============================================================\n')

fprintf('This version was last edited August 9th, 2019\n')    

fprintf('============================================================\n')
fprintf('============================================================\n')
fprintf('============================================================\n')

    %% Setup
    
    if ~exist('manualthreshs','var')
        manualthreshs = [];
    end
    
    % decide if using CWT (wavelet xform) or filter-Hilbert (band pass w/ Hilbert xform)
    if strcmp(sigproc,'filter-hilbert')%~doWavelet    %exist('bpparams','var') && isnumeric(bpparams) && length(bpparams) == 4
        disp('Using filter-Hilbert method for signal processing...')
        fprintf('Pass band will be %1.1f to %1.1f Hz,... \nStop bands below %1.1fHz and above %1.1fHz...\n', bpparms(1), bpparms(2), bpparms(3), bpparms(4))
            fc = (bpparms(1)+bpparms(2)) / 2;
            fb = nan;
    elseif strcmp(sigproc,'wavelet')
        disp('Using continuous wavelet xform method for signal processing...')                  
            fc = 13.5; % Mother wavelet center frequency
            fb = .5; % Mother wavelet bandwidth            
        fprintf('Center frequency of wavelet will be %1.1fHz, \nband width parameter set at fb=%1.1f...\n', fc, fb )
    else
        error('Must specify ''filter-hilbert'' or ''wavelet'' for signal processing.')
    end
    
    % catch wavelet settings, since this input field used to be "tthresh" 
    % in deprecated copies of this function
    if ~exist('doWavelet','var')
        doWavelet = false;
        fprintf('No wavelet parameter specified, defaulting to NO wavelet analyses...\n');
    elseif exist('doWavelet','var') && ~(doWavelet==1 || doWavelet==0)
        doWavelet = false;
        fprintf('''doWavelet parameter'' found but is not boolean. Wavelet analysis set to default:  %s\n', doWavelet );
    else
        fprintf('Wavelet (characterization) analysis set to:  %s\n', num2str(doWavelet) );
    end
    
    % ensure data dimensions correct
    dims = size(Data);
    if dims(2) > dims(1)
        Data = Data';
    end

    % pre-allocation ---
    % for complex-valued hilbert transform of data time series...
    Data_cplx                = double(nan(size(Data))); 
    
    % for square of the amplitude envelope from hilbert transform, from
    % which the amplitude threshold will be derived and on which 
    % detection will operate...
    Data_sqenvff             = double(nan(size(Data)));

    % smoothing window for later...
    window_size_in_millis    = 100; %create 100ms window to convolve with
    window_size_in_samples   = round((window_size_in_millis/1000)*fs);
    mov_avg_filter_kernel    = (1/window_size_in_samples)*ones( window_size_in_samples, 1 ); %the scaling coefficients make the kernel have unit area

    
    %% Compute Separately for Each Channel...
   
    for ch=1:size(Data,2) % Loop for each signal
        
        fprintf('\n====> Working on Channel %s.\n',ch_names{ch});
        
        
        %============================================================   
        % ------- filter-Hilbert signal processing
        if strcmp(sigproc,'filter-hilbert')
            
            % Get complex-valued time series from hilbert-transformed, sigma band-passed signal 
            fprintf('Applying bandpass and hilbert transform transform...\n  [')
            Data_cplx(:,ch)      = hilbert(firbandpass2( fs, double(Data(:,ch)), bpparms(1), bpparms(2), bpparms(3), bpparms(4) )); fprintf('.') % instead of wavelet tranformation

            % Apply 100ms moving-average filter kernel with zero-phase lag two-way filtfilt
            fprintf('Applying %dms moving average filter to \n squared envelope of wavelet-filtered data...\n  [',window_size_in_millis)
            Data_sqenvff(:,ch)   = filtfilt( mov_avg_filter_kernel, 1, double(abs( Data_cplx(:,ch) ).^2) ); fprintf('.')
        
        % --------------------------------------------   
        % ------- wavelet-based signal processing    
        elseif strcmp(sigproc,'wavelet')

            % The appropriate scale depends on sampling frequency -- Use function 'SCAL2FRQ' to figure out the relationship between frequency and scale at your sampling rate
            % scale = 7.4; % Specifies scale to examine at 200Hz -- 

            % Feeding in the frequency gives you the scale, feeding in the scale gives you the frequency
            scale = scal2frq(fc,['cmor' num2str(fc) '-' num2str(fb)],1/fs);
            %   freq = scal2frq(scale,['cmor' num2str(fc) '-' num2str(fb)],1/Fs);
            
            current_data       = Data(:,ch);
            EEGWave            = cwt(current_data, scale, ['cmor' num2str(fc) '-' num2str(fb)] ); % conducts the wavelet tranformation on workspace variable "EEG"
%             EEGChannel         = real(EEGWave.^2); % Takes only the real component of the coefficients
            EEGChannel         = real(EEGWave); % Takes only the real component of the coefficients            
            Data_sqenvff(:,ch) = EEGChannel;

            % Set the Sampling Frequency and take Moving Average
            Data_sqenvff(:,ch) = Data_sqenvff(:,ch).^2;
%             window             = ones((fs/10),1)/(fs/10); % create 100ms window to convolve with
%             Data2              = filter( window, 1, Data_sqenvff(:,ch) ); % take the moving average using the above window
%             Data3              = filtfilt( window, 1, Data_sqenvff(:,ch) ); % take the moving average using the above window
            Data4              = filtfilt( mov_avg_filter_kernel, 1, Data_sqenvff(:,ch) ); % take the moving average using the above window
            
            Data_sqenvff(:,ch) = Data4;
            Data_cplx(:,ch)    = EEGWave;
            
            %cleanup
            clear EEGChannel Data1 Data2 Data3 Data4
            
        end
        % ------------------------------
        %============================================================   
        
        % Calculate some relevant numbers...
        signalmedian(ch)     = median(Data_sqenvff(:,ch));    %#ok compute median amplitude 
        signalmean(ch)       = mean(Data_sqenvff(:,ch));      %#ok compute mean amplitude 
        signalstd(ch)        = std(Data_sqenvff(:,ch));       %#ok compute std
        
        % Here's the most important one...
        if isempty(manualthreshs)
            threshold(ch)    = signalmedian(ch).*ampl_factor; %#ok % defines the threshold
            fprintf('Using auto thresholds... %2.3fuV\n',threshold(ch))
        else
            threshold(ch)    = manualthreshs(ch); %#ok % defines the threshold
            fprintf('Using manually set thresholds... %2.3fuV\n',threshold(ch))
        end
        
        
        %% Detect spindles
        
        % grab data from one channel
        current_data = Data_sqenvff(:,ch);

        % mark channel as bad if NaN's found in signal
        if any(isnan(current_data)), spindle_det(ch).bads = 1; else spindle_det(ch).bads = 0; end %#ok
        
        %===============================    
        % ------- detection prep
        over = current_data>threshold(ch);        %#ok Mark all points over threshold as '1'
        locs = (zeros(1,length(current_data)))';  %#ok Create a vector of zeros the length of the MS signal

        verify = false; %(for debugging and visualizing each step...)
        
                if verify                       % if we want to visualize this example
                    plotwin = 1;                % milliseconds
                    plotixs(1) = 1;             % 3451; 
%                     plotixs(2) = length(current_data);
                    plotixs(2) = plotixs(1) + fs*plotwin;
                    t=0:1/fs:(plotixs(2)-plotixs(1))/fs; whos t
                    close all; figure; set(gcf,'position',[49         479        1826         567]); setPaperSize(gcf); %fs = Fs;
                    clf
                    plot(t,real(Data_cplx(plotixs(1):plotixs(2),ch))); axis tight
                    hold on; plot(t,Data_sqenvff(plotixs(1):plotixs(2),ch),'r')
                    set(gca,'fontsize',16); xlabel('Time (ms)')
%                     plot(abs(Data_bp(plotixs(1):plotixs(2),ch)),'g')
                end
                
                
        %===============================    
        % ------- Mark all points over threshold (with a numeral '1')
        over = current_data>threshold(1,ch);                
                if verify, overtmp = over(plotixs(1):plotixs(2)); %#ok (only do this if verify flag is on)
                hold on; plot(t,overtmp*280,'k'), end

                    
        %===============================    
        % ------- find signal peaks within each identified window > threshold   
        [~,locs] = findpeaks(current_data);                           
                if verify, [~,locstmp] = findpeaks(current_data(plotixs(1):plotixs(2)));                
                hold on; stem(1000*locstmp/(fs*plotwin),270*ones(length(locstmp),1),'o'), end          

                    
        %===============================    
        % ------- correct in case (first window > threshold) starts at index (1) or (end of last) is at index (end)    
        over_idxs  = find(over>0);                          %#ok ==> sample indexes where envelope exceeds threshold
        start_idxs = find(diff(over)>0)+1;                  % beginning indexes of each segment of over_idxs
        if over(1)>0, start_idxs = [1; start_idxs]; end     %#ok ==> in case first sample is over threshold
        end_idxs   = find(diff(over)<0)+1;                  %  end indexes of each segment of over_idxs
        if length(end_idxs)<length(start_idxs), end_idxs = [end_idxs; length(current_data)]; end %#ok  if end of last cand is clipped, make end of signal the last sample


        %===============================    
        % -------ignore spindles that start in first 2 secs of data
        end_idxs(start_idxs<fs*2) = [];
        start_idxs(start_idxs<fs*2) = [];  


        %===============================    
        % ------- in case at end of record, "over" starts flipping to 1,
          % but reaches end of record before flipping back to 0
        locs = [locs; length(current_data)]; %#ok:: so this is a fix for that
        %-------------------------------


        %===============================    
        % ------- find index of max peak amplitude (in case there's more than one)
        peak_idxs = nan(size(start_idxs));
        for ix = 1:length(start_idxs)            
%                     num_peaks_in_this_candidate         = length( intersect(start_idxs(ix):end_idxs(ix),locs) );
            highest_peak_in_this_candidate_ampl = max( current_data( intersect(start_idxs(ix):end_idxs(ix),locs) ) );
            highest_peak_in_this_candidate_idx  = start_idxs(ix) + find( current_data( start_idxs(ix):end_idxs(ix) ) == highest_peak_in_this_candidate_ampl, 1, 'first' );
            peak_idxs(ix)                       = highest_peak_in_this_candidate_idx;  
        end       
                
                
            %===============================    
            % ------- accounts for spindle that starts right near the end
            % of the record and may not "end" before the data file
            % ends
            if peak_idxs(end) > length(current_data), peak_idxs(end) = length(current_data); end
            %-------------------------------                
                
                
        %===============================    
        % ------- Remove spindles whose first 400ms includes any
        % dips below threshold (to mimic how the original wavelet
        % detector marked candidate events by starting a count
        % from the first point it sees above threshold and only 
        % branding it as a spindle candidate if that counter hits
        % the min duration (ex. 400ms) without ever dipping below
        % threshold).
            fprintf('Ensuring that spindle candidates start with a minimum CONSECUTIVE amount of time above threshold...\n')                   
            i = length(start_idxs);
            origlen = i;
            while i ~= 1 
                if start_idxs(i)+minmaxdurs(1)*(fs/1000) < length(current_data) %only if won't hit end of record 
                    if any(current_data(start_idxs(i):start_idxs(i)+(floor(minmaxdurs(1)*(fs/1000))-1)) < threshold(ch)) %if any of the samples in the first x-ms are below threshold, toss this one
                        start_idxs(i)     = []; 
                        end_idxs(i)       = []; 
                        peak_idxs(i)      = [];
                    end
                else
                    if any(current_data(start_idxs(i):length(current_data)) < threshold(ch)) %if any of the samples in the first x-ms are below threshold, toss this one
                        start_idxs(i)     = []; 
                        end_idxs(i)       = []; 
                        peak_idxs(i)      = [];
                    end
                end
                i = i - 1;
            end
            fprintf('Number of detections reduced from %d to %d.\n',origlen,length(start_idxs))
        %===============================
               
                
                
        %===============================    
        % -------- DURATIONS (operates on 100ms averaged signal envelope, NOT the squared envelope but the envelope from the hilbert directly)
        ix = 1; buffersize = 2; 
        for ix = 1:length(start_idxs)    

            if peak_idxs(ix)-buffersize*fs < 1, ix1 = 1; else ix1 = peak_idxs(ix)-buffersize*fs; end %#ok
            if peak_idxs(ix)+buffersize*fs > length(current_data), ix2 = length(current_data); else ix2 = peak_idxs(ix)+buffersize*fs; end

            centerindex = peak_idxs(ix); % begin search for edges at predetermined peak
            dropmax     = sqrt(current_data(centerindex))/2; %see, it takes the square-root of the squared signal...

            i=centerindex;  % find ending edge
            while sqrt(current_data(i)) > dropmax
                if i==ix2
                    break
                else i=i+1;
                end
            end
            adj_end_idxs(ix)=i; %#ok

            i=centerindex;  % find beginning edge
            while sqrt(current_data(i)) > dropmax
                if i==1
                    break
                else i=i-1;
                end
            end
            adj_start_idxs(ix)=i; %#ok

            %store DURATION of spindle defined by Full-Width Half-Max (FWHM thresholding)
            spindle_det(ch).duration(ix) = (1000/fs)*(adj_end_idxs(ix) - adj_start_idxs(ix));

        end                
        %------------------------------
        
                
        %===============================    
        % ------- Merge any overlapping spindles...
            fprintf('Merging overlapping spindle events...\n')                   
            i = length(adj_start_idxs);
            origlen = i;
            while i ~= 1 
                if ~isempty(intersect( (adj_start_idxs(i):adj_end_idxs(i)), (adj_start_idxs(i-1):adj_end_idxs(i-1)) ))
                    adj_end_idxs(i-1)   = adj_end_idxs(i);
                    start_idxs(i)     = []; 
                    end_idxs(i)       = []; 
                    peak_idxs(i)      = [];
                    adj_start_idxs(i) = [];
                    adj_end_idxs(i)   = [];
                    spindle_det(ch).duration(i) = [];
                end
                i = i - 1;
            end
            fprintf('Number of detections reduced from %d to %d.\n',origlen,length(start_idxs))
        %===============================
                
                
% -------------------                 
        for ix = 1:length(start_idxs) 
% ------------------- 

        %-----------
        %===============================

            % Peaks and Troughs    
            spindle_det(ch).filtPeakIdxAmp(ix)   = max( real(Data_cplx(adj_start_idxs(ix):adj_end_idxs(ix),ch)) );
            spindle_det(ch).filtPeakIdx(ix)      = adj_start_idxs(ix)-1 + find( real(Data_cplx(adj_start_idxs(ix):adj_end_idxs(ix),ch))==max( real(Data_cplx(adj_start_idxs(ix):adj_end_idxs(ix),ch)) ), 1, 'first'  );    
            spindle_det(ch).filtTroughIdxAmp(ix) = min( real(Data_cplx(adj_start_idxs(ix):adj_end_idxs(ix),ch)) );
            spindle_det(ch).filtTroughIdx(ix)    = adj_start_idxs(ix)-1 + find( real(Data_cplx(adj_start_idxs(ix):adj_end_idxs(ix),ch))==min( real(Data_cplx(adj_start_idxs(ix):adj_end_idxs(ix),ch)) ), 1, 'first'  );    

            spindle_det(ch).rawPeakIdxAmp(ix)    = max( real(Data(adj_start_idxs(ix):adj_end_idxs(ix),ch)) );
            spindle_det(ch).rawPeakIdx(ix)       = adj_start_idxs(ix)-1 + find( real(Data(adj_start_idxs(ix):adj_end_idxs(ix),ch))==max( real(Data(adj_start_idxs(ix):adj_end_idxs(ix),ch)) ), 1, 'first' );    
            spindle_det(ch).rawTroughIdxAmp(ix)  = min( real(Data(adj_start_idxs(ix):adj_end_idxs(ix),ch)) );
            spindle_det(ch).rawTroughIdx(ix)     = adj_start_idxs(ix)-1 + find( real(Data(adj_start_idxs(ix):adj_end_idxs(ix),ch))==min( real(Data(adj_start_idxs(ix):adj_end_idxs(ix),ch)) ), 1, 'first' );    

            spindle_det(ch).envPeakIdx(ix)       = adj_start_idxs(ix)-1 + find( real(Data_sqenvff(adj_start_idxs(ix):adj_end_idxs(ix),ch))==max( real(Data_sqenvff(adj_start_idxs(ix):adj_end_idxs(ix),ch)) ), 1, 'first' );    
            spindle_det(ch).envPeakIdxAmp(ix)    = max( real(Data_sqenvff(adj_start_idxs(ix):adj_end_idxs(ix),ch)) );    
            spindle_det(ch).envPeakIdxTraceVal(ix) = real(Data(spindle_det(ch).envPeakIdx(ix),ch));

        %-----------    
        %===============================

% -------------------                 
        end %end peaks and troughs
% ------------------- 



% ===== FFT for artifact rejection ==========
% -------------------                 
        for ix = 1:length(start_idxs) 
% ------------------- 

            % --------spectrum of RAW data (for getting non-wavelet freq and artefact rejection)

            spiRaw              = real(Data(adj_start_idxs(ix):adj_end_idxs(ix),ch)); % data(spiStartIndex:spiEndIndex);
            ns                  = length(spiRaw);                                     % number of samples analyzed

            %length for FFT
            fftLength           = 2^nextpow2((fs*minmaxdurs(2)/1000));

            %FFT with zero-padding::
            %fft generates as many (sets of real and imaginary) coefficients as
            %there are data points (including zero-padding)
            %normalize by actual spindle length (otherwise longer spindle=more
            %power)
            spiFFT              = fft(spiRaw,fftLength)/ns;          %
            spiFFTAmp           = abs(spiFFT(1:floor(fftLength/2))); % the magnitude of HALF the window length constitute the amplitude coefficients

            spiFFTPow           = spiFFTAmp.^2;

            %frequencies vary from 0 to half the sampling freq (nyquist) with steps half the number of samples in window.
            spindleFreqs        = linspace(0,1,fftLength/2)*(fs/2);
            lowFrInd            = find(spindleFreqs>=bpparms(1),1,'first');
            highFrInd           = find(spindleFreqs<=bpparms(2),1,'last');

            if highFrInd < lowFrInd
                highFrInd       = lowFrInd;
            end

            %artifact frequencies
            %broadband frequency range that should not show power increase during
            %spindles
            lofreq_noise        = 20;
            hifreq_noise        = 80;
            low_art             = find(spindleFreqs>(lofreq_noise),1,'first');
            high_art            = find(spindleFreqs<(hifreq_noise),1,'last');

            %power at all frequencies in broadband range
            artpow              = (spiFFTPow(low_art:high_art-1));                    

            [max_spipow, spiMaxPowInd]  = max(spiFFTPow(lowFrInd:highFrInd));
            spindle_det(ch).freqFFT(ix) = spindleFreqs(lowFrInd+spiMaxPowInd-1);                    

            too_much_artifact(ix) = any(artpow>max_spipow); %#ok
            %too_much_artifact = too_much_artifact';
            
        % -----------
        % ------------------------------
        %===============================
        
        %Inter-event intervals (added July 12, 2019)
            if length(spindle_det(ch).filtPeakIdxAmp) ~= 1   
                if ix ~= 1 && ix ~= length(start_idxs)
                    spindle_det(ch).isi_previous_event_in_millis(ix) = (start_idxs(ix) - start_idxs(ix-1)) / (fs/1000);
                    spindle_det(ch).isi_next_event_in_millis(ix)     = (start_idxs(ix+1) - start_idxs(ix)) / (fs/1000);
                    spindle_det(ch).isi_min_in_millis(ix)            = min(  [spindle_det(ch).isi_previous_event_in_millis(ix) spindle_det(ch).isi_next_event_in_millis(ix)] );
                    spindle_det(ch).isi_max_in_millis(ix)            = max(  [spindle_det(ch).isi_previous_event_in_millis(ix) spindle_det(ch).isi_next_event_in_millis(ix)] );
                    spindle_det(ch).isi_mean_in_millis(ix)           = mean( [spindle_det(ch).isi_previous_event_in_millis(ix) spindle_det(ch).isi_next_event_in_millis(ix)] );
                elseif ix == 1
                    spindle_det(ch).isi_previous_event_in_millis(ix) = NaN;
                    spindle_det(ch).isi_next_event_in_millis(ix)     = (start_idxs(ix+1) - start_idxs(ix)) / (fs/1000);
                    spindle_det(ch).isi_min_in_millis(ix)            = spindle_det(ch).isi_next_event_in_millis(ix);
                    spindle_det(ch).isi_max_in_millis(ix)            = spindle_det(ch).isi_next_event_in_millis(ix);
                    spindle_det(ch).isi_mean_in_millis(ix)           = spindle_det(ch).isi_next_event_in_millis(ix);
                elseif ix == length(start_idxs)
                    spindle_det(ch).isi_previous_event_in_millis(ix) = (start_idxs(ix) - start_idxs(ix-1)) / (fs/1000);
                    spindle_det(ch).isi_next_event_in_millis(ix)     = NaN;
                    spindle_det(ch).isi_min_in_millis(ix)            = spindle_det(ch).isi_previous_event_in_millis(ix);
                    spindle_det(ch).isi_max_in_millis(ix)            = spindle_det(ch).isi_previous_event_in_millis(ix);
                    spindle_det(ch).isi_mean_in_millis(ix)           = spindle_det(ch).isi_previous_event_in_millis(ix);
                end
            else
                spindle_det(ch).isi_previous_event_in_millis(ix) = nan;
                spindle_det(ch).isi_next_event_in_millis(ix)     = nan;
                spindle_det(ch).isi_min_in_millis(ix)            = nan;
                spindle_det(ch).isi_max_in_millis(ix)            = nan;
                spindle_det(ch).isi_mean_in_millis(ix)           = nan;
            end
        
        
% -------------------        
        end %end FFT & ISI (for each spindle event loop (for ix = 1:length(start_idxs))
% -------------------       
        


        %===============================
        % ---------  cull out spindles short than / longer than min/max duration settings ('minmaxdurs')
        % ---> also added removing spindles with broadband power increase greater than spindle power increase
        spins_good_lengths = (spindle_det(ch).duration > minmaxdurs(1) & spindle_det(ch).duration < minmaxdurs(2) & too_much_artifact==0);                                              
        clear too_much_artifact
        %===============================
        %is this used anymore?
        

        %=============================================================
        %=============================================================
        %--------------------wavelets-------------------------
            
            %===============================    
            % -------- SETUP: wavelet convolution for high-resolution frequency estimates needed for frequency shift
            %(this is a bottleneck in terms of time)
            if doWavelet==1

                fprintf('creating wavelets...\n')
                %---------wavelet parameters-----------------

                % frequency resolution (note: limited frequency resolution)
                wave_freq_res       = 0.25;
                wave_freq           = 5:wave_freq_res:20;
                num_frex            = length(wave_freq);

                % number of cycles (tradeoff between time and frequency precision)
                wavelet_cycle_range = [7 7] ;

                % gaussian width and time
                s                   = logspace(log10(wavelet_cycle_range(1)), log10(wavelet_cycle_range(2)), num_frex)./(2*pi.*wave_freq);

                % 4-s wavelet
                t                   = -2:1/fs:2;

                % lengths of kernel
                Ltapr               = length(t);

                % make wavelets
                wavelets            = zeros(num_frex,length(t));

                for fi=1:num_frex
                    wavelets(fi,:)  = exp(2*1i*pi*wave_freq(fi).*t).*exp(-t.^2./(2*s(fi)^2));
                end

                % data padding on each side of spindle to avoid edge artifacts
                waveletPadTime      = 0.5;
                waveletPadSample    = waveletPadTime*fs; 

            else
                waveletPadSample    = 0;
                
            end
            %------- END SETUP ------
            %===============================    
   
            
            %===============================
            %------- APPLY TO SPINDLES
            if doWavelet==1
                
                for spidx = 1:length(adj_start_idxs)

%                     expandedSpindle=data(spiStartIndex-waveletPadSample:spiEndIndex+waveletPadSample);
                    if adj_end_idxs(spidx)+waveletPadSample <= size(Data_cplx,1) && adj_start_idxs(spidx)-waveletPadSample > 0
                        expandedSpindle       = real(Data_cplx( adj_start_idxs(spidx)-waveletPadSample:adj_end_idxs(spidx)+waveletPadSample, ch ));
                    elseif adj_end_idxs(spidx)+waveletPadSample <= size(Data_cplx,1) && adj_start_idxs(spidx)-waveletPadSample <= 0
                        expandedSpindle       = real(Data_cplx( 1:adj_end_idxs(spidx)+waveletPadSample, ch ));                        
                        expandedSpindle       = [zeros( 1 - (adj_start_idxs(spidx)-waveletPadSample), 1 ); expandedSpindle ];  
                    elseif adj_end_idxs(spidx)+waveletPadSample >= size(Data_cplx,1) && adj_start_idxs(spidx)-waveletPadSample > 0
                        expandedSpindle       = real(Data_cplx( adj_start_idxs(spidx)-waveletPadSample:size(Data_cplx,1), ch ));
                        expandedSpindle       = [expandedSpindle; zeros( (adj_end_idxs(spidx)+waveletPadSample) - size(Data_cplx,1), 1 ) ];  %#ok
                    end
                    
                    %lengths of signal and kernel
                    Ldata                 =  length(expandedSpindle);
                    Lconv1                =  Ldata+Ltapr-1;
                    Lconv                 =  pow2(nextpow2(Lconv1));

                    %normalize in freq domain so results scale with original data
                    fftwavelets           = fft(wavelets,Lconv,2);
                    fftwavelets           = fftwavelets./repmat(max(fftwavelets,[],2),[1 size(fftwavelets,2)]);

                    % wavelet convolution is faster in freq domain
                    expandedSpindle_fft   = fft(expandedSpindle,Lconv);

                    filtered              = zeros(num_frex,Ldata);
                    amp                   = zeros(num_frex,Ldata);

                    %for each frequency bin
                    for fi= 1:num_frex

                        %multiply spectrum of data with wavelet kernel, inverse FFT
                        convolutionResult = ifft(expandedSpindle_fft'.*fftwavelets(fi,:),Lconv);
                        convolutionResult = convolutionResult(1:Lconv1);
                        convolutionResult = convolutionResult(floor((Ltapr-1)/2):end-1-ceil((Ltapr-1)/2));

                        %filtered signal is the real part
                        filtered(fi,:)    = real(convolutionResult);
                        
                        %amplitude is magnitude of the complex signal
                        amp(fi,:)         = abs(convolutionResult);
                        
                    end

                    timeFreqExpandedSpindle = amp;               

                    %---------wavelet-based frequency detection
                    %indices of low and hi freqs in wavelet freq bins
                    lofreq                = bpparms(1);                    % select desired spindle range based on stopband cutoffs
                    hifreq                = bpparms(2);                    % select desired spindle range based on stopband cutoffs
                    loWave                = dsearchn(wave_freq',lofreq);
                    hiWave                = dsearchn(wave_freq',hifreq);
                    newFrRange            = wave_freq(loWave):wave_freq_res:wave_freq(hiWave);

                    %remove padding to get freq with max power
                    spiTF                 = timeFreqExpandedSpindle(loWave:hiWave,waveletPadSample+1:end-waveletPadSample);

                    [~,maxTFInd]          = max(spiTF,[],1);
                    maxTFfr               = newFrRange(maxTFInd);

                    if verify
                        figure;
                        plot(real(Data_cplx( adj_start_idxs(spidx):adj_end_idxs(spidx), ch )))
                        hold on; plot(abs(Data_cplx( adj_start_idxs(spidx):adj_end_idxs(spidx), ch )),'r')
                        hold on; plot(peak_idxs(spidx)-adj_start_idxs(spidx),1e-5,'*')
                    end
                    
                  %======================================
                  %------- ASSIGN SPINDLE FREQUENCY -----
%                   spindle_det(ch).spindlesWaveletFreq(spidx) = maxTFfr(spiMaxEnvInd);
                  spiMaxEnvInd = peak_idxs(spidx)-adj_start_idxs(spidx);
                  spindle_det(ch).spindlesWaveletFreq(spidx) = maxTFfr( spiMaxEnvInd );
                  %======================================
            
                    
                    %----------wavelet-based frequency drift

                    %isolate half of spindle, centered around maximum
                    duration              = adj_end_idxs(spidx)-adj_start_idxs(spidx);
                    spiLength             = duration;%/1000*fs;
                    shortStart            = floor(spiMaxEnvInd-spiLength/3);
                    shortEnd              = ceil(spiMaxEnvInd+spiLength/3);

                    if shortStart<1
                        shortStart        = 1;
                    end
                    if shortEnd>spiLength
                        shortEnd          = spiLength;
                    end

                    shortSegInds          = round(shortStart:shortEnd);

                    %robust fit through maximum frequencies

                    YfitValues            = maxTFfr(shortSegInds);
                    XfitValues            = (1:length(YfitValues))/fs;

                    %coefficients for change in frequency (Hz) per s
                    coeffs                = robustfit(XfitValues,YfitValues);                    
                   
                  %=======================================
                  %------- ASSIGN SPINDLE FREQ SLOPE -----
                  spindle_det(ch).spindlesFrSlope(spidx) = coeffs(2);
                  %=======================================
                    
                    
                end %for spidx=1:length(adj_start_idxs)   
                %------ END SPINDLE LOOP -------
                %===============================
                
                
            %------------------------------    
            else %if ~doWavelet
                  
                for spidx = 1:length(adj_start_idxs)
                    spindle_det(ch).spindlesWaveletFreq(spidx)  = nan;
                    spindle_det(ch).spindlesFrSlope(spidx)      = nan;
                end               
                
                
            end %if doWavelet
            %------ END IF-DO COND'L -------
            %===============================
            
%             spindlesFFTfreq(spiCounter)=spindleFreqs(spiMaxPowInd+lowFrInd-1);%freq with highest amplitude
%             
%             
%             %---------------differentiate between slow and fast spindles if using broadband filter------------
%             
%             midfreq=events.midfreq;
%             
%             %use FFT-based frequency to determine slow/fast (higher
%             %granularity)
%             currSpiFreq=spindlesFFTfreq(spiCounter);          
%             
%             %slow spindles: 1, fast spindles: 2
%             if currSpiFreq<=midfreq
%                 spindlesType(spiCounter)=1;
%             elseif currSpiFreq>midfreq
%                 spindlesType(spiCounter)=2;
%             end
            
           
                

        %===============================
        % --------- transposing AND culling out by durations
        %===============================

        spindle_det(ch).filtPeakIdxAmp     = spindle_det(ch).filtPeakIdxAmp( spins_good_lengths )';     %#ok
        spindle_det(ch).filtPeakIdx        = spindle_det(ch).filtPeakIdx( spins_good_lengths )';        %#ok
        spindle_det(ch).filtTroughIdxAmp   = spindle_det(ch).filtTroughIdxAmp( spins_good_lengths )';   %#ok
        spindle_det(ch).filtTroughIdx      = spindle_det(ch).filtTroughIdx( spins_good_lengths )';      %#ok 

        spindle_det(ch).rawPeakIdxAmp      = spindle_det(ch).rawPeakIdxAmp( spins_good_lengths )';      %#ok
        spindle_det(ch).rawPeakIdx         = spindle_det(ch).rawPeakIdx( spins_good_lengths )';         %#ok
        spindle_det(ch).rawTroughIdxAmp    = spindle_det(ch).rawTroughIdxAmp( spins_good_lengths )';    %#ok
        spindle_det(ch).rawTroughIdx       = spindle_det(ch).rawTroughIdx( spins_good_lengths )';       %#ok

        spindle_det(ch).envPeakIdxAmp      = spindle_det(ch).envPeakIdxAmp( spins_good_lengths )';      %#ok
        spindle_det(ch).envPeakIdx         = spindle_det(ch).envPeakIdx( spins_good_lengths )';         %#ok
        spindle_det(ch).envPeakIdxTraceVal = spindle_det(ch).envPeakIdxTraceVal( spins_good_lengths )'; %#ok

        spindle_det(ch).spindlesFrSlope    = spindle_det(ch).spindlesFrSlope( spins_good_lengths )'; %#ok    
        spindle_det(ch).spindlesWaveletFreq= spindle_det(ch).spindlesWaveletFreq( spins_good_lengths )'; %#ok  
        
        spindle_det(ch).isi_previous_event_in_millis = spindle_det(ch).isi_previous_event_in_millis';%#ok
        spindle_det(ch).isi_next_event_in_millis     = spindle_det(ch).isi_next_event_in_millis';%#ok
        spindle_det(ch).isi_min_in_millis            = spindle_det(ch).isi_min_in_millis';%#ok
        spindle_det(ch).isi_max_in_millis            = spindle_det(ch).isi_max_in_millis';%#ok
        spindle_det(ch).isi_mean_in_millis           = spindle_det(ch).isi_mean_in_millis';%#ok
        
        %-----------
        %===============================

        spindle_det(ch).startSample        = adj_start_idxs( spins_good_lengths )';       %#ok
        spindle_det(ch).endSample          = adj_end_idxs( spins_good_lengths )';         %#ok  % whos spins_good_lengths peak_idxs start_idxs end_idxs durations
        spindle_det(ch).peakSample         = peak_idxs( spins_good_lengths );             %#ok
        spindle_det(ch).peakAmp            = current_data( spindle_det(ch).peakSample );  %#ok
        spindle_det(ch).symmetry           = (spindle_det(ch).peakSample - spindle_det(ch).startSample) ./ (spindle_det(ch).endSample - spindle_det(ch).startSample); %#ok

        %-----------
        %===============================

        spindle_det(ch).freqFFT            = spindle_det(ch).freqFFT( spins_good_lengths )';  %#ok
        spindle_det(ch).duration           = spindle_det(ch).duration( spins_good_lengths )'; %#ok

        %-----------
        %===============================
                
        spindle_det(ch).minDuration        = minmaxdurs(1);                %#ok
        spindle_det(ch).maxDuration        = minmaxdurs(2);                %#ok
        if exist('fc','var')
         spindle_det(ch).midFreq           = fc;                           %#ok
         spindle_det(ch).fc                = fc;                           %#ok
        else
         spindle_det(ch).midFreq           = round( ((bpparms(2)-bpparms(1))/2+bpparms(1)) *10)/10;  %#ok  
        end
        if exist('fb','var')
         spindle_det(ch).fb                = fb;                           %#ok
        else
         spindle_det(ch).fb                = NaN;                          %#ok
         spindle_det(ch).lowerFreq         = bpparms(1);                  %#ok
         spindle_det(ch).upperFreq         = bpparms(2);                  %#ok
        end
        spindle_det(ch).fs                 = fs;                           %#ok
        spindle_det(ch).meanSigma          = signalmean(ch);               %#ok
        spindle_det(ch).sdSigma            = signalstd(ch);                %#ok
        spindle_det(ch).medianSigma        = signalmedian(ch);             %#ok         
        spindle_det(ch).amplFactor         = ampl_factor;                  %#ok
        spindle_det(ch).thresholdType      = 'median';                     %#ok
        spindle_det(ch).threshold          = ampl_factor*signalmedian(ch); %#ok
        spindle_det(ch).number             = sum(spins_good_lengths);      %#ok
        spindle_det(ch).label              = ch_names{ch};                 %#ok      

        %===============================

        spindle_det(ch).density            = spindle_det(ch).number / ((size(Data,1)/fs)/60); %#ok
        spindle_det(ch).datalen            = size(Data,1);                 %#ok
        
        %===============================

        clear durations end_idxs peak_idxs start_idxs spins_good_lengths adj_end_idxs adj_start_idxs locs over_idxs over

        filtTraces(ch).filtereddata        = Data_cplx(:,ch);              %#ok
        filtTraces(ch).envelope            = sqrt( Data_sqenvff(:,ch) );   %#ok %NOTE THAT this was the SQUARED envelope, not the original signal envelope from abs(hilbert(signal))
        filtTraces(ch).raw                 = Data(:,ch);                   %#ok
                
                
            %
            
        
    end % End the loop for each channel

warning('on','all'); warning
  
end
  
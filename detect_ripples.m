function [ripple_det] = detect_ripples( Data,...
                                        ch_names,...
                                        Fs,...
                                        ncyclesthresh,...
                                        ampl_factor,...
                                        filt_params,...
                                        verify,...
                                        manualthreshs,...
                                        noiserej_bands,...
                                        artratio)                                        
%% Function to detect hippocampal ripples 
% 
% -------------------------------------------------------------------------
% This function detects hippocampal ripples, returning a struct with ripple
% time points, features, spread, etc. I is based on the validated detector 
% from Staba et al. 2002, with additional spike-rejection steps from 
% Helfrich et al. 2019.
% ------------------------------------------------------------------------- 
%
% Dependencies
% ---------------------------
% functions:
%
% ---------------------------
%
% Input 
% ---------------------------
%     Data               :    array of input data (TIME x CHANNELS)
%     ch_names           :    string cells with channel names  
%     Fs                 :    Sampling rate in Hz (scalar)
%     ncyclesthresh      :    # of rectified signal peaks needed 
%                             for ripple candidacy (default is '6' from 
%                             Staba et al. 2002)
%     ampl_factor        :    scalar; threshold for signal envelope to 
%                             cross for candidate detection (default is '5'
%                             from Staba et al. 2002). This is applied as 
%                             a multiplier to channel-specific standard
%                             deviations. Consider replacing with median
%                             absolute deviations as a more appropriate
%                             statistical measure of spread, but stick to
%                             standard deviations if you want to cite the
%                             Staba paper.
%     filt_params        :    struct of parameters for band pass filtering.
%                             Example:
%                                 params.lowstop          = 70; 
%                                 params.lowpassfreq      = 80;     
%                                 params.highpassfreq     = 250;
%                                 params.highstop         = 280;
%     verify             :    legacy input, for visualizing detections.
%                             Logical scalar, must be either '0' or '1'
%     manualthreshs      :    ignore if empty set ('[]'); otherwise will
%                             apply this scalar value as the envelope 
%                             threshold typically derived in a channel-
%                             specific fashion using the multiplier 
%                             specified in 'ampl_factor'
%     noiserej_bands     :    for additional spike rejection. can be any of the
%                             following strings:
%                                   'none'   :    no ratio-based detection performed
%                                   'high'   :    a high band (160-190Hz), defined in the
%                                                 detect_ripples function, is used for broadband
%                                                 artifact rejection
%                                   'low'    :    a low band (50-60 Hz), defined in the detect_ripples
%                                                 function, is used for broadband artifact rejection
%                               'highlow'    :    both bands used
%                             I recommend using 'none', since it wasn't 
%                             done in the Staba paper and in my experience 
%                             didn't work very well for artifact rejection 
%                             -- it threw out a lot of real detections along
%                             with spikes.
%     artratio           :    scalar. If positive, this is the ratio of ripple power to
%                             artifact band power, defined in the 'artrejmode' variable above, that a
%                             detection must exceed to be counted as a ripple.  So for example if this
%                             was '3', then ripple power needs to be at least 3x greater than artifact
%                             band power to be counted as a detection. Setting this to '-1', as I did
%                             here, skips ratio-based spike rejection and instead applies the "3-peak
%                             rule" implemented in the Helfrich et al. 2019 paper (i.e., at the peak of
%                             a ripple candidate, the RAW trace must exhibit a ripple band oscillatory
%                             peak on either side of the detected peak, within a 40ms window -- so it
%                             is looking to find 3 peaks within 40ms centered on the peak of the ripple
%
%     EXAMPLE FUNCTION CALL:
%           
%           [ripples] = detect_ripples( ecog_signals, channel_labels, 512, 6, 5, params, false, [], 'none', -1);
%
% ---------------------------
%
% Output
% ---------------------------
%     ripple_det         :    struct with ripple info
%
% ---------------------------
%
% Author(s):   William Coon   wcoon@mgh.harvard.edu
%
% Last edited: September 4, 2019 1304hrs
%
% -------------------------------------------------------------------------
%%                                        

warning('off','all'); warning
fprintf('============================================================\n')
fprintf('============================================================\n')
fprintf('============================================================\n')
fprintf('This version was last edited September 4th, 2019\n')    
fprintf('============================================================\n')
fprintf('============================================================\n')
fprintf('============================================================\n')
pause(3)

%% 
    % ensure data dimensions correct
    dims = size(Data);
    if dims(2) > dims(1)
        Data = Data';
    end
    
    %get/set 
    if ~exist('manualthreshs','var')
        manualthreshs = [];
    end    
    if ~exist('noiserej_bands','var')
        fprintf('No specified noise rejection parameter, defaulting to high band only...\n')
        noiserej_bands = 'high';
    end
    if ~( strcmp(noiserej_bands,'high') || strcmp(noiserej_bands,'highlow') || strcmp(noiserej_bands,'low') || strcmp(noiserej_bands,'none') )
        error('Last input parameter, ''noiserej'', must be either ''high'' or ''highlow'' or ''low'' or ''none''')
    end
    if ~exist('artratio','var')
        error('No ratio set for artifact rejection.')
    end
    
    
    %% Calculate the (bandpass) for each channel

    %do bandpassing and hilbert instead
    Data_bp = nan(size(Data));
    for x=1:size(Data,2)
        EEGWave       = firbandpass2(Fs, double( Data(:,x) ), filt_params.highpassfreq, filt_params.lowpassfreq, filt_params.highstop, filt_params.lowstop);
        Data_bp(:,x)  = real(EEGWave); clear EEGWave
    end
%     Data_bp = Data_bp';
    
    
    %% Envelope

    Data_env    = abs( hilbert(Data_bp) ); % take the moving average using the above window


    %% Compute Threshhold (Separately for Each Channel)
%     signalmean = mean(Data_env); % compute mean amplitude of rectified signal
    signalstd  = std(Data_env); % compute std dev of the amplitude of rectified signal
    if isempty(manualthreshs)
        threshold  = signalstd.*ampl_factor; % % defines the threshold
        fprintf('Using auto thresholds...\n')                   
    else
        threshold  = manualthreshs;
        fprintf('Using manually set thresholds...\n')
    end
   
    for ch=1:size(Data,2) % Loop for each signal
        
        fprintf('Working on Channel %s.\n',ch_names{ch});
        
        %% Detect ripples
        
        current_data = Data_env(:,ch);
        rectif       = abs(Data_bp(:,ch));                        % rectified, band passed signal
        fs=Fs;
        
        if verify
            vidx = 1801+fs*16;
            close all; figure; set(gcf,'position',[49         479        1826         567]); setPaperSize(gcf); fs = Fs;
            clf
            plot(Data_bp(vidx:vidx+Fs,ch)); axis tight
            hold on; plot(Data_env(vidx:vidx+Fs,ch),'r')
            plot(abs(Data_bp(vidx:vidx+Fs,ch)),'g')
        end
        
        if any(isnan(current_data))
            ripple_det(ch).bads = 1;
        else
            ripple_det(ch).bads = 0;
        end
        
        
        %===============================  
        % Mark all points over threshold as '1'
        over = current_data>threshold(1,ch);                
                if verify, overtmp = over(vidx:vidx+fs);
                hold on; plot(overtmp*40,'k'), end
        %-------------------------------     
        
        
        %===============================    
        % rectified signal peaks (need at least six)    
        [~,locs] = findpeaks(rectif);                           
                if verify, [~,locstmp] = findpeaks(rectif(vidx:vidx+Fs));                
                hold on; stem(locstmp,30*ones(length(locstmp),1),'o'), end
        %-------------------------------
        
               
        %===============================      
        over_idxs  = find(over>0);                          % sample indexes where envelope exceeds threshold
        start_idxs = find(diff(over)>0)+1;                  % beginning indexes of each segment of over_idxs
        if over(1)>0, start_idxs = [1; start_idxs]; end     % in case first sample is over threshold
        end_idxs   = find(diff(over)<0)+1;                  % end indexes of each segment of over_idxs
        if length(end_idxs)<length(start_idxs), end_idxs = [end_idxs; length(current_data)]; end %#ok  if end of last cand is clipped, make end of signal the last sample
        %-------------------------------   
        
        
        % ===== LESS THAN 'n' REJECTION==========
        % ------------------------ 
        %remove case where less than 'n' oscillations happen
        for ix = 1:length(start_idxs)  
            clc
            fprintf('total = %d  (Less than ''n''-cycle rejection)\n', length(start_idxs))
            fprintf('%d  ',ix)
            num_peaks_in_this_candidate = length( intersect(start_idxs(ix):end_idxs(ix),locs) );
            if num_peaks_in_this_candidate < ncyclesthresh
                over( start_idxs(ix):end_idxs(ix) ) = zeros( length(start_idxs(ix):end_idxs(ix)), 1 );
            end                                                          
        end; fprintf('\n')       
                if verify, overtmp = over(vidx:vidx+fs);
                hold on; plot(overtmp*45,'y'), end
        %-------------------------------     
        
        
        %===============================  
        %re-tally
        start_idxs = find(diff(over)>0);                    % beginning indexes of each segment of over_idxs
        if over(1)>0, start_idxs = [1; start_idxs]; end     % in case first sample is over threshold
        end_idxs   = find(diff(over)<0);                    % end indexes of each segment of over_idxs
        if length(end_idxs)<length(start_idxs), end_idxs = [end_idxs; length(current_data)]; end %#ok  if end of last cand is clipped, make end of signal the last sample
        %------------------------------- 
                           
        
        % ===== SPIKE REJECTION==========
        %===============================    
        % raw signal peaks (need at least 3 within +/- 20ms from time zero...    
        [~,locs2] = findpeaks(Data(:,ch));                           
        %-------------------------------
        % ------------------------ 
        %remove case where less than 3 spikes occur within +/-20ms from ripple peak           
        spike_count = 0;
      if artratio<0  
        for ix = 1:length(start_idxs)
            clc
            fprintf('total = %d  (Spike rejection...)\n', length(start_idxs))
            fprintf('%d  ',ix)
            highest_peak_in_this_candidate_ampl = max( current_data( start_idxs(ix):end_idxs(ix) ) );
            highest_peak_in_this_candidate_idx  = start_idxs(ix) + find( current_data( start_idxs(ix):end_idxs(ix) ) == highest_peak_in_this_candidate_ampl, 1, 'first' );
            %40-ms window, for rejecting spikes instead of true oscillations
            win = highest_peak_in_this_candidate_idx-round(Fs*(40/1000))/2:highest_peak_in_this_candidate_idx+round(Fs*(40/1000))/2;
            num_raw_trace_peaks = length( intersect(win,locs2) );
            if num_raw_trace_peaks < 3
                over( start_idxs(ix):end_idxs(ix) ) = zeros( length(start_idxs(ix):end_idxs(ix)), 1 );
                spike_count = spike_count+1;
            end  
        end
            
        %===============================  
        %re-tally
        start_idxs = find(diff(over)>0);                    % beginning indexes of each segment of over_idxs
        if over(1)>0, start_idxs = [1; start_idxs]; end     % in case first sample is over threshold
        end_idxs   = find(diff(over)<0);                    % end indexes of each segment of over_idxs
        if length(end_idxs)<length(start_idxs), end_idxs = [end_idxs; length(current_data)]; end %#ok  if end of last cand is clipped, make end of signal the last sample
        %------------------------------- 
      end
        
        % ===== FFT for ARTIFACT REJECTION==========
        % ------------------- 
                fprintf('Rejecting events with broadband increase (160-190Hz) greater than (80-150Hz)...\n') 
                fft_artrejratio_count = 0; fft_artrejabs_count = 0;
                for ix = 1:length(start_idxs) 
        % ------------------- 

                    % --------spectrum of RAW data (for getting non-wavelet freq and artefact rejection)

                    spiRaw              = real(Data(start_idxs(ix):end_idxs(ix),ch)); % data(spiStartIndex:spiEndIndex);
                    ns                  = length(spiRaw);                                     % number of samples analyzed
                    fs=Fs;
                    
                    %length for FFT
                    fftLength           = 2^nextpow2((fs*500/1000));

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
                    lowFrInd            = find(spindleFreqs>=80,1,'first');
                    highFrInd           = find(spindleFreqs<=150,1,'last');

                    if highFrInd < lowFrInd
                        highFrInd       = lowFrInd;
                    end

                    %artifact frequencies
                    %broadband frequency range that should not show power increase during
                    %ripples
                    lofreq_noise        = 50;%160;
                    hifreq_noise        = 60;%190;
                    lofreq_noise1       = 160;
                    hifreq_noise1       = 190;
                    low_art             = find(spindleFreqs>(lofreq_noise),1,'first');
                    high_art            = find(spindleFreqs<(hifreq_noise),1,'last');
                    low_art1            = find(spindleFreqs>(lofreq_noise1),1,'first');
                    high_art1           = find(spindleFreqs<(hifreq_noise1),1,'last');

                    %power at all frequencies in broadband range
                    if strcmp(noiserej_bands,'highlow')
                        artpow              = (spiFFTPow([low_art:high_art-1 low_art1:high_art1-1])); 
                    elseif strcmp(noiserej_bands,'high')
                        artpow              = (spiFFTPow(low_art1:high_art1-1)); 
                    elseif strcmp(noiserej_bands,'low')
                        artpow              = (spiFFTPow(low_art:high_art-1)); 
                    elseif strcmp(noiserej_bands,'none')
                        artpow              = [];     
                    end

                    [max_spipow, spiMaxPowInd]  = max(spiFFTPow(lowFrInd:highFrInd));
                     ripple_det(ch).freqFFT(ix) = spindleFreqs(lowFrInd+spiMaxPowInd-1);                    

                    %reject any power increases absolutely greater than ripple band 
                    if ~strcmp(noiserej_bands,'none')
                        if any(artpow>max_spipow)
                            over( start_idxs(ix):end_idxs(ix) ) = zeros( length(start_idxs(ix):end_idxs(ix)), 1 );
                            fft_artrejabs_count = fft_artrejabs_count + 1;
                        end                    

                        %reject any where ripple band increase isn't on average X-times higher than the artifact band
                        if max_spipow/max(artpow) <= artratio
                            over( start_idxs(ix):end_idxs(ix) ) = zeros( length(start_idxs(ix):end_idxs(ix)), 1 );
                            fft_artrejratio_count = fft_artrejratio_count + 1;
                        end
                    end

        % -------------------        
                end %end FFT  for ARTIFACT REJECTION (for each spindle event loop (for ix = 1:length(start_idxs))
        % -------------------    
        
        
        %===============================  
        % -------------------
        %remove cases where artifact is present (ripple band pass filtered trace exceeds 50uV)
        artthresh = 50;
        fprintf('Removing artifact-laden detections (greater than %duV in filtered signal)...\n',artthresh); pause(1);
        for ix = 1:length(start_idxs)  
            clc
            fprintf('total = %d  (Amp threshold rejection)\n', length(start_idxs))
            fprintf('%d  ',ix)
            sig_max = max(abs(Data_bp(start_idxs(ix):end_idxs(ix), ch )));
            if sig_max > artthresh
                over( start_idxs(ix):end_idxs(ix) ) = zeros( length(start_idxs(ix):end_idxs(ix)), 1 );
            end          
        end; fprintf('\n')                    
        % -------------------        
        % end Amplitude Threshold
        % -------------------    

        
        %===============================  
        %re-tally
        start_idxs = find(diff(over)>0);                    % beginning indexes of each segment of over_idxs
        if over(1)>0, start_idxs = [1; start_idxs]; end     % in case first sample is over threshold
        end_idxs   = find(diff(over)<0);                    % end indexes of each segment of over_idxs
        if length(end_idxs)<length(start_idxs), end_idxs = [end_idxs; length(current_data)]; end %#ok  if end of last cand is clipped, make end of signal the last sample
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
        %-------------------------------
        
            %===============================    
            % ------- accounts for ripple that starts right near the end
            % of the record and may not "end" before the data file
            % ends

if ~isempty(peak_idxs)
                if peak_idxs(end) > length(current_data), peak_idxs(end) = length(current_data); end

            %-------------------------------                
                
        
        %===============================    
        % -------- DURATIONS 
        % NOTE: adj_start/end_ixs are JUST for FWHM duration definitions.
        % The Staba et al. algorithm should use the same definitions as
        % above, I think
        ix = 1; buffersize = 2; fs=Fs; clear adj_start_idxs adj_end_idxs
        for ix = 1:length(start_idxs)    

            if peak_idxs(ix)-buffersize*fs < 1, ix1 = 1; else ix1 = peak_idxs(ix)-buffersize*fs; end %#ok
            if peak_idxs(ix)+buffersize*fs > length(current_data), ix2 = length(current_data); else ix2 = peak_idxs(ix)+buffersize*fs; end

            centerindex = peak_idxs(ix); % begin search for edges at predetermined peak
            dropmax     = sqrt(current_data(centerindex))/2;

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

            %store DURATION of ripple defined by Full-Width Half-Max (FWHM thresholding)
            ripple_det(ch).FWHMduration(ix) = (1000/fs)*(adj_end_idxs(ix) - adj_start_idxs(ix));

        end                
        %------------------------------
        
                
        %===============================    
        % ------- Merge any overlapping ripples...
            fprintf('Merging overlapping ripple events...\n'); pause(1)                   
            i = length(start_idxs);
            origlen = i;
            while i ~= 1 
                if ~isempty(intersect( (start_idxs(i):end_idxs(i)), (start_idxs(i-1):end_idxs(i-1)) ))
%                     adj_end_idxs(i-1)   = adj_end_idxs(i);
                    end_idxs(i-1)     = end_idxs(i);
                    start_idxs(i)     = []; 
                    end_idxs(i)       = []; 
                    peak_idxs(i)      = [];
                    adj_start_idxs(i) = [];
                    adj_end_idxs(i)   = [];
                    ripple_det(ch).FWHMduration(i) = [];
                    ripple_det(ch).freqFFT = [];
                end
                i = i - 1;
            end
            fprintf('Number of detections reduced from %d to %d.\n',origlen,length(start_idxs)); pause(2)
        %===============================
        
        
        % ===== FFT for CHARACTERIZATION ==========
        % ------------------- 
                fprintf('Finding peak frequency (FFT)...\n')
                ripple_det(ch).freqFFT = zeros( length(start_idxs), 1 ); 
                for ix = 1:length(start_idxs) 
        % ------------------- 

                    % --------spectrum of RAW data (for getting non-wavelet freq and artefact rejection)

                    spiRaw              = real(Data(start_idxs(ix):end_idxs(ix),ch)); % data(spiStartIndex:spiEndIndex);
                    ns                  = length(spiRaw);                                     % number of samples analyzed
                    
                    %length for FFT
                    fftLength           = 2^nextpow2((fs*500/1000));

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
                    lowFrInd            = find(spindleFreqs>=80,1,'first');
                    highFrInd           = find(spindleFreqs<=150,1,'last');

                    if highFrInd < lowFrInd
                        highFrInd       = lowFrInd;
                    end
                    
                    [~,spiMaxPowInd]   = max(spiFFTPow(lowFrInd:highFrInd));
                    ripple_det(ch).freqFFT(ix) = spindleFreqs(lowFrInd+spiMaxPowInd-1);                    

        % -------------------        
                end %end FFT (for each spindle event loop (for ix = 1:length(start_idxs))
        % ------------------- 
        
        
        %==========================================
        %====== CHARACTERIZATION DETAILS ==========
        %----------------------------------------
            ctr = 0;
            for ix = 1:length(start_idxs)  
                clc
                fprintf('total: = %d\n', length(start_idxs))
                fprintf('%d  ',ix)
                ctr=ctr+1;
                ripple_det(ch).num_peaks(ctr) = length( intersect(start_idxs(ix):end_idxs(ix),locs) ) / 2;

                % Peaks and Troughs    
                ripple_det(ch).filtPeakIdxAmp(ix)   = max( real(Data_bp(start_idxs(ix):end_idxs(ix),ch)) );
                ripple_det(ch).filtPeakIdx(ix)      = start_idxs(ix)-1 + find( real(Data_bp(start_idxs(ix):end_idxs(ix),ch))==max( real(Data_bp(start_idxs(ix):end_idxs(ix),ch)) ), 1, 'first'  );    
                ripple_det(ch).filtTroughIdxAmp(ix) = min( real(Data_bp(start_idxs(ix):end_idxs(ix),ch)) );
                ripple_det(ch).filtTroughIdx(ix)    = start_idxs(ix)-1 + find( real(Data_bp(start_idxs(ix):end_idxs(ix),ch))==min( real(Data_bp(start_idxs(ix):end_idxs(ix),ch)) ), 1, 'first'  );    

                ripple_det(ch).rawPeakIdxAmp(ix)    = max( real(Data(start_idxs(ix):end_idxs(ix),ch)) );
                ripple_det(ch).rawPeakIdx(ix)       = start_idxs(ix)-1 + find( real(Data(start_idxs(ix):end_idxs(ix),ch))==max( real(Data(start_idxs(ix):end_idxs(ix),ch)) ), 1, 'first' );    
                ripple_det(ch).rawTroughIdxAmp(ix)  = min( real(Data(start_idxs(ix):end_idxs(ix),ch)) );
                ripple_det(ch).rawTroughIdx(ix)     = start_idxs(ix)-1 + find( real(Data(start_idxs(ix):end_idxs(ix),ch))==min( real(Data(start_idxs(ix):end_idxs(ix),ch)) ), 1, 'first' );    

                ripple_det(ch).envPeakIdx(ix)       = start_idxs(ix)-1 + find( real(Data_env(start_idxs(ix):end_idxs(ix),ch))==max( real(Data_env(start_idxs(ix):end_idxs(ix),ch)) ), 1, 'first' );    
                ripple_det(ch).envPeakIdxAmp(ix)    = max( real(Data_env(start_idxs(ix):end_idxs(ix),ch)) );    
                ripple_det(ch).envPeakIdxTraceVal(ix) = real(Data(ripple_det(ch).envPeakIdx(ix),ch));
                
                %Inter-event intervals
                if length(ripple_det(ch).filtPeakIdxAmp) ~= 1                    
                    if ix ~= 1 && ix ~= length(start_idxs)
                        ripple_det(ch).isi_previous_event_in_millis(ix) = (start_idxs(ix) - start_idxs(ix-1)) / (fs/1000);
                        ripple_det(ch).isi_next_event_in_millis(ix)     = (start_idxs(ix+1) - start_idxs(ix)) / (fs/1000);
                        ripple_det(ch).isi_min_in_millis(ix)            = min(  [ripple_det(ch).isi_previous_event_in_millis(ix) ripple_det(ch).isi_next_event_in_millis(ix)] );
                        ripple_det(ch).isi_max_in_millis(ix)            = max(  [ripple_det(ch).isi_previous_event_in_millis(ix) ripple_det(ch).isi_next_event_in_millis(ix)] );
                        ripple_det(ch).isi_mean_in_millis(ix)           = mean( [ripple_det(ch).isi_previous_event_in_millis(ix) ripple_det(ch).isi_next_event_in_millis(ix)] );
                    elseif ix == 1
                        ripple_det(ch).isi_previous_event_in_millis(ix) = NaN;
                        ripple_det(ch).isi_next_event_in_millis(ix)     = (start_idxs(ix+1) - start_idxs(ix)) / (fs/1000);
                        ripple_det(ch).isi_min_in_millis(ix)            = ripple_det(ch).isi_next_event_in_millis(ix);
                        ripple_det(ch).isi_max_in_millis(ix)            = ripple_det(ch).isi_next_event_in_millis(ix);
                        ripple_det(ch).isi_mean_in_millis(ix)           = ripple_det(ch).isi_next_event_in_millis(ix);
                    elseif ix == length(start_idxs)
                        ripple_det(ch).isi_previous_event_in_millis(ix) = (start_idxs(ix) - start_idxs(ix-1)) / (fs/1000);
                        ripple_det(ch).isi_next_event_in_millis(ix)     = NaN;
                        ripple_det(ch).isi_min_in_millis(ix)            = ripple_det(ch).isi_previous_event_in_millis(ix);
                        ripple_det(ch).isi_max_in_millis(ix)            = ripple_det(ch).isi_previous_event_in_millis(ix);
                        ripple_det(ch).isi_mean_in_millis(ix)           = ripple_det(ch).isi_previous_event_in_millis(ix);
                    end
                else
                    ripple_det(ch).isi_previous_event_in_millis(ix) = nan;
                    ripple_det(ch).isi_next_event_in_millis(ix)     = nan;
                    ripple_det(ch).isi_min_in_millis(ix)            = nan;
                    ripple_det(ch).isi_max_in_millis(ix)            = nan;
                    ripple_det(ch).isi_mean_in_millis(ix)           = nan;
                end
                    
            end

            if isfield(ripple_det(ch),'num_peaks')
                ripple_det(ch).num_peaks = ripple_det(ch).num_peaks;
            else
                ripple_det(ch).num_peaks = 0;
            end
            spin = start_idxs;

            ripple_det(ch).label         = ch_names{ch};
            ripple_det(ch).sample        = spin;
            ripple_det(ch).ripple_count = length(spin);
            ripple_det(ch).backgr_std    = signalstd(ch);
            ripple_det(ch).datalen       = size(Data,1);
            ripple_det(ch).startSample   = spin;
            ripple_det(ch).endSample     = end_idxs;
            ripple_det(ch).FWHMstartSample   = adj_start_idxs';
            ripple_det(ch).FWHMendSample     = adj_end_idxs';
            ripple_det(ch).duration      = (ripple_det(ch).endSample-ripple_det(ch).startSample).*(1/Fs)';
            ripple_det(ch).number        = length(spin);
            ripple_det(ch).num_peaks     = ripple_det(ch).num_peaks';
            ripple_det(ch).FWHMduration  = ripple_det(ch).FWHMduration';
            ripple_det(ch).num_spikes_rejected = spike_count;
            ripple_det(ch).num_abs_art_rejected = fft_artrejabs_count;
            ripple_det(ch).num_ratio_art_rejected = fft_artrejratio_count;
        
        %-----------
        %===============================

            % Peaks and Troughs    
            ripple_det(ch).filtPeakIdxAmp   = ripple_det(ch).filtPeakIdxAmp';
            ripple_det(ch).filtPeakIdx      = ripple_det(ch).filtPeakIdx';    
            ripple_det(ch).filtTroughIdxAmp = ripple_det(ch).filtTroughIdxAmp';
            ripple_det(ch).filtTroughIdx    = ripple_det(ch).filtTroughIdx';    

            ripple_det(ch).rawPeakIdxAmp    = ripple_det(ch).rawPeakIdxAmp';
            ripple_det(ch).rawPeakIdx       = ripple_det(ch).rawPeakIdx';    
            ripple_det(ch).rawTroughIdxAmp  = ripple_det(ch).rawTroughIdxAmp';
            ripple_det(ch).rawTroughIdx     = ripple_det(ch).rawTroughIdx';    

            ripple_det(ch).envPeakIdx       = ripple_det(ch).envPeakIdx';    
            ripple_det(ch).envPeakIdxAmp    = ripple_det(ch).envPeakIdxAmp';    
            ripple_det(ch).envPeakIdxTraceVal = ripple_det(ch).envPeakIdxTraceVal';
            ripple_det(ch).freqFFT          = ripple_det(ch).freqFFT;
            
            ripple_det(ch).isi_previous_event_in_millis = ripple_det(ch).isi_previous_event_in_millis';
            ripple_det(ch).isi_next_event_in_millis     = ripple_det(ch).isi_next_event_in_millis';
            ripple_det(ch).isi_min_in_millis            = ripple_det(ch).isi_min_in_millis';
            ripple_det(ch).isi_max_in_millis            = ripple_det(ch).isi_max_in_millis';
            ripple_det(ch).isi_mean_in_millis           = ripple_det(ch).isi_mean_in_millis';

        %-----------    
        %===============================
        
else %if isempty (no detections), fill with empty sets as placeholders
    
    ripple_det(ch).filtPeakIdxAmp   = [];
    ripple_det(ch).filtPeakIdx      = [];
    ripple_det(ch).filtTroughIdxAmp = [];
    ripple_det(ch).filtTroughIdx    = [];

    ripple_det(ch).rawPeakIdxAmp    = [];
    ripple_det(ch).rawPeakIdx       = [];
    ripple_det(ch).rawTroughIdxAmp  = [];
    ripple_det(ch).rawTroughIdx     = [];

    ripple_det(ch).envPeakIdx       = [];
    ripple_det(ch).envPeakIdxAmp    = [];
    ripple_det(ch).envPeakIdxTraceVal = [];
    ripple_det(ch).freqFFT          = [];
    
    ripple_det(ch).num_peaks        = [];
    
    ripple_det(ch).label            = ch_names{ch};
    ripple_det(ch).sample           = [];
    ripple_det(ch).ripple_count     = [];
    ripple_det(ch).backgr_std       = [];
    ripple_det(ch).datalen          = size(Data,1);
    ripple_det(ch).startSample      = [];
    ripple_det(ch).endSample        = [];
    ripple_det(ch).FWHMstartSample  = [];
    ripple_det(ch).FWHMendSample    = [];
    ripple_det(ch).duration         = [];
    ripple_det(ch).number           = 0;
    ripple_det(ch).num_peaks        = [];
    ripple_det(ch).FWHMduration     = [];
    
    ripple_det(ch).isi_previous_event_in_millis = [];
    ripple_det(ch).isi_next_event_in_millis     = [];
    ripple_det(ch).isi_min_in_millis            = [];
    ripple_det(ch).isi_max_in_millis            = [];
    ripple_det(ch).isi_mean_in_millis           = [];
    ripple_det(ch).num_spikes_rejected          = [];
    ripple_det(ch).num_abs_art_rejected = [];
    ripple_det(ch).num_ratio_art_rejected = [];
               
end % end if ~isempty

            
    end  % End the loop for each channel
    
    fprintf('\ndone.\n')
    warning('on','all'); warning

    
    
    function [signalbp] = firbandpass2(fs, signal, highpassfreq, lowpassfreq, highstop, lowstop)
    %% FIR band pass
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
    %     Example = firbandpass2( 100, ecog_data, 10, 15, 8, 17 );
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


    end %end firbandpass2


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


end % end function

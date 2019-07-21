function [spindle_det] = detect_ripples(Data,...
                                                            ch_names,...
                                                            Fs,...
                                                            ncyclesthresh,...
                                                            ampl_factor,...
                                                            filt_params,...
                                                            verify)
%========================================================================
%% RIPPLE detection
%========================================================================
%---------------------------
% 
% ncyclesthresh       = 6;
% 
% % Factor to multiply the mean of the signal in the spindle detection
% rip_ampl_factor = 5.0; % default from Staba et al. is 5
% 
% if fs<257
%     params.highpassfreq = 90;
%     params.lowpassfreq  = 120;
%     params.highstop     = 80;
%     params.lowstop      = 127;
% else
%     params.highpassfreq = 90;
%     params.lowpassfreq  = 250;
%     params.highstop     = 80;
%     params.lowstop      = 280;
% end
% verify = true;
% % movingwinsize       = 0.003; %seconds, width of moving average window for RMS filtering
% 
% [allRipples] = fun3_ripple_detection_wav_altered( double(n2ecog)',...
%                                                     srtd_hipp_ch_strs,...
%                                                     fs,...
%                                                     ncyclesthresh,...
%                                                     rip_ampl_factor,...
%                                                     params,...
%                                                     ~verify );

%========================================================            
                                                        
fprintf('============================================================\n')
fprintf('============================================================\n')
fprintf('============================================================\n')

fprintf('This version was last edited Feb. 6th, 2019\n')    

fprintf('============================================================\n')
fprintf('============================================================\n')
fprintf('============================================================\n')                                         
%% 
    % ensure data dimensions correct
    dims = size(Data);
    if dims(2) > dims(1)
        Data = Data';
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
    threshold  = signalstd.*ampl_factor; % % defines the threshold
   
    for ch=1:size(Data,2) % Loop for each signal
        
        fprintf('Working on Channel %s.\n',ch_names{ch});
        
        %% Detect spindles
        
        current_data = Data_env(:,ch);
        rectif       = abs(Data_bp(:,ch));                        % rectified, band passed signal

        if verify
            vidx = 1801+fs*16;
            close all; figure; set(gcf,'position',[49         479        1826         567]); setPaperSize(gcf); fs = Fs;
            clf
            plot(Data_bp(vidx:vidx+Fs,ch)); axis tight
            hold on; plot(Data_env(vidx:vidx+Fs,ch),'r')
            plot(abs(Data_bp(vidx:vidx+Fs,ch)),'g')
        end
        
        if any(isnan(current_data))
            spindle_det(ch).bads = 1;
        else
            spindle_det(ch).bads = 0;
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
                if verify, [~,locstmp] = findpeaks(rectif(vidx:vidx+fs));                
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
            fprintf('total = %d  (Less than ''n''-cycle rejection\n', length(start_idxs))
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
        
        
        % ===== FFT for ARTIFACT REJECTION==========
        % ------------------- 
                fprintf('Rejecting events with broadband increase (160-190Hz) greater than (80-150Hz)...\n') 
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
                    %spindles
                    lofreq_noise        = 160;
                    hifreq_noise        = 190;
                    low_art             = find(spindleFreqs>(lofreq_noise),1,'first');
                    high_art            = find(spindleFreqs<(hifreq_noise),1,'last');

                    %power at all frequencies in broadband range
                    artpow              = (spiFFTPow(low_art:high_art-1));                    

                    [max_spipow, spiMaxPowInd]  = max(spiFFTPow(lowFrInd:highFrInd));
                     spindle_det(ch).freqFFT(ix) = spindleFreqs(lowFrInd+spiMaxPowInd-1);                    

                    if any(artpow>max_spipow)
                    	over( start_idxs(ix):end_idxs(ix) ) = zeros( length(start_idxs(ix):end_idxs(ix)), 1 );
                    end

        % -------------------        
                end %end FFT  for ARTIFACT REJECTION (for each spindle event loop (for ix = 1:length(start_idxs))
        % -------------------    
        
        
        %===============================  
        % -------------------
        %remove cases where artifact is present (raw trace exceeds 50uV)
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
            % ------- accounts for spindle that starts right near the end
            % of the record and may not "end" before the data file
            % ends
            if peak_idxs(end) > length(current_data), peak_idxs(end) = length(current_data); end
            %-------------------------------                
                
        
        %===============================    
        % -------- DURATIONS (operates on 100ms averaged signal envelope, NOT the squared envelope but the envelope from the hilbert directly)
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

            %store DURATION of spindle defined by Full-Width Half-Max (FWHM thresholding)
            spindle_det(ch).FWHMduration(ix) = (1000/fs)*(adj_end_idxs(ix) - adj_start_idxs(ix));

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
                    spindle_det(ch).FWHMduration(i) = [];
                    spindle_det(ch).freqFFT = [];
                end
                i = i - 1;
            end
            fprintf('Number of detections reduced from %d to %d.\n',origlen,length(start_idxs)); pause(2)
        %===============================
        
        
        % ===== FFT for CHARACTERIZATION ==========
        % ------------------- 
                fprintf('Finding peak frequency (FFT)...\n')
                spindle_det(ch).freqFFT = zeros( length(start_idxs), 1 ); 
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
                    spindle_det(ch).freqFFT(ix) = spindleFreqs(lowFrInd+spiMaxPowInd-1);                    

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
                spindle_det(ch).num_peaks(ctr) = length( intersect(start_idxs(ix):end_idxs(ix),locs) ) / 2;

                % Peaks and Troughs    
                spindle_det(ch).filtPeakIdxAmp(ix)   = max( real(Data_bp(start_idxs(ix):end_idxs(ix),ch)) );
                spindle_det(ch).filtPeakIdx(ix)      = start_idxs(ix)-1 + find( real(Data_bp(start_idxs(ix):end_idxs(ix),ch))==max( real(Data_bp(start_idxs(ix):end_idxs(ix),ch)) ), 1, 'first'  );    
                spindle_det(ch).filtTroughIdxAmp(ix) = min( real(Data_bp(start_idxs(ix):end_idxs(ix),ch)) );
                spindle_det(ch).filtTroughIdx(ix)    = start_idxs(ix)-1 + find( real(Data_bp(start_idxs(ix):end_idxs(ix),ch))==min( real(Data_bp(start_idxs(ix):end_idxs(ix),ch)) ), 1, 'first'  );    

                spindle_det(ch).rawPeakIdxAmp(ix)    = max( real(Data(start_idxs(ix):end_idxs(ix),ch)) );
                spindle_det(ch).rawPeakIdx(ix)       = start_idxs(ix)-1 + find( real(Data(start_idxs(ix):end_idxs(ix),ch))==max( real(Data(start_idxs(ix):end_idxs(ix),ch)) ), 1, 'first' );    
                spindle_det(ch).rawTroughIdxAmp(ix)  = min( real(Data(start_idxs(ix):end_idxs(ix),ch)) );
                spindle_det(ch).rawTroughIdx(ix)     = start_idxs(ix)-1 + find( real(Data(start_idxs(ix):end_idxs(ix),ch))==min( real(Data(start_idxs(ix):end_idxs(ix),ch)) ), 1, 'first' );    

                spindle_det(ch).envPeakIdx(ix)       = start_idxs(ix)-1 + find( real(Data_env(start_idxs(ix):end_idxs(ix),ch))==max( real(Data_env(start_idxs(ix):end_idxs(ix),ch)) ), 1, 'first' );    
                spindle_det(ch).envPeakIdxAmp(ix)    = max( real(Data_env(start_idxs(ix):end_idxs(ix),ch)) );    
                spindle_det(ch).envPeakIdxTraceVal(ix) = real(Data(spindle_det(ch).envPeakIdx(ix),ch));

            end

            if isfield(spindle_det(ch),'num_peaks')
                spindle_det(ch).num_peaks = spindle_det(ch).num_peaks;
            else
                spindle_det(ch).num_peaks = 0;
            end
            spin = start_idxs;

            spindle_det(ch).label         = ch_names{ch};
            spindle_det(ch).sample        = spin;
            spindle_det(ch).ripple_count = length(spin);
            spindle_det(ch).backgr_std    = signalstd(ch);
            spindle_det(ch).datalen       = size(Data,1);
            spindle_det(ch).startSample   = spin;
            spindle_det(ch).endSample     = end_idxs;
            spindle_det(ch).FWHMstartSample   = adj_start_idxs';
            spindle_det(ch).FWHMendSample     = adj_end_idxs';
            spindle_det(ch).duration      = (spindle_det(ch).endSample-spindle_det(ch).startSample).*(1/Fs)';
            spindle_det(ch).number        = length(spin);
            spindle_det(ch).num_peaks     = spindle_det(ch).num_peaks';
            spindle_det(ch).FWHMduration  = spindle_det(ch).FWHMduration';
        
        %-----------
        %===============================

            % Peaks and Troughs    
            spindle_det(ch).filtPeakIdxAmp   = spindle_det(ch).filtPeakIdxAmp';
            spindle_det(ch).filtPeakIdx      = spindle_det(ch).filtPeakIdx';    
            spindle_det(ch).filtTroughIdxAmp = spindle_det(ch).filtTroughIdxAmp';
            spindle_det(ch).filtTroughIdx    = spindle_det(ch).filtTroughIdx';    

            spindle_det(ch).rawPeakIdxAmp    = spindle_det(ch).rawPeakIdxAmp';
            spindle_det(ch).rawPeakIdx       = spindle_det(ch).rawPeakIdx';    
            spindle_det(ch).rawTroughIdxAmp  = spindle_det(ch).rawTroughIdxAmp';
            spindle_det(ch).rawTroughIdx     = spindle_det(ch).rawTroughIdx';    

            spindle_det(ch).envPeakIdx       = spindle_det(ch).envPeakIdx';    
            spindle_det(ch).envPeakIdxAmp    = spindle_det(ch).envPeakIdxAmp';    
            spindle_det(ch).envPeakIdxTraceVal = spindle_det(ch).envPeakIdxTraceVal';
            spindle_det(ch).freqFFT          = spindle_det(ch).freqFFT;

        %-----------    
        %===============================
               
            
            
    end  % End the loop for each channel
    
    fprintf('\ndone.\n')

end % end function

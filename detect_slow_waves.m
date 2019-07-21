function [sw] = detect_slow_waves(eeg_slow, label, Fs)

% fprintf('============================================================\n')
% fprintf('============================================================\n')
% fprintf('============================================================\n')
% 
% fprintf('This version was last edited Feb. 5th, 2019\n')    
% 
% fprintf('============================================================\n')
% fprintf('============================================================\n')
% fprintf('============================================================\n')

    % Slow wave detection algorithm from "Slow oscillations orchestrating 
    % fast oscillations and memory consolidation
    % by Matthias Molle and and Jan Born
    % implemented by CD on Feb 2, 2015.
    
    % Input: eeg_slow -> one channel of filtered EEG signal in the SW band, 
    % and sampling frequency Fs 
    % Returns sw a structure with detected slow waves and their parameters
    % and sig, a signal containing the SW-filtered EEG signal and zero
    % values otherwise. 

    % 1. As a first step, filter the EEG signal in the slow wave band
    % eeg_slow = eegfilt_mine(eeg_sig,Fs,band(1),band(2)); 

    % 2. Compute all points of positive to negative zero crossings 
    s         = diff(sign(eeg_slow(1,:)));
    indx_down = find(s<0);     % positive to negative
    indx_up   = find(s>0) + 1; % negative to positive (just to define positive slope)
    
    % Make sure that the first zero crossing is from + to -
    if indx_up(1) < indx_down(1)
        indx_up(1) = [];
    end
    
    % 3. Measure the length of all intervals of postive to negative zero
    % crossings (t) is measured
    t = diff(indx_down)/Fs; % sec
    
    % 4.For intervals with a length of 0.8<=t<=2s, the averages of
    % their negative peak amplitudes (x), the positive peak amplitudes (y)
    % and their difference (y-x) are calculated.
    ev = find(t>=0.8 & t<=2); %This means we are looking at SWs between 0.5 and 1.25Hz

    for i = 1:length(ev)

        [x(i),ix(i)]  = min(eeg_slow(1,indx_down(ev(i)):indx_down(ev(i)+1)));  %#ok; negative peak amplitudes 
        [y(i),iy(i)]  = max(eeg_slow(1,indx_down(ev(i)):indx_down(ev(i)+1)));  %#ok; positive peak amplitudes
        pp(i)         = abs(y(i))+abs(x(i));                                   %#ok; peak to peak amplitude   
        dur(i)        = t(ev(i));                                              %#ok; duration of SW in seconds
        indx(i)       = ix(i) + indx_down(ev(i)) -1;                           %#ok; index [in samples] of SW min 
        indy(i)       = iy(i) + indx_down(ev(i)) -1;                           %#ok; index [in samples] of SW max
        neg_slope(i)  = -x(i)/(((ix(i) - 1)/Fs)*1000);                         %#ok; Negative slope [uV/msec]
        pos_slope(i)  = -x(i)/(((indx_up(ev(i)) - indx_down(ev(i)) - ix(i) + 1)/Fs)*1000); %#ok; Positive slope (two zero crossings) 
        pp_slope(i)   = pp(i)/(((indy(i) - indx(i))/Fs)*1000);                 %#ok; Peak-to-Peak slope
        
    end

    % 5. Intervals of positive to negative zero crossings in which 
    %   the corresponding negative peak amplitude x is lower than a threshold
    %   (thr) multiplied by the average of all x, and the corresponding 
    %   amplitude difference (positive peak minus negative peak) is
    %   larger than thr times the average of all (y-x), are marked as 
    %   slow oscillation epochs, with thr set to >0.75.
    % 
    % Decided to remove peak-to-peak amplitude threshold as well as the
    % amplitude threshold for the negative deflection. An amplitude
    % threshold can be defined more flexibly after the initial detection of
    % all the candidate SOs
    
    % store SW parameters in a structure
    sw.startSample = indx_down(ev)';                   % Start of SW (sample)
    sw.endSample   = indx_down(ev+1)';                 % End of SW (sample)
    sw.indx        = indx';                            % index at min (sample)
    sw.indy        = indy';                            % index at max (sample)
    sw.neg_peak    = x';                               % min negative peak in uV
    sw.pos_peak    = y';                               % max positive peak in uV
    sw.dur         = dur';                             % SW duration in sec
    sw.pp          = pp';                              % SW peak to peak amplitude 
    sw.num         = length(sw.startSample)';          % count SWs
    sw.den         = sw.num/(length(eeg_slow)/Fs/60)'; % SW density
    sw.neg_slope   = neg_slope';                       % Negative slope 
    sw.pos_slope   = pos_slope';                       % Positive slope  
    sw.pp_slope    = pp_slope';                        % Peak-to-Peak slope
    sw.label       = label';                           % EEG electrode
    sw.datalen     = length(eeg_slow)';                % length of data in samples

% flip for consistency with other detectors
% sw.startSample = sw.startSample(:);                  % Start of SW (sample)
% sw.endSample   = sw.endSample(:);                    % End of SW (sample)
% sw.indx        = sw.indx(:);                         % index at min (sample)
% sw.indy        = sw.indy(:);                         % index at max (sample)
% sw.neg_peak    = sw.neg_peak(:);                     % min negative peak in uV
% sw.pos_peak    = sw.pos_peak(:);                     % max positive peak in uV
% sw.dur         = sw.dur(:);                          % SW duration in sec
% sw.pp          = sw.pp(:);                           % SW peak to peak amplitude 
% sw.neg_slope   = sw.neg_slope(:);                    % Negative slope 
% sw.pos_slope   = sw.pos_slope(:);                    % Positive slope  
% sw.pp_slope    = sw.pp_slope(:);                     % Peak-to-Peak slope
% sw.label       = sw.label(:);                        % EEG electrode
    
% Note from paper: Importantly, based essentially on zero crossings of the filtered signal, 
% this algorithm takes into account the fact that the depolarizing positive 
% half-wave of the slow oscillations, although lower in amplitude, 
% typically exhibits a distinctly longer duration than the relatively 
% short-lived and sharper hyperpolarizing negative half-wave. *take into account? not really*

end



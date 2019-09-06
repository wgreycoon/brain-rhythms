function [signal] = iirnotchfiltering( fs, signal, linefreq )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%Set IIR Filter Parameters
notchfilter.fcenter = [linefreq,linefreq*2,linefreq*3,linefreq*4];
notchfilter.fcenter = notchfilter.fcenter( find( notchfilter.fcenter<(fs/2) ) ); %#ok
notchfilter.bw      = ones(1,length(notchfilter.fcenter)).*0.001;
notchfilter.ID      = 'notch'    ;    % Identifier 
notchfilter.perform = true       ;    % toggle to apply or not apply filter  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LINE-NOISE FILTERING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the IIR-peak filter coefficients in a,b format 
for idx = 1:length(notchfilter.fcenter),
    notch{idx}.wo = notchfilter.fcenter(idx)/(fs/2);  
    notch{idx}.bw = notchfilter.bw(idx);
    [notch{idx}.b,notch{idx}.a] = iirnotch(notch{idx}.wo,notch{idx}.bw);  
end

% do notch filtering...
fprintf(1, '> Notch filtering signal at %d Hz...\n',linefreq);
fprintf(1,'[');
% for each channel
parfor idx_channel=1:size(signal,2),
    % get the signal for this channel
    signal_preliminary = double(signal(:,idx_channel));
    % remove all harmonics of line-noise
    for idx = 1:length(notchfilter.fcenter), %#ok<PFBNS>
        signal_preliminary = filtfilt(notch{idx}.b,notch{idx}.a,signal_preliminary); %#ok<PFBNS>
    end 
    % return the signal
    signal(:,idx_channel) = single(signal_preliminary);
    fprintf(1,'.');
end
fprintf(1,'] done\n');

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
end
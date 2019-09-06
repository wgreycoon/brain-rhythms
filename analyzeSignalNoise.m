function [ settings ] = analyzeSignalNoise( fs, signal, chs2exclude, linefreq, parswitch )
% 
% -------------------------------------------------------------------------
% analyzeSignalNoise is a function designed to look at ECoG data shared with
% MGH in a collaboration with MISTI in Colombia and MIT in Boston. It
% assesses noise in an array of channels x time by estimating line noise in
% each channel using an IIR peak filter, removing common noise with a
% common average reference, and reporting the reduction in line noise that
% results from this spatial filtering. It also identifies and reports
% channels whose noise levels are > 1.5 standard deviations from the mean.
% ------------------------------------------------------------------------- 
%
% Dependencies
% ---------------------------
% functions:
%     N/A
%
% ---------------------------
%
% Input 
% ---------------------------
%     fs           :    sampling frequency (in Hz)
%     signal       :    array of signal data (time x channels)
%     chs2exclude  :    optional list of channels to exclude from analyses,
%                       including mean/std calculation and spatial filter.
%                       Leave as an empty set ('[]') to ignore.
%     linefreq     :    frequency of power grid AC (ex. 50Hz for Europe,
%                       60Hz for USA)
%     parswitch    :    logical '1' or '0' to control parallel for loop use 
%
%     Example = analyzeNoiseMGH( 100, ecog_data, [], 60 );
%
% ---------------------------
%
% Output
% ---------------------------
%     settings     :    struct with new fields:
%          - channels
%          - channels_selected
%          - channels_noise
%          - noisefiltering.notch.fcenter
%          - noisefiltering.notch.bw
%          - NoiseResults.LineNoiseBefore
%          - NoiseResults.LineNoiseAfter
%
% ---------------------------
%
% Author(s):   William Coon   wcoon@mgh.harvard.edu
%
% Last edited: July 28, 2018 1904hrs
%
% -------------------------------------------------------------------------
%%

% set any a priori channel exclusion
settings.channels2exclude = chs2exclude;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE LINE-NOISE FILTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define the line noise frequency and bandwidth
peak.fcenter = linefreq;
peak.bw      = 0.001;

% calculate the IIR-peak filter coefficients in a,b format 
peak.wo = peak.fcenter/(fs/2);  
peak.bw = peak.bw;
[peak.b,peak.a] = iirpeak(peak.wo,peak.bw);  


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEASSURE LINE-NOISE POWER BEFORE SIGNAL PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, '> Measuring %d Hz noise power before signal processing \n',linefreq);

fprintf(1,'[');
if parswitch
parfor idx_channel=1:size(signal,2),
    % calculate average root-mean-square of the line-noise
    signal_noise_before(idx_channel) = mean(sqrt(filtfilt(peak.b,peak.a,double(signal(:,idx_channel))).^2),1); %#ok<PFBNS>
    fprintf(1,'.');    
end
else
for idx_channel=1:size(signal,2),
    % calculate average root-mean-square of the line-noise
    signal_noise_before(idx_channel) = mean(sqrt(filtfilt(peak.b,peak.a,double(signal(:,idx_channel))).^2),1); %#ok<PFBNS>
    fprintf(1,'.');    
end
end
fprintf(1,'] done\n');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND DATA CHANNELS WITH SIGNIFICANT LINE-NOISE POWER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, '> Searching for channels with significant %d Hz noise \n',linefreq);

if isfield(settings,'channels_noise'),    settings = rmfield(settings,'channels_noise');      end
if isfield(settings,'channels_selected'), settings = rmfield(settings,'channels_selected');   end
if isfield(settings,'channels'),          settings = rmfield(settings,'channels');            end

settings.channels = sort(setdiff(1:size(signal,2),settings.channels2exclude));

fprintf(1,'[');

    % check if any channels are left and 
    if ~isempty(settings.channels),
        % calculate the common average reference signal 
        signal_mean = mean(signal(:,settings.channels),2);
    else
        % if no channels is left then the common average reference signal is zero
        signal_mean = zeros(size(signal,1),1);
    end

    % for each channel on this amp
    for idx_ch=settings.channels,

        % subtract the common average signal from each channel of this amp          
        signal_preliminary = double(signal(:,idx_ch)) - double(signal_mean);

        % calculate the residual line-noise root-mean-square power
        signal_noise(idx_ch) = mean(sqrt(filtfilt(peak.b,peak.a,signal_preliminary).^2),1); 
        fprintf(1,'.');

    end

fprintf(1,'] done\n');

% find those channels for which the line-noise power is 1.5 standard deviations higher than the average 
settings.channels_noise    = find(signal_noise > (mean(signal_noise(settings.channels))+1.5*std(signal_noise(settings.channels))));
settings.channels_noise    = intersect(settings.channels_noise,settings.channels);
settings.channels_selected = setdiff(settings.channels,settings.channels_noise);

fprintf(1, '> Found %d channels with significant %d Hz noise: ',length(settings.channels_noise),linefreq);
fprintf(1, '%d ',settings.channels_noise);
fprintf(1, '\n');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REMOVE COMMON NOISE USING A COMMON-AVERAGE REFERENCE FILTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf(1, '> Common average filtering signal \n');

% now calculate the common-average reference signal again, without those channels that have significant line-noise
fprintf(1,'[');

    % exclude the channels that had signifiant line-noise
    list_channels = settings.channels_selected;

    % check if any channels are left and 
    if ~isempty(list_channels),   

        % calculate the common average reference signal 
        signal_mean = mean(signal(:,intersect(settings.channels,settings.channels_selected)),2);

        % subtract the common average signal from each channel of this amp
        for idx_ch=settings.channels,
            signal(:,idx_ch) = signal(:,idx_ch) - signal_mean;
            fprintf(1,'.');
        end
    else
        % if no channel is left then the signal is not filtered
        for idx_ch=settings.channels,
            signal(:,idx_ch) = signal(:,idx_ch);
            fprintf(1,'.');
        end
    end

fprintf(1,'] done\n');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEASURE LINE-NOISE POWER AFTER SIGNAL PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, '> Measuring %d Hz noise power after signal processing \n',linefreq);

fprintf(1,'[');
% for each channels calculate the root-mean-square line-noise power
if parswitch
parfor idx_channel=1:size(signal,2),
    signal_noise_after(idx_channel) = mean(sqrt(filtfilt(peak.b,peak.a,double(signal(:,idx_channel))).^2),1); %#ok<PFBNS>
    fprintf(1,'.');    
end
else
for idx_channel=1:size(signal,2),
    signal_noise_after(idx_channel) = mean(sqrt(filtfilt(peak.b,peak.a,double(signal(:,idx_channel))).^2),1); %#ok<PFBNS>
    fprintf(1,'.');    
end
end
fprintf(1,'] done\n');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REPORT LINE-NOISE POWER AFTER SIGNAL PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1, '> Reduced %d Hz noise from %.2f to %.2f uV',linefreq,mean(signal_noise_before(settings.channels_selected)),mean(signal_noise_after(settings.channels_selected)));
fprintf(1, '\n');

settings.NoiseResults.LineNoiseBefore = sprintf('%.2f uV',mean(signal_noise_before(settings.channels_selected)));
settings.NoiseResults.LineNoiseAfter  = sprintf('%.2f uV',mean(signal_noise_after(settings.channels_selected)));
    
end
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
    n     = cfg.n;                                      % number of cycles in wavelet
    time  = -cfg.halfwidth:1/srate:cfg.halfwidth;       % time, from -t to t second in steps of 1/sampling-rate
    s     = n/(2*pi*f);                                 % s = standard deviation of Gaussian, derive from 'n' number of cycles

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
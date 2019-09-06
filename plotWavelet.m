function [] = plotWavelet( wavelet_in, cfg )

    figure
    subplot(221)
    % show the projection onto the real axis
    plot3(cfg.time,real(wavelet_in),imag(wavelet_in),'m')
    xlabel('Time (s)'), ylabel('real axis')
    view(0,90)
    title('Projection onto real and time axes')

    % show the projection onto the imaginary axis
    subplot(222)
    plot3(cfg.time,real(wavelet_in),imag(wavelet_in),'g')
    xlabel('Time (s)'), ylabel('imaginary axis')
    view(0,0)
    title('Projection onto imaginary and time axes') 

    % plot projection onto real and imaginary axes
    subplot(223)
    plot3(cfg.time,real(wavelet_in),imag(wavelet_in),'k')
    ylabel('real axis'), zlabel('imag axis')
    view(90,0)
    title('Projection onto imaginary and time axes')

    % plot real and imaginary projections simultaneously
    subplot(224)
    plot(cfg.time,real(wavelet_in),'b')
    hold on
    plot(cfg.time,imag(wavelet_in),'b:')
    legend({'real part';'imaginary part'})


end
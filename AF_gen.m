function [AF_complex, AF_abs, AF_log] = AF_gen(signal, Tp, fs, Ba, Nf, up_sample, lpf, flag)
    lpf_n = (length(lpf)-1)/2;
    fsr = fs * up_sample;
    Nr = ceil(Tp*fsr);
    dt = 1/fsr;
    tr = (0:Nr-1)*dt;
    t = [(Nr-1:-1:1)*(-dt),0,(1:Nr-1)*dt];
    df = Ba/Nf;
    f = [(Nf-1:-1:1)*(-df),0,(1:Nf-1)*df]; 
    % 插值
    signal_interprt = zeros(1, size(signal, 2)*up_sample);
    signal_interprt(1:up_sample:end) = signal;
    signal_interprt_tmp = conv(signal_interprt, lpf);
    signal_interprt_fnl = signal_interprt_tmp(lpf_n+1:end-lpf_n);
    if flag
        signal_interprt_fnl = signal;
    end
    % figure;
    % fre = 0 : fsr/length(signal_interprt_fnl) : fsr - fsr/length(signal_interprt_fnl);
    % plot(fre, 20*log10(abs(fft(signal_interprt_fnl))));
    % 计算模糊函数
    AF = zeros(2*Nr-1, 2*Nf-1);
    for n = 1:2*Nf-1
        AF(:,n) = xcorr(signal_interprt_fnl.*exp(1i*2*pi*f(n).*tr), signal_interprt_fnl).';
    end
    AF_complex = AF;
    AF_abs = abs(AF) ./ max(max(abs(AF)));
    AF_log = 20*log10(AF_abs);
    % 
    % figure;
    % plot(f, AF_log(Nr,:))
    % hold on;
    % ylim([-50 0]);
    % 0多普勒-时延剖面
    % figure;
    plot(t, AF_log(:,Nf));
    hold on;
end
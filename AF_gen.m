function [AF_complex, AF_abs, AF_log] = AF_gen(signal, Tp, fsr, Ba, Nf, up_sample)
    fsr = fsr*up_sample;
    Nr = ceil(Tp*fsr);
    dt = 1/fsr;
    tr = (0:Nr-1)*dt;
    df = Ba/Nf;
    t = [(Nr-1:-1:1)*(-dt),0,(1:Nr-1)*dt];
    f = [(Nf-1:-1:1)*(-df),0,(1:Nf-1)*df]; 
    % 计算模糊函数
    AF = zeros(2*Nr-1, 2*Nf-1);
    for n = 1:2*Nf-1
        AF(:,n) = xcorr(signal.*exp(1i*2*pi*f(n).*tr), signal).';
    end
    AF_complex = AF;
    AF_abs = abs(AF) ./ max(max(abs(AF)));
    AF_log = 20*log10(AF_abs);
    % 0时延-多普勒剖面
    % figure;
    plot(f/100e6, AF_log(Nr,:))
    hold on;
    ylim([-50 0]);
    % 0多普勒-时延剖面
    % figure;
    % Nr = ceil(Tp*fsr*up_sample);
    % dt = 1/fsr/up_sample;
    % t = [(Nr-1:-1:1)*(-dt),0,(1:Nr-1)*dt];
    % AF_time_seg = signal_interpret(AF_complex(:,Nf), up_sample);
    % plot(20*log10(abs(AF_time_seg)./max(abs(AF_time_seg))));

    % plot(t/Tp, AF_log(:,Nf));
    % hold on;
end
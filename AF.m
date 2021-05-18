%% 模糊函数仿真
% OFDM参数
% OFDM块时长5us，带宽100MHz，10个OFDM块，采样率128MHz
close all; clear all;
B = 100e6;
T_OFDM = 5e-6;
Num_OFDM = 10;
% radar参数(基带，未加载频)
C = 3e8;
Br = B; 
Tr = T_OFDM*Num_OFDM;
fsr = 128e6; % 根据OFDM时长和点数计算 5us，512点IFFT
fd = 50e3; % 多普勒是50kHz
Ba = fd; % 多普勒带宽是50kHz
%% 时间域
Nr = ceil(Tr*fsr);
dt = 1/fsr;
% t = (0:2*Nr-1)*dt + (-Tr)+dt; % ??? 问题1：这个t取不到0
t = [(Nr-1:-1:1)*(-dt),0,(1:Nr-1)*dt];
tr = (0:Nr-1)*dt;
Nf = 1000;
df = Ba/Nf;
f = [(Nf-1:-1:1)*(-df),0,(1:Nf-1)*df];;
%% 生成发射信号
st = Gen_OFDM(Num_OFDM,B,T_OFDM);
comp = fftshift(ifft(fft(st,2*Nr-1).*conj(fft(st,2*Nr-1))));
comp_log = 20*log10(abs(comp)./max(abs(comp)));
% figure,plot(t,comp_log);title('点扩展函数'),xlabel('时延'),ylabel('幅度/dB'); 
% figure,plot(abs(fft(st))./max(abs(fft(st))));title('发射信号频域'); xlabel('采样点');ylabel('幅度')% ???点扩展函数怎么插值
%% 生成模糊函数
AF_OFDM = zeros(2*Nr-1,2*Nf-1);
for n = 1: 2*Nf-1
    st_time = st.*exp(1i*2*pi*f(n).*tr); % ?????就这块有问题,很好，没问题了
    st_time_fft = fft(st_time,2*Nr-1);
    st_freq = fftshift(ifft(st_time_fft.*conj(fft(st,2*Nr-1))));
    AF_OFDM(:,n) = st_freq.';
end

% figure; plot(20*log10(abs(fftshift(fft(st .* conj(st))))))

%% 距离速度二维图
 AF_OFDM = abs(AF_OFDM)./max(max(abs(AF_OFDM)));
 AF_OFDM_log = 20*log10(AF_OFDM);
 figure,mesh(f,t,AF_OFDM); title('ofdm的模糊函数');xlabel('多普勒频移');ylabel('时延');zlabel('幅度(dB)');%xlim([-1000,1000]),zlim([-100,0])
%% 0时延—多普勒剖面
figure,imagesc(AF_OFDM);
figure,plot(f/1e3,AF_OFDM_log(Nr,:));title('0时延-多普勒剖面');xlabel('多普勒/KHz');ylabel('幅度（归一化）');
%% 0多普勒-时延剖面
figure,plot(t/1e3,AF_OFDM_log(:,Nf));title('0多普勒-时延剖面');xlabel('时延/us');ylabel('幅度（归一化）');
%% 分辨率计算
% 速度分辨率   单个脉冲的多普勒分辨率是时宽的倒数 【参考文献：雷达通信一体化波形设计与处理方法研究2.3.1节】
fprintf('多普勒分辨率：')
fprintf('Dat_f = %d KHz\n',1/Tr/1e3);
% 距离分辨率  距离分辨率等于c/2B 【参考文献：雷达通信一体化波形设计与处理方法研究2.3.1节】
fprintf('距离分辨率：');
fprintf('Dat_R = %1.2f 米\n',C/2/B);
%% 多普勒容限计算
% 多普勒容限是速度模糊函数珠峰的3dB带宽  【参考文献：雷达通信一体化波形设计与处理方法研究2.3.1节】
ThreedB_duppoler = ThreedB(AF_OFDM(Nr,:).',4);
fprintf('多普勒容限：')
fprintf('Doppuler_tolera = %1.3f KHz\n',ThreedB_duppoler*df/1e3);
%% 峰值旁瓣比计算
[IRW,PSLR,ISLR] = IRW_PSLR_ISLR_zxx(AF_OFDM(:,Nf),10,0);
fprintf('峰值旁瓣比：');
fprintf('PSLR = %1.3f dB\n',PSLR);
IRW
%% 画波形脉压和线性调频信号的对比(相同带宽下)
load('s_abs_log.mat')
tau = 0.5e-3;
td = tau/Nr:tau/Nr:tau;
k = B/tau;
LFM = exp(1i*pi*k*td.^2);
LFM_compres = ifft(fft(LFM,2*Nr-1).*conj(fft(LFM,2*Nr-1)));
LFM_compres_interp = echopc_inter(LFM_compres.',10);
LFM_compres_interp_log = 20*log10(abs(LFM_compres_interp)./max(abs(LFM_compres_interp)));
% 截取线性调频和波形脉压的对比
winNum = 1200;
len = length(LFM_compres_interp_log);
x_ra = len/2-winNum/2+1:len/2+winNum/2;
figure,plot(LFM_compres_interp_log(x_ra+5));hold on,plot(s_abs_log(x_ra));
xlabel('采样点'),ylabel('幅度（dB）');ylim([-80,0])
legend('线性调频信号','OFDM信号');

A = abs(sin(pi*t*v)./(pi*t*v));
A(isnan(A)) = 1;




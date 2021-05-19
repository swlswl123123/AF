close all; clear all;
%% 模糊函数仿真
% 信号时长为50us, 带宽为100MHz
B = 100e6;
T_OFDM = 5e-6;
Num_OFDM = 10;
% 插值参数
up_sample = 8;
% radar参数(基带，未加载频)
C = 3e8;
Tr = T_OFDM*Num_OFDM; % 脉冲时长
fsr = 128e6; % 采用率
Ba = 50e3; % 多普勒带宽是50kHz
Nf = 1000; % 频偏采样点数
Nr = ceil(Tr*fsr); % 时间采样点数
tr = (0:Nr-1)/fsr; % 时长
%% 生成OFDM发射信号
st = Gen_OFDM(Num_OFDM,B,T_OFDM);
% 插值
st = signal_interpret(st, 8);
%% 生成模糊函数
[AF_OFDM_cmp, AF_OFDM_abs, AF_OFDM_log] = AF_gen(st, Tr, fsr, Ba, Nf, 8);
%% 生成LFM发射信号
k = B/Tr;
% fsr = fsr * 10;
% Nr = ceil(Tr*fsr);
% tr = (0:Nr-1)/fsr;
LFM = exp(1i*pi*k*tr.^2);
% 插值
LFM = signal_interpret(LFM, 8);
[AF_LFM_cmp, AF_LFM_abs, AF_LFM_log] = AF_gen(LFM, Tr, fsr, Ba, Nf, 8);

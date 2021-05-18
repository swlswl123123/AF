%% ģ����������
% OFDM����
% OFDM��ʱ��5us������100MHz��10��OFDM�飬������128MHz
close all; clear all;
B = 100e6;
T_OFDM = 5e-6;
Num_OFDM = 10;
% radar����(������δ����Ƶ)
C = 3e8;
Br = B; 
Tr = T_OFDM*Num_OFDM;
fsr = 128e6; % ����OFDMʱ���͵������� 5us��512��IFFT
fd = 50e3; % ��������50kHz
Ba = fd; % �����մ�����50kHz
%% ʱ����
Nr = ceil(Tr*fsr);
dt = 1/fsr;
% t = (0:2*Nr-1)*dt + (-Tr)+dt; % ??? ����1�����tȡ����0
t = [(Nr-1:-1:1)*(-dt),0,(1:Nr-1)*dt];
tr = (0:Nr-1)*dt;
Nf = 1000;
df = Ba/Nf;
f = [(Nf-1:-1:1)*(-df),0,(1:Nf-1)*df];;
%% ���ɷ����ź�
st = Gen_OFDM(Num_OFDM,B,T_OFDM);
comp = fftshift(ifft(fft(st,2*Nr-1).*conj(fft(st,2*Nr-1))));
comp_log = 20*log10(abs(comp)./max(abs(comp)));
% figure,plot(t,comp_log);title('����չ����'),xlabel('ʱ��'),ylabel('����/dB'); 
% figure,plot(abs(fft(st))./max(abs(fft(st))));title('�����ź�Ƶ��'); xlabel('������');ylabel('����')% ???����չ������ô��ֵ
%% ����ģ������
AF_OFDM = zeros(2*Nr-1,2*Nf-1);
for n = 1: 2*Nf-1
    st_time = st.*exp(1i*2*pi*f(n).*tr); % ?????�����������,�ܺã�û������
    st_time_fft = fft(st_time,2*Nr-1);
    st_freq = fftshift(ifft(st_time_fft.*conj(fft(st,2*Nr-1))));
    AF_OFDM(:,n) = st_freq.';
end

% figure; plot(20*log10(abs(fftshift(fft(st .* conj(st))))))

%% �����ٶȶ�άͼ
 AF_OFDM = abs(AF_OFDM)./max(max(abs(AF_OFDM)));
 AF_OFDM_log = 20*log10(AF_OFDM);
 figure,mesh(f,t,AF_OFDM); title('ofdm��ģ������');xlabel('������Ƶ��');ylabel('ʱ��');zlabel('����(dB)');%xlim([-1000,1000]),zlim([-100,0])
%% 0ʱ�ӡ�����������
figure,imagesc(AF_OFDM);
figure,plot(f/1e3,AF_OFDM_log(Nr,:));title('0ʱ��-����������');xlabel('������/KHz');ylabel('���ȣ���һ����');
%% 0������-ʱ������
figure,plot(t/1e3,AF_OFDM_log(:,Nf));title('0������-ʱ������');xlabel('ʱ��/us');ylabel('���ȣ���һ����');
%% �ֱ��ʼ���
% �ٶȷֱ���   ��������Ķ����շֱ�����ʱ��ĵ��� ���ο����ף��״�ͨ��һ�廯��������봦�����о�2.3.1�ڡ�
fprintf('�����շֱ��ʣ�')
fprintf('Dat_f = %d KHz\n',1/Tr/1e3);
% ����ֱ���  ����ֱ��ʵ���c/2B ���ο����ף��״�ͨ��һ�廯��������봦�����о�2.3.1�ڡ�
fprintf('����ֱ��ʣ�');
fprintf('Dat_R = %1.2f ��\n',C/2/B);
%% ���������޼���
% �������������ٶ�ģ����������3dB����  ���ο����ף��״�ͨ��һ�廯��������봦�����о�2.3.1�ڡ�
ThreedB_duppoler = ThreedB(AF_OFDM(Nr,:).',4);
fprintf('���������ޣ�')
fprintf('Doppuler_tolera = %1.3f KHz\n',ThreedB_duppoler*df/1e3);
%% ��ֵ�԰�ȼ���
[IRW,PSLR,ISLR] = IRW_PSLR_ISLR_zxx(AF_OFDM(:,Nf),10,0);
fprintf('��ֵ�԰�ȣ�');
fprintf('PSLR = %1.3f dB\n',PSLR);
IRW
%% ��������ѹ�����Ե�Ƶ�źŵĶԱ�(��ͬ������)
load('s_abs_log.mat')
tau = 0.5e-3;
td = tau/Nr:tau/Nr:tau;
k = B/tau;
LFM = exp(1i*pi*k*td.^2);
LFM_compres = ifft(fft(LFM,2*Nr-1).*conj(fft(LFM,2*Nr-1)));
LFM_compres_interp = echopc_inter(LFM_compres.',10);
LFM_compres_interp_log = 20*log10(abs(LFM_compres_interp)./max(abs(LFM_compres_interp)));
% ��ȡ���Ե�Ƶ�Ͳ�����ѹ�ĶԱ�
winNum = 1200;
len = length(LFM_compres_interp_log);
x_ra = len/2-winNum/2+1:len/2+winNum/2;
figure,plot(LFM_compres_interp_log(x_ra+5));hold on,plot(s_abs_log(x_ra));
xlabel('������'),ylabel('���ȣ�dB��');ylim([-80,0])
legend('���Ե�Ƶ�ź�','OFDM�ź�');

A = abs(sin(pi*t*v)./(pi*t*v));
A(isnan(A)) = 1;




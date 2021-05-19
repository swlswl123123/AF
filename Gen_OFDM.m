%% 
% ����Ҫ�� 100MHz��50us
% ��һ�Σ� 
% ����10�����ţ�ÿ������ʱ��5us����Чʱ��4us��ѭ��ǰ׺����1us����400����Ч���ز���512��IFFT��16QAM����Ϣ��������400Mbit/s
function frame = Gen_OFDM(Num_symbol_data,Br,T_subframe)
% Br = 100MHz
% T_subframe = 5us
% Num_sysmbol_data = 10

% Br = 1e8;
% T_subframe = 5e-6;
% Num_symbol_data = 10;
B = Br;
Ts = T_subframe*4/5; %��Чʱ��
Tg = T_subframe/5;%ѭ��ǰ׺
T = Ts+Tg;
carrier_f = 1/Ts;%���ز����
Num_carrier = floor(B*Ts); %��Ч�ز���
ifft_length = 2^nextpow2(Num_carrier);%IFFT����
cp_length = floor(ifft_length/4); 
Nr = ifft_length + cp_length;
dt = Ts/ifft_length;
Fsr = 1/dt;
modulate_bit = 4;    %1ΪBPSK��2ΪQPSK,TD-LTE���õ���16QAM����Ϊ4
%% �������
fprintf('����');
fprintf('Br = %d MHz\n',B/1e6);
fprintf('�����ʣ�');
fprintf('fsr = %d MHz\n',int32(Fsr/1e6));
fprintf('OFDM��ʱ����');
fprintf('Tg = %d us\n',int32(T/1e-6));
fprintf('�ز�������');
fprintf('Nc = %d\n',Num_carrier);
%% ����bit��
inforSource = randi([0 1],1,Num_carrier*modulate_bit*Num_symbol_data);
% inforSource = ones(1,Num_carrier*modulate_bit*Num_symbol_data);
% save inforSource inforSource;
data_temp1 = reshape(inforSource,modulate_bit,[]).';
lteSymMap = [11 10 14 15 9 8 12 13 1 0 4 5 3 2 6 7];
%% ����ͼӳ��
modulate_data = qammod(bi2de(data_temp1), 2^(modulate_bit), lteSymMap, 'UnitAveragePower', true);
% scatterplot(modulate_data);title('����ͼ');
%% ���ز�ӳ��
carrier_data = reshape(modulate_data,Num_carrier,Num_symbol_data); 
carrier_data = [carrier_data(1:Num_carrier/2,:);zeros(1,Num_symbol_data);carrier_data(Num_carrier/2+1:end,:)]; %����ֱ�����ز�
carrier_data_temp1 = zeros(ifft_length,Num_symbol_data);
carrier_data_temp1(ifft_length/2-Num_carrier/2+1:ifft_length/2+Num_carrier/2+1,:) = carrier_data;%��ӱ�������
carrier_data =fftshift(carrier_data_temp1);
%% OFDM����
time_data_ifft = ifft(carrier_data);
%% ��CP
time_data_cp = [time_data_ifft(end-cp_length+1:end,:);time_data_ifft];
%% ��������
frame = reshape(time_data_cp,1,[]);
%% ��ͼ
% figure,plot(dt*(0:length(frame)-1),abs(frame));
% xlabel('ʱ��/s');ylabel('����');title('���䲨��');

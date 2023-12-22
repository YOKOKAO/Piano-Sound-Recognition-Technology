% 选择要读入的音频文件，并根据不同文件设置不同的BPM及窗口布局
n=1;
switch n
    case 1
        filename="音阶\M1_i1.wav";
        BPM=100;i_range=3;j_range=5;
    case 2
        filename="音阶\H1_i1.wav";
        BPM=100;i_range=3;j_range=5;
    case 3
        filename="音阶\L1_i1.wav";
        BPM=100;i_range=3;j_range=5;    
    case 4
        filename="小星星\star_69.m4a";
        BPM=69;i_range=5;j_range=9;
end

% 读入音频文件，取单声道，归一化
[music,fs]=audioread(filename);
music=music(:,1)./max(music(:,1));
t=(0:length(music)-1)/fs;

% 使用帧峰检测法提取音频包络
Frame_Num=100;
envelope1=[];
for i=1:floor(length(music)/Frame_Num)-1
    temp=music((i-1)*Frame_Num+1:i*Frame_Num);
    envelope1=[envelope1 max(temp)-min(temp)];
end
t_=t(1:Frame_Num:length(t)-1);
t_=t_(1:length(envelope1));
figure(1);subplot(3,1,1);
plot(t,music);hold on;
plot(t_,envelope1);hold off;
title('帧峰包络');ylim([-1.2,1.2]);

% 包络线滤波
envelope2=smooth(envelope1,11);
envelope2=medfilt1(envelope2,11);
envelope2=smooth(envelope2,21);
subplot(3,1,2);
plot(t,music);hold on;
plot(t_,envelope2);hold off;
title('滤波后包络');ylim([-1.2,1.2]);

% 包络线归一化
envelope3=envelope2./max(envelope2);
subplot(3,1,3);
plot(t,music);hold on;
plot(t_,envelope3);hold off;
title('包络归一化');ylim([-1.2,1.2]);

% 寻找包络峰值点
[Amp1,Loc1]=findpeaks(envelope3);
figure(2);subplot(3,1,1);
plot(t,music);hold on;
plot(t_,envelope3,t_(Loc1),Amp1,'o');hold off;
title('寻找包络峰值点');ylim([-1.2,1.2]);

% 设定幅度限制
[Amp2,Loc2]=findpeaks(envelope3,'MinPeakHeight',0.35);
subplot(3,1,2);
plot(t,music);hold on;
plot(t_,envelope3,t_(Loc2),Amp2,'o');hold off;
title('设定峰值点幅度限制');ylim([-1.2,1.2]);

% 设定最小峰间隔
Min_sub_beat=4;
Min_Gap=floor((60/BPM)*fs/Min_sub_beat/Frame_Num);
[Amp3,Loc3]=findpeaks(envelope3,'MinPeakHeight',0.35,'MinPeakDistance',Min_Gap);
subplot(3,1,3);
plot(t,music);hold on;
plot(t_,envelope3,t_(Loc3),Amp3,'o');hold off;
title('设定峰值点最小间隔');ylim([-1.2,1.2]);

% 预加重处理
music=filter([1,-0.9375],1,music);

% 取音程段
Amp_MS=[];F0=0.1;N=fs/F0;
for i=1:length(Loc3)
    if i<length(Loc3)
        X_temp=music(Loc3(i)*Frame_Num:floor(Loc3(i)+(Loc3(i+1)-Loc3(i))*2/3)*Frame_Num);
    else
        X_temp=music(Loc3(i)*Frame_Num:floor(Loc3(i)+(Loc3(i)-Loc3(i-1))*2/3)*Frame_Num);
    end
    X_temp=X_temp.*hanning(length(X_temp));
    Amp_temp=fft(X_temp,N);
    Amp_MS=[Amp_MS abs(Amp_temp)];
end

% 逐个音程画图
figure(3);k=0:N-1;
for i=1:i_range
    for j=1:j_range
        if((i-1)*j_range+j <= length(Loc3))
            subplot(i_range,j_range,(i-1)*j_range+j);
            plot(k*fs/N,Amp_MS(:,(i-1)*j_range+j)/max(Amp_MS(:,(i-1)*j_range+j)));
            axis([0 2000 0 1]);
        end
    end
end

% 画整个音频的短时傅里叶图
M_T=ceil(length(music)/fs);
test_Num=4000;F0_=1;
frame_N=floor(length(music)/test_Num);
frame_t=[0:frame_N-1]*test_Num/fs;
fft_N=fs/F0_;
frame_f=[0:floor(fft_N-1)]*F0_;
frame_A=zeros(fft_N,frame_N);
for i=1:frame_N-1
    frame_A(:,i)=abs(fft(music((i-1)*test_Num+1:i*test_Num).*hamming(test_Num),fft_N));
end
figure(4);subplot(2,1,1);
imagesc(frame_t,frame_f,frame_A);
title('短时傅里叶谱图(test_Num=4000)','Interpreter','none');
axis xy;axis([0 M_T 0 2000]);

frame_t_ms=[Loc3*Frame_Num/fs];
frame_N=length(Loc3);
fft_N=fs/F0_;
frame_f_ms=[0:floor(fft_N-1)]*F0_;
frame_A_ms=zeros(fft_N,frame_N);
for i=2:frame_N
    test_ms=music(Loc3(i-1)*Frame_Num+1:Loc3(i)*Frame_Num);
    frame_A_ms(:,i-1)=abs(fft(test_ms.*hamming(length(test_ms)),fft_N));
end
test_ms=music(Loc3(i)*Frame_Num+1:Loc3(i)*Frame_Num+length(test_ms));
frame_A_ms(:,i)=abs(fft(test_ms.*hamming(length(test_ms)),fft_N));
figure(4);subplot(2,1,2);
imagesc(frame_t_ms,frame_f_ms,frame_A_ms);
title('短时傅里叶谱图(按照音程绘制)');
axis xy;axis([0 M_T 0 2000]);

% 将音频幅度矩阵Amp_MS映射到按键矩阵KeyEnergy_MS中，并逐个绘制键号对应的图及瀑布图
Piano_F=xlsread('Piano_Key_F.xlsx');
KeyEnergy_MS=[];
for j=1:88
    for i=1:length(Loc3)
        KeyEnergy_MS(j,i)=Amp_MS(round(Piano_F(j)/F0),i);
    end
end
KeyEnergy_MS=KeyEnergy_MS./max(KeyEnergy_MS);
figure(5);key=1:88;
for i=1:i_range
    for j=1:j_range
        if((i-1)*j_range+j <= length(Loc3))
            subplot(i_range,j_range,(i-1)*j_range+j);
            plot(key,KeyEnergy_MS(:,(i-1)*j_range+j));
            axis([0 88 0 1]);
        end
    end
end
figure(6);
imagesc(1:length(Loc3),1:88,KeyEnergy_MS);
axis xy;axis([0 M_T 0 88]);

% 建立识别策略，绘制识别结果及瀑布图
Key_Rec=zeros(88,length(Loc3));
for i=1:length(Loc3)
    for j=1:88
        if KeyEnergy_MS(j,i)>0.6
            Key_Rec(j,i)=1;
            break;
        end
    end
end
figure(7);key=1:88;
for i=1:i_range
    for j=1:j_range
        if((i-1)*j_range+j <= length(Loc3))
            subplot(i_range,j_range,(i-1)*j_range+j);
            stem(key+20,Key_Rec(:,(i-1)*j_range+j));
            axis([0 88 0 1.2]);
        end
    end
end
figure(8);
imagesc(1:length(Loc3),1:88,Key_Rec);
axis xy;axis([0 M_T 0 88]);

% 音频节奏评价
ideal_Beat_T=[ones(1,6) 2 ones(1,6) 2 ones(1,6) 2 ones(1,6) 2 ones(1,7) 2 ones(1,6)];
Real_Beat_T=zeros(1,length(Loc3)-1);
Err_Beat_T=zeros(1,length(Loc3)-1);
for i=1:length(Loc3)-1
    X_temp=music(Loc3(i)*Frame_Num:floor(Loc3(i)+(Loc3(i+1)-Loc3(i)))*Frame_Num);
    Real_Beat_T(i)=length(X_temp)/fs*BPM/60;
    Err_Beat_T(i)=Real_Beat_T(i)-ideal_Beat_T(i);
end
figure(9);subplot(2,1,1);
plot(1:length(Loc3)-1,Real_Beat_T,'bo');hold on;
plot(1:length(Loc3)-1,ideal_Beat_T,'r*');grid on;
axis([1 length(Loc3)-1 0 2.2])
subplot(2,1,2);
plot(1:length(Loc3)-1,Err_Beat_T,'bx-');grid on;
axis([1 length(Loc3)-1 -0.5 0.5])
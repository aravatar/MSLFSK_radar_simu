clear;
close all;
c=3e8;
f0=24e9;
lambda=c/f0;
rmax=120;
vmax=240/3.6;
derta_r=1;
snr=10^(-15/20);
sigma_r=0.5;
sigma_v=1.8/3.6;

%% 参数计算
b=c/(2*derta_r);
tcpi=2*derta_r/vmax;
miu=2*b/tcpi;
sigma_v_r=c*sigma_v/(miu*lambda);
sigma_r=min(sigma_r,sigma_v_r);
n0=tcpi*(miu*rmax/c+vmax/lambda);

n=512;%大于n0的2次幂

fincr=b/(n-1);
dert=0.5;
sigm=2/(snr*n*sinc(dert)*sinc(dert));
alph=0.05;
x = norminv(alph/2,0,sigm);
ddfai=abs(x);
rmin_uamb=2*sigma_r*pi/ddfai;
rmax_uamb=rmax/(1-ddfai/(pi));
m=ceil(log2(rmax_uamb/rmin_uamb)+2);
tstep=0.5*tcpi/(m*n);

r_uamb=zeros(1,m-1);
f_shift=zeros(1,m-1);
for k=1:m-1
    r_uamb(k)=rmax_uamb/(2^(k-1));
    f_shift(k)=c/(2*r_uamb(k))-miu*tstep;
end

%% 仿真
r0=76.9;%
v0=25.5;%目标
t=linspace(0,tcpi,m*n*2*10);
r_t=zeros(1,length(t));
td=zeros(1,length(t));
ft=zeros(1,length(t));
fn=f0;
for i=1:n
    ft_t=(fn)*ones(1,10);
    for j=1:length(f_shift)
        ft_t=[ft_t (fn-f_shift(j))*ones(1,10)];
    end
    ft(1+(i-1)*length(ft_t):length(ft_t)+(i-1)*length(ft_t))=ft_t;
    fn=fn+fincr;
end
for i=1:n
    ft_t=(fn)*ones(1,10);
    for j=1:length(f_shift)
        ft_t=[(fn-f_shift(j))*ones(1,10) ft_t];
    end
    ft(m*n*10+1+(i-1)*length(ft_t):m*n*10+length(ft_t)+(i-1)*length(ft_t))=ft_t;
    fn=fn-fincr;
end
figure;
plot(t,ft);
title("发射时频图");

tx=zeros(1,length(t)); %发射信号
rx=zeros(1,length(t)); %接收信号
mix=zeros(1,length(t)); %差频、差拍、拍频、中频信号
r_t = r0 + v0*t; % 距离更新
td = 2*r_t/c;    % 延迟时间
tx=cos(2*pi*ft.*t);
rx=cos(2*pi*ft.*(t-td));
rx=awgn(rx,5,'measured');%接收信号含噪
mix=tx.*rx;
mix=awgn(mix,5,'measured');%混频信号含噪
figure;
plot(t(1:512),tx(1:512));
title("发射波形（部分）");
figure;
plot(t(1:512),rx(1:512));
title("接收波形（部分）");
figure;
plot(t(1:512),mix(1:512));
title("混频波形（部分）");

figure;
plot(t,ft);
hold on;
plot(t+td,ft);
title("收发时频图");
legend ('发','收');

% mix_0=zeros(1,length(t)/10);
% for i=1:length(t)/10
%     mix_0(i)=mix(i*10);
% end
fs=1/(tstep);
mix_5=zeros(m,n);
for i=1:n
    for j=1:m
        mix_5(j,i)=mix((j-1)*10+1+(i-1)*m*10);
    end
end
mix_5_fft=zeros(m,n);
for i=1:m
    mix_5_fft(i,:)=fft(mix_5(i,:));
end
figure;
plot(abs(mix_5_fft(1,:)));
title("1通道LFSK频谱");
wsum=zeros(1,n);
for i=1:m
    wsum(1,:)=wsum(1,:)+abs(mix_5_fft(i,:)).*abs(mix_5_fft(i,:));
end
figure;
plot(wsum);
title("频谱周期图积累");
[a,p]=max(wsum(1:length(wsum)/2));
% p=83;
fmslfsk=fs*p/(m*n);
derta_fai=zeros(1,m-1);
% for i=1:m-1
%     derta_fai(i)=(atan(real(mix_5_fft(i,p))/imag(mix_5_fft(i,p)))-atan(real(mix_5_fft(i+1,p))/imag(mix_5_fft(i+1,p))));
% end
for i=1:m-1
    derta_fai(i)=(angle(mix_5_fft(i,p))-angle(mix_5_fft(i+1,p)));
end
% derta_fai_real=2*pi*(f_shift*2*r0/c-2*v0*f0*tstep/c);
% for i=1:m-1
%     k=round((derta_fai_real(i)-derta_fai(i))/pi);
%     derta_fai(i)=derta_fai(i)+k*pi;
% end

rout=r_uamb(1)*(derta_fai(1)/(2*pi)+fmslfsk*tstep-floor(derta_fai(1)/(2*pi)+fmslfsk*tstep));
derta_fai=2*pi*(rout./r_uamb-fmslfsk*tstep);

r_amb=zeros(1,m-1);
for k=1:m-1
    r_amb(k)=r_uamb(k)*(derta_fai(k)/(2*pi)+fmslfsk*tstep-floor(derta_fai(k)/(2*pi)+fmslfsk*tstep));
end
%解模糊
q=zeros(1,m-1);
for k=2:m-1
    q(k)=round((rout-r_amb(k))/r_uamb(k));
    rout=r_amb(k)+q(k)*r_uamb(k);
end

vout=(fmslfsk*c-miu*2*rout)/(2*f0);

figure;
scatter(r0,v0,'b');
hold on
scatter(rout,vout,'r');
scatter(0,0,'w');%参考点
scatter(rmax,vmax,'w');%参考点
xlabel("距离(m)");
ylabel("速度(m/s)");
title("测距测速结果");
legend ('真实值','测量值');









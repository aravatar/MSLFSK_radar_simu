clear;
c=3e8;
f0=24e9;
lambda=c/f0;
rmax=120;
vmax=240/3.6;
snr=10^(-15/20);
sigma_r=2;
sigma_v=3/3.6;
derta_r=1;
%%
b=c/(2*derta_r);
tcpi=2*derta_r/vmax;
miu=2*b/tcpi;
sigma_v_r=c*sigma_v/(miu*lambda);
sigma_r=min(sigma_r,sigma_v_r);
n0=tcpi*(miu*rmax/c+2*vmax/lambda);

n=512;
fincr=b/(n-1);

dert=0.5;
sigm=2/(snr*n*sinc(dert)*sinc(dert));
alph=0.05;
x = norminv(alph/2,0,sigm);
ddfai=abs(x);
rmin_uamb=2*sigma_r*pi/ddfai;
rmax_uamb=rmax/(1-ddfai/pi);
m=ceil(log2(rmax_uamb/rmin_uamb)+2);
tstep=0.5*tcpi/(m*n);

r_uamb=zeros(1,m-1);
f_shift=zeros(1,m-1);
for k=1:m-1
    r_uamb(k)=rmax_uamb/(2^(k-1));
    f_shift(k)=c/(2*r_uamb(k))+miu*tstep;
end







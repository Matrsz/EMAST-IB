pkg load control 

% Satelite 

is = 1069;
ns=1;
ds=[is 0 0];

% Modo elastico de antena 

i0=1025;
f1=2.0; 
psi1 =0.001;     

w1 =2*pi*f1;
psi1m=psi1*sqrt(is/i0); 
w1m=w1*sqrt(is/i0);

n1=is/i0*[1 2*psi1*w1   w1*w1];
d1=[1 2*psi1m*w1m w1m*w1m];
 

% Modo elastico de los paneles solares
 
i0=300;
f2=0.3; 
psi2 =0.001;     

w2 =2*pi*f2;
psi2m=psi2*sqrt(is/i0); 
w2m=w2*sqrt(is/i0);

n2=is/i0*[1 2*psi2*w2   w2*w2];
d2=[1 2*psi2m*w2m w2m*w2m];

nt=conv(n1,n2); 
dt=conv(d1,d2);

nt=conv(nt,ns); 
dt=conv(dt,ds);

sys=tf(nt,dt);
sysCL=sys/(1+sys);

figure 1; margin(sys);
figure 2; nichols(sys,{0.001,10});
figure 3; step(sysCL,1000,1);

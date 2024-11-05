pkg load control 
pkg load signal

FLEX=1;
FILTRO=0;
DELAY=1;

% Satelite

is = 1069;
ns=1;
ds=[is 0 0];

nt=ns; 
dt=ds;

% La proporcion de Inercia aportada por el apendice flexible se regulara con la constante rmi (i=1,2,...)
% rmi = Inercia del Apendice flexible / Inercia total del satelite = Ia/Is
% Por otro lado para calcular la frecuencia modal ( denominador de la funcion de transferencia ) se necesitan Is y I0=Is-Ia
% I0=Is-Ia=Is-rmi*Is=Is*(1-rmi)


% Modo elastico de antena 

rm1=0.04; 
i0=is*(1-rm1);
f1=2.0; 
psi1 =0.001;     

w1 =2*pi*f1;
psi1m=psi1*sqrt(is/i0); 
w1m=w1*sqrt(is/i0);

n1=is/i0*[1 2*psi1*w1   w1*w1];
d1=[1 2*psi1m*w1m w1m*w1m];
 
% Modo elastico de los paneles solares
 
rm2=0.34; 
i0=is*(1-rm2);
f2=0.3; 
psi2 =0.001;     

w2 =2*pi*f2;
psi2m=psi2*sqrt(is/i0); 
w2m=w2*sqrt(is/i0);

n2=is/i0*[1 2*psi2*w2   w2*w2];
d2=[1 2*psi2m*w2m w2m*w2m];

if (FLEX )
 nt=conv(nt,n1); 
 dt=conv(dt,d1);
 
 nt=conv(nt,n2); 
 dt=conv(dt,d2);
endif

% LEAD-LAG

k=3.17;
psi=1.;

Wn = sqrt ( k/is );
kd = 2*psi*Wn*is;
Td = kd/k;

nc=[k*Td k];
dc=[Td/10 1];

nt=conv(nt,nc); 
dt=conv(dt,dc);

% Delay

Tdelay=1;
ndelay=1;
ddelay=[Tdelay 1];

if(DELAY)
 nt=conv(nt,ndelay); 
 dt=conv(dt,ddelay);
endif

sys=tf(nt,dt);

% Filtro

[N1,D1]=ellip(6,2.,10,1.5,'low','s');
sys_F1=tf(N1,D1);

if(FILTRO)
 sys=sys*sys_F1;
endif

sysCL=sys/(1+sys);

figure 1; margin(sys);
figure 2; nichols(sys,{0.001,10});
figure 3; step(sysCL);

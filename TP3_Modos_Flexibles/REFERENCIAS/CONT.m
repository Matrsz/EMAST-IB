pkg load control 


% Las ganancias k y kd que siguen son para error medido en angulo
% Cuando se usa con error de cuaternion, se deben multiplicar k y kd resultantes de este analisis por 2

% LEAD-LAG

is=1069;
k=1.*10;
psi=1.;

Wn = sqrt ( k/is );
kd = 2*psi*Wn*is;
Td = kd/k;

nc=[k*Td k];
dc=[Td/10 1];

nt=nc;
dt=dc;

% Delay

Tdelay=1;
ndelay=1;
ddelay=[Tdelay 1];

nt=conv(nt,ndelay); 
dt=conv(dt,ddelay);

%
sys=tf(nt,dt);
sysCL=sys/(1+sys);

figure 1; 
margin(sys);


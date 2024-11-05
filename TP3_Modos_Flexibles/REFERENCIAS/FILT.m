pkg load control 
pkg load signal


[N1,D1]=ellip(4,2.,30,1.5,'low','s');
sys_F1=tf(N1,D1);

[N2,D2]=butter(1,2.5,'low','s');
sys_F2=tf(N2,D2);

% Filtro NOTCH
wc=1.44;
psiz=0.;
psip=1.;
sys_F3=tf([1 2*psiz*wc wc^2],[1 2*psip*wc wc^2]); 

figure 1; margin(sys_F1);
figure 2; margin(sys_F2);
figure 3; margin(sys_F3);


clc; clear;
format long

M_1=0.7;  
M_2=1.5;
a=0.1;
e=0.3;
iota=0.6;
t0=0;
omega=deg2rad(25);
phi0=deg2rad(62);
G=4*(pi)^2;          
M=M_1+M_2;
E=0;
%%Calculation of initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
energy=-G*M/(2*a);
strof=sqrt(G*M*a*(1-e^2));
r=a*(1-e*cos(E));
Vy=strof/r;

xold1=[a,0,0,0,Vy,0];

R_iota=[1 0 0; 0 cos(iota) sin(iota); 0 -sin(iota) cos(iota)];
R_omega=[cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1];
R_total=R_omega'*R_iota';
%%Analytical-Newton
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dt=1/(365*24);
dt=0.00001;            
T=1;                  %2 years period
nstep=fix(T/dt);      %integer of T/dt
xnew1=zeros(nstep,6);


for h=1:1:nstep
    t=dt*h;

Ean=zeros(1,10);
Ean(1)=sqrt(G*M/a^3)*(t-t0);

for i=1:1:9
    Ean(i+1)=Ean(i)-(Ean(i)-e*sin(Ean(i))-(sqrt(G*M/a^3))*(t-t0))/(1-e*cos(Ean(i)));
end

Ean=mod(Ean(10),2*pi);
radius=a*(1-e*cos(Ean));
cosf2=(((a*(1-e^2))/radius)-1)/e;

if Ean<=pi
    f2=acos(cosf2);
elseif Ean>pi
    f2=2*pi-acos(cosf2);
end

phi2=mod(phi0+f2,2*pi);

xnew1(h,1)=radius*cos(phi2);
xnew1(h,2)=radius*sin(phi2);

stroformi=sqrt(G*M*a*(1-e^2));
Vt=stroformi/radius;

if Ean<=pi && Ean>=0
    Vr=sqrt(-(stroformi^2)/radius^2-(G*M)/a+(2*G*M)/radius);
elseif Ean>pi && Ean<=2*pi 
    Vr=-sqrt(-(stroformi^2)/radius^2-(G*M)/a+(2*G*M)/radius);
end
xnew1(h,4)=Vr*xnew1(h,1)/radius-Vt*xnew1(h,2)/radius;
xnew1(h,5)=Vr*xnew1(h,2)/radius+Vt*xnew1(h,1)/radius;

end 

xnew1_t=xnew1';
xnew1_tx=[xnew1_t(1,:);xnew1_t(2,:);xnew1_t(3,:)];
xnew1_tv=[xnew1_t(4,:);xnew1_t(5,:);xnew1_t(6,:)];

ksi1=R_total*xnew1_tx;
ksi1v=R_total*xnew1_tv;

r_m1=M_2*ksi1/M;
r_m2=-M_1*ksi1/M;
v_m1=M_2*ksi1v/M;
v_m2=-M_1*ksi1v/M;

vz_m1=v_m1(3,:)*1.5*10^8/31557600; %in km/sec
vz_m2=v_m2(3,:)*1.5*10^8/31557600;


figure(1)
plot(xnew1(:,4),xnew1(:,5),'.b')
grid on
xlabel('Ux (AU/yr)')
ylabel('Uy (AU/yr)')
title('Ηλιοκεντρικό σύστημα στο επίπεδο της έλλειψης - Ταχύτητες')
hold on

figure(2)
plot(xnew1(:,1),xnew1(:,2),'.b')
grid on
xlabel('x (AU)')
ylabel('y (AU)')
title('Ηλιοκεντρικό σύστημα στο επίπεδο της έλλειψης - Τροχιά')


figure(3)
plot3(ksi1v(1,:),ksi1v(2,:),ksi1v(3,:),'.b')
grid on
xlabel('Ux (AU/yr)')
ylabel('Uy (AU/yr)') 
zlabel('Uz (AU/yr)')
title('Παρατηρούμενες ταχύτητες στο ηλιοκεντρικό σύστημα')


figure(4)
plot3(ksi1(1,:),ksi1(2,:),ksi1(3,:),'.b')
grid on
xlabel('x (AU)')
ylabel('y (AU)')
zlabel('z (AU)')
title('Παρατηρούμενη τροχιά στο ηλιοκεντρικό σύστημα')

figure(5)
plot3(v_m1(1,:),v_m1(2,:),v_m1(3,:),'.b')
grid on
xlabel('Uξ (AU/yr)')
ylabel('Uη (AU/yr)')
zlabel('Uζ (AU/yr)')
title('Παρατηρούμενες ταχύτητες σωμάτων')
hold on
plot3(v_m2(1,:),v_m2(2,:),v_m2(3,:),'.r')
legend('M1','M2')
hold off


figure(6)
plot3(r_m1(1,:),r_m1(2,:),r_m1(3,:),'.b')
grid on
xlabel('ξ (AU)')
ylabel('η (AU)')
zlabel('ζ (AU)')
title('Παρατηρούμενες τροχιές σωμάτων')
hold on
plot3(r_m2(1,:),r_m2(2,:),r_m2(3,:),'.r')
legend('M1','M2')
hold off

figure(7)
plot(r_m1(1,:),r_m1(2,:),'.b')
grid on
xlabel('ξ')
ylabel('η')
title('...')
hold on
plot(r_m2(1,:),r_m2(2,:),'.r')
legend('M1','M2')
hold off

figure(8)
plot(v_m1(1,:),v_m1(2,:),'.b')
grid on
xlabel('vξ')
ylabel('vη')
title('...')
hold on
plot(v_m2(1,:),v_m2(2,:),'.r')
legend('M1','M2')
hold off


time=zeros(1,100000);
time(1,1)=0;
for k=2:1:100000
    time(1,k)=time(1,k-1)+0.00365;
end
figure(9)
plot(time(1:7800),vz_m1(1:7800),'.b')
hold on
plot(time(1:7800),vz_m2(1:7800),'.r')
xlabel('t(days)')
ylabel('Vζ(km/sec)')
title('Καμπύλη φασματοσκοπικού συστήματος')
legend('M1','M2')
grid on

%vz_cm=(vz_m1(1,1)+vz_m2(1,1))/2

%figure(10)
%plot(time(1:7800),vz_m1(1:7800)-10,'.b')
%hold on
%plot(time(1:7800),vz_m2(1:7800)+10,'.r')
%xlabel('t(days)')
%ylabel('Vζ(km/sec)')
%title('Καμπύλη φασματοσκοπικού συστήματος')
%legend('M1','M2')
%grid on

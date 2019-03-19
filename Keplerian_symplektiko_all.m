clc; clear;
format long

R=1;                  %earth-sun distance 
M=1;                  %sun mass in Mo
G=4*(pi)^2;         %gravitational constant in AU^3/(yr^2*Mo)
dt=0.0001;            
T=2;                  %2 years period
nstep=fix(T/dt);      %integer of T/dt
xold=[1,-0.945,0,1.023,2*pi-0.8546,0]; 
v=zeros(1,nstep);
r=zeros(1,nstep);
En=zeros(1,nstep);
time=zeros(1,nstep);
res=zeros(nstep,6);

%Ellipse characteristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Energy=(xold(4)^2+xold(5)^2)/2-G*M/sqrt(xold(1)^2+xold(2)^2);
if Energy<0 
    a=-G*M/(2*Energy);
else fprintf('Error E>0')
    return; %stop running
end

l=xold(1)*xold(5)-xold(2)*xold(4); %angular momentum
e=sqrt(1-l^2/(G*M*a));

cosf=((a*(1-e^2)/sqrt(xold(1)^2+xold(2)^2))-1)/e;

if (xold(1)*xold(4)+xold(2)*xold(5))/sqrt(xold(1)^2+xold(2)^2)>0
    f=acos(cosf); 
elseif (xold(1)*xold(4)+xold(2)*xold(5))/sqrt(xold(1)^2+xold(2)^2)<=0
    f=2*pi-acos(cosf);  
end

if xold(2)>=0
    phi=acos(xold(1)/sqrt(xold(1)^2+xold(2)^2)); 
elseif xold(2)<0 
    phi=2*pi-acos(xold(1)/sqrt(xold(1)^2+xold(2)^2)); 
end
phi0=mod(phi-f,2*pi); %in rad

cosE=(1-(sqrt(xold(1)^2+xold(2)^2)/a))/e;

if xold(2)>=0 
    E=acos(cosE);
elseif xold(2)<0
    E=2*pi-acos(cosE);
end
t0=(e*sin(E)-E)/(sqrt(G*M/(a^3)));

%%%Euler - Arithmetic method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
for i=1:1:nstep 
   v=(xold(4)^2+xold(5)^2+xold(6)^2)^(1/2);
   r=(xold(1)^2+xold(2)^2+xold(3)^2)^(1/2);
   En(1,i)=(v^2)/2-(G*M)/r;
   time(1,i)=i*dt;
    
       xnew(4)=xold(4)-(G*M*xold(1)/((xold(1)^2+xold(2)^2+xold(3)^2)^(3/2)))*dt;
       xnew(5)=xold(5)-(G*M*xold(2)/((xold(1)^2+xold(2)^2+xold(3)^2)^(3/2)))*dt;
       xnew(6)=xold(6)-(G*M*xold(3)/((xold(1)^2+xold(2)^2+xold(3)^2)^(3/2)))*dt;
       xnew(1)=xold(1)+xnew(4)*dt;
       xnew(2)=xold(2)+xnew(5)*dt;
       xnew(3)=xold(3)+xnew(6)*dt;
       res(i,:)=[xnew(1),xnew(2),xnew(3),xnew(4),xnew(5),xnew(6)];
        
  xold=xnew;  
end
figure(1)
plot(res(:,1),res(:,2),'.b')
xlabel('x')
ylabel('y')
title('Keplerian orbit of Earth')
grid on
axis([-1.5,1.5,-1.5,1.5])
hold on
plot(0,0,'*r')
legend('Earth','Sun')

figure(2)
plot(time,En,'.b')
xlabel('t')
ylabel('E')
title('Energy in Earths motion')
grid on

%Newton's Method and trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xnew2=zeros(nstep,6);
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

xnew2(h,1)=radius*cos(phi2);
xnew2(h,2)=radius*sin(phi2);

stroformi=sqrt(G*M*a*(1-e^2));
Vt=stroformi/radius;

if Ean<=pi && Ean>=0 
   Vr=sqrt(-(stroformi^2)/radius^2-(G*M)/a+(2*G*M)/radius);
else
    Vr=-sqrt(-(stroformi^2)/radius^2-(G*M)/a+(2*G*M)/radius);
end

xnew2(h,4)=Vr*xnew2(h,1)/radius-Vt*xnew2(h,2)/radius;
xnew2(h,5)=Vr*xnew2(h,2)/radius+Vt*xnew2(h,1)/radius;

end 
figure(3)
plot(xnew2(:,4),xnew2(:,5),'.b')
grid on
xlabel('Ux')
ylabel('Uy')
title('Keplerian orbit of Earth')
hold on
plot(0,0,'*r')
plot(res(:,4),res(:,5),'--r')
legend('Analytical','Sun','Euler')
hold off

figure(4)
plot(xnew2(:,1),xnew2(:,2),'.b')
grid on
xlabel('x')
ylabel('y')
title('Keplerian orbit of Earth')
hold on
plot(0,0,'*r')
plot(res(:,1),res(:,2),'--r')
legend('Analytical','Sun','Euler')

clc; clear;
format long
R=1;                  %earth-sun distance 
u=2*pi;             %linear velocity in AU/yr
M=1;                  %sun mass in Mo
G=4*(pi)^2;         %gravitational constant in AU^3/(yr^2*Mo)
dt=0.001;            %
T=2;                  %2 years period
nstep=fix(T/dt);      %integer of T/dt
xold=[1,0.043,0,1.023,2*pi-0.8546,0]; %xold=[R,0,0,0,u,0];
v=zeros(1,nstep);
r=zeros(1,nstep);
En=zeros(1,nstep);
time=zeros(1,nstep);
res=zeros(nstep,6);

Energy=(xold(4)^2+xold(5)^2)/2-G*M/sqrt(xold(1)^2+xold(2)^2);
if Energy<0 
    a=-G*M/(2*Energy);
else printf('Error E>0')
end

l=xold(1)*xold(5)-xold(2)*xold(4); %στροφορμή
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



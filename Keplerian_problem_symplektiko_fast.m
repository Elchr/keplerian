clc; clear;
k=1;
R=1;                  %earth-sun distance 
u=2*pi;             %linear velocity in AU/yr
M=1;                  %sun mass in Mo
G=4*(pi)^2;         %gravitational constant in AU^3/(yr^2*Mo)
dE=zeros(8,1);
time2=zeros(8,1);

for dt=[0.001,0.003,0.0001,0.0003,0.00001,0.00003,0.000001,0.000003]
T=sqrt(2);                  %2 years period
nstep=fix(T/dt);            %integer of T/dt
xold=[R,0,0,0,u,0];
v=zeros(1,nstep);
r=zeros(1,nstep);
E=zeros(1,nstep);
time=zeros(1,nstep);

for i=1:1:nstep 
   v=(xold(4)^2+xold(5)^2+xold(6)^2)^(1/2);
   r=(xold(1)^2+xold(2)^2+xold(3)^2)^(1/2);
   E(1,i)=(v^2)/2-(G*M)/r;
   time=i*dt;
    
       xnew(4)=xold(4)-(G*M*xold(1)/((xold(1)^2+xold(2)^2+xold(3)^2)^(3/2)))*dt;
       xnew(5)=xold(5)-(G*M*xold(2)/((xold(1)^2+xold(2)^2+xold(3)^2)^(3/2)))*dt;
       xnew(6)=xold(6)-(G*M*xold(3)/((xold(1)^2+xold(2)^2+xold(3)^2)^(3/2)))*dt;
       xnew(1)=xold(1)+xnew(4)*dt;
       xnew(2)=xold(2)+xnew(5)*dt;
       xnew(3)=xold(3)+xnew(6)*dt;
        
 xold=xnew;  
end

dE(k,1)=abs(E(1,1)-E(1,nstep));
time2(k,1)=dt;
k=k+1;
end

plot(log10(time2),log10(dE),'*')
%set(gca,'xscale','log','yscale','log')
mdl = fitlm(log10(time2),log10(dE))
%hl = lsline;  %least square line
%B = [ones(size(hl.XData(:))), hl.XData(:)]\hl.YData(:);
%Slope = B(2)
%axis([1,10,1,10])
grid on
ylabel('log|ÄE|')
xlabel('logÄt')
title('Óõìðëåêôéêü Ó÷Þìá Euler T=sqrt(2)')

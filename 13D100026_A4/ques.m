clear all;
clc;
close all;
epsilon=0.001;
gridtype= 1;
if gridtype == 1
	m=9; 
	n=4;
	dzi=1/(m-1); 
	dzeta=1/(n-1);
	x_old=ones(n,m); 
	y_old=ones(n,m);
	x=ones(n,m); 
	y=ones(n,m); 
	for i=1:m-1
		x(1,i)=5+2*cos(pi*(i-1)/4);
 		y(1,i)=5+2*sin(2*pi-pi*(i-1)/4)
 	end
	x(n,1)=10; y(n,1)=5;
	x(n,2)=10; y(n,2)=0;
	x(n,3)=5; y(n,3)=0;
	x(n,4)=0; y(n,4)=0;
	x(n,5)=0; y(n,5)=5;
	x(n,6)=0; y(n,6)=10;
	x(n,7)=5; y(n,7)=10;
	x(n,8)=10; y(n,8)=10;
	x(n,9)=10; y(n,9)=5;
	x(1,1)=7; x(2,1)=8; x(3,1)=9; x(4,1)=10; y(1:4,1)=5;
	x(1,m)=7; x(2,m)=8; x(3,m)=9; x(4,m)=10; y(1:4,m)=5;
end

if gridtype== 2
	m=16; 
	n=5;
	dzi=1/(m-1); 
	dzeta=1/(n-1);
	x_old=rand(n,m); 
	y_old=rand(n,m);
	x=rand(n,m); 
	y=rand(n,m); 
	for i=4:13
		x(1,i)=5+2*cos(2*pi*(i-4)/9);
		y(1,i)=5+2*sin(2*pi-2*pi*(i-4)/9);
	end
	x(1:5,1)=10;
	for j=1:5
		y(j,1)=5-10*(j-1)/8;
	end
	x(1:5,m)=10;
	for j=1:5
		y(j,m)=5+10*(j-1)/8;
	end
	y(n,1:6)=0;
	for i=1:6
		x(n,i)=10-10*(i-1)/5;
	end
	x(n,6:11)=0;
	for i=6:11
		y(n,i)=10*(i-6)/5;
	end
	y(n,11:16)=10;
	for i=11:16
		x(n,i)=10*(i-11)/5;
 	end
	y(1,1:4)=5;
	for i=1:4
		x(1,i)=11-i;
	end
	y(1,13:16)=5;
	for i=13:16
		x(1,i)=7+(i-13);
	end
end

if gridtype==3;
	m=11; 
	n=5;
	dzi=1/(m-1); 
	dzeta=1/(n-1);
	x_old=rand(n,m); 
	y_old=rand(n,m);
	x=rand(n,m); 
	y=rand(n,m); 
	for i=4:8
		x(1,i)=5+2*cos(pi*(8-i)/4); 
		y(1,i)=2*sin(pi*(8-i)/4); 
	end
	x(1:n,1)=0;
	for j=1:n
		y(j,1)=5*(j-1)/4;
	end
	y(n,1:m)=5;
	for i=1:m
		x(n,i)=i-1;
	end
	x(1:n,m)=10;
	for j=1:n
		y(j,m)=5*(j-1)/4;
	end
	y(1,1:4)=0;
	for i=1:4
		x(1,i)=i-1;
	end
	y(1,8:11)=0;
	for i=8:11
    		x(1,i)=i-1;
	end
end

Error=1;
while Error>=epsilon
   x_old=x;
   y_old=y;
   for i=2:m-1
       for j=2:n-1
           a(j,i)=((x_old(j+1,i)-x_old(j-1,i))/(2*dzeta))^2+((y_old(j+1,i)-y_old(j-1,i))/(2*dzeta))^2;
           b(j,i)=((x_old(j,i+1)-x_old(j,i-1))*(x_old(j+1,i)-x_old(j-1,i))+(y_old(j,i+1)-y_old(j,i-1))*(y_old(j+1,i)-y_old(j-1,i)))/(4*dzi*dzeta);
           c(j,i)=((x_old(j,i+1)-x_old(j,i-1))/(2*dzi))^2+((y_old(j,i+1)-y_old(j,i-1))/(2*dzi))^2;
           Dx(j,i)=-2*b(j,i)*(x_old(j+1,i+1)+x_old(j-1,i-1)-x_old(j+1,i-1)-x_old(j-1,i+1))/(4*dzi*dzeta);
           Dy(j,i)=-2*b(j,i)*(y_old(j+1,i+1)+y_old(j-1,i-1)-y_old(j+1,i-1)-y_old(j-1,i+1))/(4*dzi*dzeta);
           x(j,i)=(a(j,i)*(x(j,i+1)+x(j,i-1))*(dzeta*dzeta)+c(j,i)*(x(j+1,i)+x(j-1,i))*(dzi*dzi)+Dx(j,i)*dzi*dzi*dzeta*dzeta)/(2*a(j,i)*(dzeta*dzeta)+2*c(j,i)*(dzi*dzi));
           y(j,i)=(a(j,i)*(y(j,i+1)+y(j,i-1))*(dzeta*dzeta)+c(j,i)*(y(j+1,i)+y(j-1,i))*(dzi*dzi)+Dy(j,i)*dzi*dzi*dzeta*dzeta)/(2*a(j,i)*(dzeta*dzeta)+2*c(j,i)*(dzi*dzi));
       end
   end
   Error=max(max(abs(x(j,i)-x_old(j,i)),abs(y(j,i)-y_old(j,i))));
   disp(Error)
end

figure;
if (gridtype==3)
    for k=2:n-1
      plot (x(k,:),y(k,:)+5,'b-o','linewidth',1.5,'markeredgecolor','k','markersize',10,'markerfacecolor','y'); hold on;
    end
    for k=2:m-1
     plot (x(:,k),y(:,k)+5,'b-o','linewidth',1.5,'markeredgecolor','k','markersize',10,'markerfacecolor','y'); hold on;
    end
    plot (x(1,:),((-1)*y(1,:))+5,'b-','linewidth',1.5,'markeredgecolor','k','markersize',10,'markerfacecolor','y'); hold on;
    for k=2:n-1
     plot (x(k,:),((-1)*y(k,:))+5,'b-','linewidth',1.5,'markeredgecolor','k','markersize',10,'markerfacecolor','y'); hold on;
    end
    for k=2:m-1
     plot (x(:,k),((-1)*y(:,k))+5,'b-','linewidth',1.5,'markeredgecolor','k','markersize',10,'markerfacecolor','y'); hold on;
    end
    plot (x(1,:),y(1,:)+5,'k-o','linewidth',2,'markeredgecolor','k','markersize',10,'markerfacecolor','r'); 
    hold on;
    plot (x(n,:),y(n,:)+5,'k-o','linewidth',2,'markeredgecolor','k','markersize',10,'markerfacecolor','r'); 
    hold on;
    plot (x(:,1),y(:,1)+5,'k-o','linewidth',2,'markeredgecolor','k','markersize',10,'markerfacecolor','r'); 
    hold on;
    plot (x(:,m),y(:,m)+5,'k-o','linewidth',2,'markeredgecolor','k','markersize',10,'markerfacecolor','r'); 
    hold on;
    hold off;
else
    for k=2:n-1
    plot (x(k,:),y(k,:),'b-o','linewidth',1.5,'markeredgecolor','k','markersize',10,'markerfacecolor','y'); 
    hold on;
    end
    for k=2:m-1
    plot (x(:,k),y(:,k),'b-o','linewidth',1.5,'markeredgecolor','k','markersize',10,'markerfacecolor','y'); 
    hold on;
    end
    plot (x(1,:),y(1,:),'k-o','linewidth',2,'markeredgecolor','k','markersize',10,'markerfacecolor','r'); 
    hold on;
    plot (x(n,:),y(n,:),'k-o','linewidth',2,'markeredgecolor','k','markersize',10,'markerfacecolor','r'); 
    hold on;
    plot (x(:,1),y(:,1),'k-o','linewidth',2,'markeredgecolor','k','markersize',10,'markerfacecolor','r'); 
    hold on;
    plot (x(:,m),y(:,m),'k-o','linewidth',2,'markeredgecolor','k','markersize',10,'markerfacecolor','r'); 
    hold on;
    hold off; 
end
L=10; R=2; 
rho=7750; cp=500;  k=16.2; 
T_in=1; T_out=0;
dt=400; 
Q_gen=0; 
n=9; 
m=4; 
epsilon=0.000001; 
dzeta=1/(n-1);  
dzi=1/(m-1);
x_f=x; yf=y;
for i=2:m
    for j=2:n
        xp(i,j)=(x(i-1,j)+x(i-1,j-1)+x(i,j-1)+x(i,j))/4;
        yp(i,j)=(y(i-1,j)+y(i-1,j-1)+y(i,j-1)+y(i,j))/4;
    end
end
for j=2:n
    xp(1,j)=(x(1,j-1)+x(1,j))/2;
    yp(1,j)=(y(1,j-1)+y(1,j))/2;
end
for j=2:n
    xp(m+1,j)=(x(m,j-1)+x(m,j))/2;
    yp(m+1,j)=(y(m,j-1)+y(m,j))/2;
end
xp(:,1)=xp(:,n); 
yp(:,1)=yp(:,n);
xp(:,10)=xp(:,2); 
yp(:,10)=yp(:,2);

T=zeros(m+1,n+1);
T(1,:)=T_in; 
T(m+1,:)=T_out;
Told=T;
Tnew=Told;
m=4; n=9;
% Vertical faces Parameters.
dzi_x=zeros(m,n); 
dzeta_y=zeros(m,n);
dse_zi=zeros(m,n); 
dse_zeta=zeros(m,n);
for i=2:m
    for j=1:n
        dzi_xx_=xp(i,j+1)-xp(i,j); 
        dzi_xy_=yp(i,j+1)-yp(i,j);
        dzi_x(i,j)=sqrt(dzi_xx_^2+dzi_xy_^2);
        zi1_e=dzi_xx_/dzi_x(i,j); 
        zi2_e=dzi_xy_/dzi_x(i,j);
        dzeta_yx=x_f(i,j)-x_f(i-1,j);
        dzeta_yy=yf(i,j)-yf(i-1,j); 
        dzeta_y(i,j)=sqrt(dzeta_yx^2+dzeta_yy^2);
        zeta1_e=dzeta_yx/dzeta_y(i,j);
        zeta2_e=dzeta_yy/dzeta_y(i,j); 
        dse_y=-dzeta_yx;
        dse_x=dzeta_yy;
        dse_zi(i,j)=(-dse_y*zeta1_e+dse_x*zeta2_e)/(zi1_e*zeta2_e-zi2_e*zeta1_e);
        dse_zeta(i,j)=(dse_y*zi1_e-dse_x*zi2_e)/(zi1_e*zeta2_e-zi2_e*zeta1_e);
    end
end

dzi_n=zeros(m,n); 
dzeta_n=zeros(m,n); 
dsn_zi=zeros(m,n);
dsn_zeta=zeros(m,n);
for i=1:m
    for j=2:n
        dzi_nx=xp(i+1,j)-xp(i,j); 
        dzi_ny=yp(i+1,j)-yp(i,j);
        dzi_n(i,j)=sqrt(dzi_nx^2+dzi_ny^2);
        zi1_n=dzi_nx/dzi_n(i,j); 
        zi2_n=dzi_ny/dzi_n(i,j);
        dzeta_nx=x_f(i,j)-x_f(i,j-1); 
        dzeta_ny=yf(i,j)-yf(i,j-1); 
        dzeta_n(i,j)=sqrt(dzeta_nx^2+dzeta_ny^2);
        zeta1_n=dzeta_nx/dzeta_n(i,j); 
        zeta2_n=dzeta_ny/dzeta_n(i,j);
        dsn_y=dzeta_nx; 
        dsn_x=-dzeta_ny; 
        dsn_zi(i,j)=(-dsn_y*zeta1_n+dsn_x*zeta2_n)/(zi1_n*zeta2_n-zi2_n*zeta1_n);
        dsn_zeta(i,j)=(dsn_y*zi1_n-dsn_x*zi2_n)/(zi1_n*zeta2_n-zi2_n*zeta1_n);
    end
end
% Volume.
dv=zeros(m+1,n+1);
for i=2:m
    for j=2:n
        d1=sqrt((x_f(i-1,j)-x_f(i,j-1))^2+(yf(i-1,j)-yf(i,j-1))^2);
        d2=sqrt((x_f(i-1,j-1)-x_f(i,j))^2+(yf(i-1,j-1)-yf(i,j))^2);
        dv(i,j)=(1/4)*sqrt(4*(d1^2)*(d2^2)-(dzeta_y(i,j)^2+dzeta_y(i,j-1)^2-dzeta_n(i,j)^2-dzeta_n(i-1,j)^2));
    end
end
dv(:,1)=dv(:,9);dv(:,10)=dv(:,2);
error=1; count=0;
while error>=epsilon
for i=2:m-1
    for j=1:n
        T_c(i,j)=(dv(i,j+1)*Told(i,j+1)+dv(i,j)*Told(i,j)+dv(i+1,j)*Told(i+1,j)+dv(i+1,j+1)*Told(i+1,j+1))/(dv(i,j+1)+dv(i,j)+dv(i+1,j)+dv(i+1,j+1));
    end
end
T_c(1,:)=T_in; T_c(m,:)=T_out;
for i=2:m
    for j=2:n
        Q1_zeta=k*(((Told(i,j+1)-Told(i,j))*dse_zi(i,j)/(dzi_x(i,j)))+((Told(i+1,j)-Told(i,j))*dsn_zi(i,j)/(dzi_n(i,j)))-((Told(i,j)-Told(i,j-1))*dse_zi(i,j-1)/(dzi_x(i,j-1)))-((Told(i,j)-Told(i-1,j))*dsn_zi(i-1,j)/(dzi_n(i-1,j))));
        Q1_zi=k*(((T_c(i,j)-T_c(i-1,j))*dse_zeta(i,j)/dzeta_y(i,j))+((T_c(i,j)-T_c(i,j-1))*dsn_zeta(i,j)/dzeta_n(i,j))-((T_c(i,j-1)-T_c(i-1,j-1))*dse_zeta(i,j-1)/dzeta_y(i,j-1))-((T_c(i-1,j)-T_c(i-1,j-1))*dsn_zeta(i-1,j)/dzeta_n(i-1,j)));
        Q_cond=Q1_zeta+Q1_zi;
        Q_gen=Q_gen*dv(i,j);
        Tnew(i,j)=Told(i,j)+(dt/(rho*dv(i,j)*cp))*(Q_cond+Q_gen);
    end
end 
Tnew(:,1)=Tnew(:,n); Tnew(:,n+1)=Tnew(:,2);
error=max(max(Tnew-Told));

Told=Tnew;
end

for i=2:m-1
    for j=1:n
        T_c(i,j)=(dv(i,j+1)*Told(i,j+1)+dv(i,j)*Told(i,j)+dv(i+1,j)*Told(i+1,j)+dv(i+1,j+1)*Told(i+1,j+1))/(dv(i,j+1)+dv(i,j)+dv(i+1,j)+dv(i+1,j+1));
    end
end
T_c(1,:)=T_in; T_c(m,:)=T_out;
figure;
v=linspace(0,1,11);
contourf(x_f,yf,T_c);
[c,h]=contourf(x_f,yf,T_c,v);
clabel(c,h);
colorbar;

Q1=0;
for j=2:n
    Q1_zeta=-k*((Tnew(2,j)-Tnew(1,j))*dsn_zi(1,j)/(dzi_n(1,j)));
    Q1_zi=-k*((T_c(1,j)-T_c(1,j-1))*dsn_zeta(1,j)/dzeta_n(1,j));
    Q1=Q1+Q1_zeta+Q1_zi;
end
disp(Q1);
Q2=0;
for j=2:n
    Q1_zeta=-k*(Tnew(m+1,j)-Tnew(m,j))*dsn_zi(m,j)/(dzi_n(m,j));
    Q1_zi=-k*((T_c(m,j)-T_c(m,j-1))*dsn_zeta(m,j)/dzeta_n(m,j));
    Q2=Q2+Q1_zeta+Q1_zi;
end
disp(Q2);

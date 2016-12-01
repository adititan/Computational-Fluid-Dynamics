clear
clc
rho = 1;
uo = 1;
L=1; 
Re=400; 
mu=1/Re;
imax=42; 
jmax=42;
epsilon_st= 0.0001;
epsilon = 0.00000001
dx = L/(imax-2); 
dy = L/(jmax-2); 
dt=0.001;

u=zeros(jmax,imax-1);  
v=zeros(jmax-1,imax);  
p=zeros(jmax,imax);  
p_correc=zeros(jmax,imax); 
u_approx=zeros(jmax,imax-1);  
v_approx=zeros(jmax-1,imax);

u_approx(jmax,:)=uo; 

e_c=(dt*rho*dy)/dx;  
w_c=(dt*rho*dy)/dx;
n_c=(dt*rho*dx)/dy;  
s_c=(dt*rho*dx)/dy;
ap=e_c+w_c+n_c+s_c;

unsteadiness_nd=1; 

while unsteadiness_nd>=epsilon_st
    u_old=u_approx;
    v_old=v_approx;
    
    for j=2:jmax-1                                
       for i=1:imax-2 
           m_ux(j,i)=((u_old(j,i)+u_old(j,i+1))*rho)/2;
     end 
    end
    
    for j=2:jmax-1                                
       for i=1:imax-2 
           a_ux(j,i)=(max(m_ux(j,i),0)*u_old(j,i))+(min(m_ux(j,i),0)*u_old(j,i+1));
     end 
    end
    
    for j=2:jmax-1                                
       for i=1:imax-2 
         d_ux(j,i)=((u_old(j,i+1)-u_old(j,i))/dx)*mu;
       end 
    end
    
     for j=1:jmax-1                               
       for i=2:imax-2 
                m_uy(j,i)=((v_old(j,i)+v_old(j,i+1))*rho)/2;
       end 
     end 
    
     
     for j=1:jmax-1                               
       for i=2:imax-2 
            a_uy(j,i)=(max(m_uy(j,i),0)*u_old(j,i))+(min(m_uy(j,i),0)*u_old(j+1,i));
       end 
     end 
    
    for j=1:jmax-1                               
       for i=2:imax-2 
         d_uy(j,i)=((u_old(j+1,i)-u_old(j,i))/dy)*mu;
       end 
    end
    
    for j=2:jmax-2                                
       for i=1:imax-1 
           m_vx(j,i)=((u_old(j,i)+u_old(j+1,i))*rho)/2;
       end 
    end
    
     for j=2:jmax-2                                
       for i=1:imax-1 
           a_vx(j,i)=(max(m_vx(j,i),0)*v_old(j,i))+(min(m_vx(j,i),0)*v_old(j,i+1));
       end 
    end
           
    for j=2:jmax-2                                
       for i=1:imax-1   
         d_vx(j,i)=((v_old(j,i+1)-v_old(j,i))/dx)*mu;
       end 
    end
     
    for j=1:jmax-2                              
       for i=2:imax-1
            m_vy(j,i)=((v_old(j,i)+v_old(j+1,i))*rho)/2;
       end 
    end  
    
    for j=1:jmax-2                              
       for i=2:imax-1 
           a_vy(j,i)=(max(m_vy(j,i),0)*v_old(j,i))+(min(m_vy(j,i),0)*v_old(j+1,i));
       end 
    end
    for j=1:jmax-2                              
       for i=2:imax-1 
         d_vy(j,i)=((v_old(j+1,i)-v_old(j,i))/dy)*mu;
       end 
    end
    
    for j=2:jmax-1                               
       for i=2:imax-2 
         Ax(j,i)=((a_ux(j,i)-a_ux(j,i-1))*dy)+((a_uy(j,i)-a_uy(j-1,i))*dx);
         Dx(j,i)=((d_ux(j,i)-d_ux(j,i-1))*dy)+((d_uy(j,i)-d_uy(j-1,i))*dx);
         Sx(j,i)=(p(j,i)-p(j,i+1))*dy;
         u_approx(j,i)=u_old(j,i)+((dt/(rho*dx*dy))*(Dx(j,i)-Ax(j,i)+Sx(j,i))) ;
       end 
    end
    
    for j=2:jmax-2                              
       for i=2:imax-1 
         Ay(j,i)=((a_vx(j,i)-a_vx(j,i-1))*dy)+((a_vy(j,i)-a_vy(j-1,i))*dx);
         Dy(j,i)=((d_vx(j,i)-d_vx(j,i-1))*dy)+((d_vy(j,i)-d_vy(j-1,i))*dx);
         Sy(j,i)=(p(j,i)-p(j+1,i))*dx;
         v_approx(j,i)=v_old(j,i)+((dt/(rho*dx*dy))*(Dy(j,i)-Ay(j,i)+Sy(j,i))) ;
       end 
    end
         for j=2:jmax-1
    for i=2:imax-1
        S_mp(j,i)=((u_approx(j,i)-u_approx(j,i-1))*dy)+((v_approx(j,i)-v_approx(j-1,i))*dx);
    end    
  end  

 while (max(max(abs(S_mp(j,i)))) > epsilon)      
  for j=2:jmax-1
    for i=2:imax-1
        S_mp(j,i)=((u_approx(j,i)-u_approx(j,i-1))*dy)+((v_approx(j,i)-v_approx(j-1,i))*dx);
    end    
  end  

  
      Error=1;
    while (Error>=epsilon_st)
    P_old=p_correc; 
     for j=2:jmax-1
       for i=2:imax-1
        p_correc(j,i)=((e_c*p_correc(j,i+1))+(w_c*p_correc(j,i-1))+(n_c*p_correc(j+1,i))+(s_c*p_correc(j-1,i))-(S_mp(j,i)))/ap;
       end    
     end 
        p_correc(:,imax)= p_correc(:,imax-1); 
        p_correc(jmax,:)= p_correc(jmax-1,:);
        p_correc(:,1)= p_correc(:,2); 
        p_correc(1,:)= p_correc(2,:); 
        Error=max(max(abs(p_correc-P_old)));
    end
    for j=2:jmax-1                              
       for i=2:imax-2 
          u_approx(j,i)=u_approx(j,i)+ ((dt/(rho*dx))*(p_correc(j,i)-p_correc(j,i+1))) ;
       end 
     end
    
    for j=2:jmax-2                              
       for i=2:imax-1 
          v_approx(j,i)=v_approx(j,i) + ((dt/(rho*dy))*(p_correc(j,i)-p_correc(j+1,i))) ;
       end 
    end  
      
 end 
for j=1:jmax
    for i=1:imax
        p(j,i) = p(j,i) + p_correc(j,i);
    end
end 
u_approx(:,imax-1)=0;  
v_approx(jmax-1,:)=0; 
unsteadiness_nd=max((max(max(abs(u_approx-u_old)))),(max(max(abs(v_approx-v_old)))));  
unsteadiness_nd 
end


x=linspace(0,L,imax-1);
y=linspace(0,L,jmax-1); 


xc(1)=x(1); xc(imax)=x(imax-1); 
for i=2:imax-1
xc(i)=(x(i)+x(i-1))/2;
end
yc(1)=y(1); yc(jmax)=y(jmax-1); 
for i=2:jmax-1
yc(i)=(y(i)+y(i-1))/2;
end 
 
v=linspace(-1,1,51); 
w=linspace(-1,1,26);

figure; 

[C,h]=contourf(x,yc,u_approx,w); 
clabel(C,h);  
title('U Velocity Contour Plot');

figure;
[C,h]=contourf(xc,yc,p,w); 
clabel(C,h); 
title('Pressure Contour Plot');

u_approx(:,imax)=0;
v_approx(jmax,:)=0;
figure;
quiver(xc,yc,u_approx,v_approx,1.25); 
title('U velocity vector plot'); 


%benchmark Results
yp = [1 0.9766 0.9688 0.9609 0.9531 0.8516 0.7344 0.6172 0.5 0.4531 0.2813 0.1719 0.1016 0.0703 0.0625 0.0547 0];          
uy_100 = [1 0.84123 0.78871 0.73722 0.68717 0.23151 0.00332 -0.13641 -0.20581 -0.2109 -0.15662 -0.1015 -0.06434 -0.04775 -0.04192 -0.03717 0];
vx_100 = [0 -0.05906 -0.07391 -0.08864 -0.10313 -0.16914 -0.22445 -0.24533 0.05454 0.17527 0.17507 0.16077 0.12317 0.1089 0.100091 0.09233 0];
xp = [1 0.9688 0.9609 0.9531 0.9453 0.9063 0.8594 0.8047 0.5 0.2344 0.2266 0.1563 0.0938 0.0781 0.0703 0.0625 0];

figure;
u_mid = (u_approx(:,imax/2)+u_approx(:,(imax/2)+1))/2
plot(u_mid,yc,'b-s'); 
xlabel('U Velocity','FontWeight','bold'); 
ylabel('Y','FontWeight','bold');  
title('Variation of U-Velocity along the Verticle Centerline'); 
grid on; 
hold on;
if (Re==100)
    plot(uy_100,yp,'k-o'); hold off;
    legend('using Code','Benchmark') ;
end
figure;
v_mid = (v_approx(jmax/2,:)+v_approx((jmax/2)+1,:))/2
plot(xc,v_mid,'r-*'); 
xlabel('X','FontWeight','bold'); 
ylabel('V','FontWeight','bold'); 
title('Variation of V-Velocity along the Horizontal Centerline'); 
grid on; 
hold on;

if (Re==100) 
    plot(xp,vx_100,'k-o'); hold off;
    legend('using Code','Benchmark') ;
end
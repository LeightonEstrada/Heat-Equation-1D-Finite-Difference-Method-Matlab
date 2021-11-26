clear all;clc;
%animacion 1d sobre la ecuacion del calor
%leighton estrada rayme
%04140403 fcm cc 14.5
%---------------------------------------------
%entrada de datos
disp('ecuacion calor u_t=k u_xx; time L=2')
disp('ingrese numero de particiones espaciales')
%n=input('n=');
n=500;
disp('ingrese numero de particiones temporales')
%m=input('m=');
m=500;
hx=2/n;ht=2/m;
f=inline('12*x-26*x.^2+20*x.^3-5*x.^4');
k=1;r=(k*ht)/(2*hx^2);
%---------------------------------------------
%condiciones iniciales
for i=1:m+1 u(1,i)=f((i-1)*ht);end
%---------------------------------------------
%condiciones de frontera
for i=1:n+1 u(i,1)=0;u(i,m+1)=0;end
%---------------------------------------------
%metodo de crank nicholson
for i=1:n-1 A(i,i)=1+2*r;B(i,i)=1-2*r;end
for i=1:n-2 A(i,i+1)=-r;B(i,i+1)=r;end
for i=2:n-1 A(i,i-1)=-r;B(i,i-1)=r;end
for i=2:n+1
    w3=[r*u(i-1,1);zeros(m-3,1);r*u(i-1,m+1)];
    w4=[-r*u(i,1);zeros(m-3,1);-r*u(i,m+1)];
    b=B*u(i-1,2:m)'+w4-w3;
    u(i,2:m)=(A\b)';
end
%---------------------------------------------
%ploteado 3d
for i=1:n+1 for j=1:m+1 t(i,j)=(j-1)*ht;end;end
for j=1:n+1 for i=1:m+1 x(i,j)=(i-1)*hx;end;end
%surf(t,x,u)
%fin
%---------------------------------------------
%animacion 3d
for i=1:length(x)
    h=sprintf('%5.2f',t(1,i));
    plot(t(i,:),u(i,:));
    title({'Metodo Crank-Nicholson 1D';'Hecho por Leighton Estrada R. UNMSM FCM Computacion Cientifica';['Tiempo t=',num2str(h),' seg']});%grid
    xlabel('t');ylabel('u');
    grid
    axis([0 2 0 2])
    M(i)=getframe;
end

clc
clear
format long
k=50; %原模型阶数，即电路结点数，在二次线性化之后模型的阶数为2*N阶
n=2*k;%二次双线性系统的阶数
%降阶模型的阶数大小
nr=2;
%选取插值点的调整系数
L=1;
h = L/(k+1);

%控制最后取值边界
w=0.01; 

%构建单位矩阵
I=eye(n,n);


A=zeros(n,n);
% left boundary
A(1,1) = 1-2/h^2;
A(1,2) = 1/h^2;
A(k+1,k+1) = 2-4/h^2;

% right boundary 
A(k,k-1) = 1/h^2;
A(k,k) = 1-1/h^2;
A(2*k,2*k) = 2-2/h^2;

for i = 2:k-1
   A(i,i-1) = 1/h^2;
   A(i,i) = 1-2/h^2;
   A(i,i+1) = 1/h^2;
   
   A(k+i,k+i) = 2-4/h^2;
end
%A=A-0.0095*eye(n);
%P=real(eig(A));

%构建H矩阵
H=zeros(n,n^2);
for i=1:n
    H(1,(i-1)*n+1)=10;
    %H(:,(i-1)*n+1:i*n)=eye(n);
    for j=2:n
      
     %H(i,(j-1)*n+i)=1+0.1*i;
     H(j,(i-1)*n+j)=1+(-1)^j*0.1*j;
     H(j-1,(i-1)*n+j)=0.1;
     H(j,(i-1)*n+j-1)=0.2;

    end
end
H3=H;

N=zeros(n,n^2);

%
for i=1:n
    N(k+1,1) = 2/h^2;
    %N(:,(i-1)*n+1:i*n)=eye(n);
end
gam=1;
N=gam*N;
N3=N;

B=sparse(n,1);
B(1,1)=1;B(2,1)=1;B(k+1,1)=40;B(k+2,1)=40;

c=zeros(n,1);
c(1,1)=1;

In=linspace(1,1,n);
In=In.';
a0=zeros(n,1);
tmax=5;
tspan=0:0.02:tmax;% 0:0.4:40;%0:0.01:2

% [t x]=ode15s(@(t,x) A*x+H1*kron(x,x)+N*x*(1/gam)*cos(5*t)*sin(0.2*pi*t)+B*(1/gam)*cos(5*t)*sin(0.2*pi*t),0:0.05:5,a0);



%[t x]=ode15s(@(t,x) A*x+H1*kron(x,x)+N*x*(cos(0.5*pi*t)+0.1)+B*(cos(0.5*pi*t)+0.1),0:0.02:15,a0);
%[t x]=ode15s(@(t,x) A*x+H3*kron(x,x)+N*kron(20*(0.1*exp(-t)+0.1)*In,x)+B*20*(0.1*exp(-t)+0.1),0:0.02:tmax,a0);
[t x]=ode15s(@(t,x) A*x+H3*kron(x,x)+N*kron((cos(0.5*pi*t)+0.1)*In,x)+B*(cos(0.5*pi*t)+0.1),0:0.02:tmax,a0);
%[t x]=ode15s(@(t,x) A*x+H3*kron(x,x)+B*(cos(0.5*pi*t)+0.1),0:0.02:15,a0);
y=gam*x*c;

plot(tspan,y,'-k','LineWidth',2)
hold on
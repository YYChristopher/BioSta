theta_A_ori=0.6;
theta_B_ori=0.5;
pi_A_ori=0.5;
pi_B_ori=0.5;%设定初始的参数空间中的4个参数
E_last=1;%设定初始的期望值
error=0.000001;%设定迭代终止的误差值
n=1000;%设定选择一次硬币后抛掷硬币的次数
k=1000;%设定选择硬币的次数，即实验组数(获得的结果矩阵应该为k*n阶矩阵)
iter=1;%设定初始迭代次数为1
result=zeros(k,n);%设定初始的结果矩阵全为0.记硬币正面朝上为1，反面朝上为0.
choice=zeros(k,1);%设定选择硬币的向量，记1为选择硬币A，0为选择硬币B.
num=zeros(1,k);%设定表示每组实验正面朝上次数矩阵，即抛n次实验多少次正面朝上。
max=100000;
gamma_A=zeros(1,k);
gamma_B=zeros(1,k);
theta_A=zeros(1,max);
theta_B=zeros(1,max);
pi_A=zeros(1,max);
pi_B=zeros(1,max);
theta_A(1,1)=theta_A_ori;
theta_B(1,1)=theta_B_ori;
pi_A(1,1)=pi_A_ori;
pi_B(1,1)=pi_B_ori;
for i=1:k
    choice(i,1)=binornd(1,0.50001,1,1);
end%设定实际上会选到硬币A的概率为0.6，即pi_A=0.6;对应的，pi_B=0.4.
for i=1:k
    if(choice(i,1)==1)
        result(i,1:n)=binornd(1,0.8,1,n);
    else
        result(i,1:n)=binornd(1,0.6,1,n);%在这里设定硬币A实际正面朝上概率为0.7，即theta_A=0.7.相应地，theta_B=0.45.
    end
end
for i=1:k
    for j=1:n
        if(result(i,j)==1)
            num(1,i)=num(1,i)+1;
        end
    end
end
stop=0;
for t=2:max
   for i=1:k
       gamma_A(1,i)=pi_A(1,t-1)*(theta_A(1,t-1)^(num(1,i)))*((1-theta_A(1,t-1))^(n-num(1,i)))/(pi_A(1,t-1)*(theta_A(1,t-1)^(num(1,i)))*((1-theta_A(1,t-1))^(n-num(1,i)))+pi_B(1,t-1)*(theta_B(1,t-1)^(num(1,i)))*((1-theta_B(1,t-1))^(n-num(1,i))));
       gamma_B(1,i)=pi_B(1,t-1)*(theta_B(1,t-1)^(num(1,i)))*((1-theta_B(1,t-1))^(n-num(1,i)))/(pi_A(1,t-1)*(theta_A(1,t-1)^(num(1,i)))*((1-theta_A(1,t-1))^(n-num(1,i)))+pi_B(1,t-1)*(theta_B(1,t-1)^(num(1,i)))*((1-theta_B(1,t-1))^(n-num(1,i))));
   end
   theta_A(1,t)=sum(num.*gamma_A)/(n*sum(gamma_A)); 
   theta_B(1,t)=sum(num.*gamma_B)/(n*sum(gamma_B));
   pi_A(1,t)=sum(gamma_A)/k;
   pi_B(1,t)=sum(gamma_B)/k;
   diff=abs(theta_A(1,t)-theta_A(1,t-1))+abs(theta_B(1,t)-theta_B(1,t-1))+abs(pi_A(1,t)-pi_A(1,t-1));
   if(diff<=error)
       stop=iter+1;
       break;
   else
       iter=iter+1;
   end
end%迭代EM算法
figure(1);
plot(1:stop,theta_A(1,1:stop),'r-','linewidth',2);
hold on;
plot(1:stop,theta_B(1,1:stop),'b-','linewidth',2);
hold on;
plot(1:stop,pi_A(1,1:stop),'m-','linewidth',2);
hold on;
plot(1:stop,pi_B(1,1:stop),'g-','linewidth',2);
hold on;
axis([1 stop 0 1]);
legend('theta_A','theta_B','pi_A','pi_B');
theta_A_ori=0.6;
theta_B_ori=0.5;
pi_A_ori=0.5;
pi_B_ori=0.5;%�趨��ʼ�Ĳ����ռ��е�4������
E_last=1;%�趨��ʼ������ֵ
error=0.000001;%�趨������ֹ�����ֵ
n=1000;%�趨ѡ��һ��Ӳ�Һ�����Ӳ�ҵĴ���
k=1000;%�趨ѡ��Ӳ�ҵĴ�������ʵ������(��õĽ������Ӧ��Ϊk*n�׾���)
iter=1;%�趨��ʼ��������Ϊ1
result=zeros(k,n);%�趨��ʼ�Ľ������ȫΪ0.��Ӳ�����泯��Ϊ1�����泯��Ϊ0.
choice=zeros(k,1);%�趨ѡ��Ӳ�ҵ���������1Ϊѡ��Ӳ��A��0Ϊѡ��Ӳ��B.
num=zeros(1,k);%�趨��ʾÿ��ʵ�����泯�ϴ������󣬼���n��ʵ����ٴ����泯�ϡ�
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
end%�趨ʵ���ϻ�ѡ��Ӳ��A�ĸ���Ϊ0.6����pi_A=0.6;��Ӧ�ģ�pi_B=0.4.
for i=1:k
    if(choice(i,1)==1)
        result(i,1:n)=binornd(1,0.8,1,n);
    else
        result(i,1:n)=binornd(1,0.6,1,n);%�������趨Ӳ��Aʵ�����泯�ϸ���Ϊ0.7����theta_A=0.7.��Ӧ�أ�theta_B=0.45.
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
end%����EM�㷨
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
N=5000;
pi_real=[3/10,5/10,2/10];%设定真实pi值
mu_real=[-2,1;-1,-2;1,1];%设定真实期望
cov_real(:,:,1)=[1,0;0,1];
cov_real(:,:,2)=[3,1;1,3];
cov_real(:,:,3)=[3,1;1,3];%设定真实方差
X_1=mvnrnd(mu_real(1,:),cov_real(:,:,1),N*pi_real(1));
X_2=mvnrnd(mu_real(2,:),cov_real(:,:,2),N*pi_real(2));
X_3=mvnrnd(mu_real(3,:),cov_real(:,:,3),N*pi_real(3));%生成真实数据点集
X=[X_1;X_2;X_3];
X=X(randperm(size(X,1)),:);

x=-10:0.5:10;
y=-10:0.5:10;
[x y]=meshgrid(x,y);
mesh=[x(:),y(:)];
nor_real=pi_real(1)*mvnpdf(mesh,mu_real(1,:),cov_real(:,:,1))+pi_real(2)*mvnpdf(mesh,mu_real(2,:),cov_real(:,:,2))+pi_real(3)*mvnpdf(mesh,mu_real(3,:),cov_real(:,:,3));
figure(1);
plot(X_1(:,1),X_1(:,2),'rx',X_2(:,1),X_2(:,2),'bo',X_3(:,1),X_3(:,2),'g<');
title('混合高斯分布的真实分布');
legend('高斯分布1','高斯分布2','高斯分布3');
hold on;
confiEllipse(X_1,0.95);
hold on;
confiEllipse(X_2,0.95);
hold on;
confiEllipse(X_3,0.95);
hold on;
figure(2);
contour(x,y,reshape(nor_real,size(x,2),size(y,2)));
figure(3);
surf(x,y,reshape(nor_real,size(x,2),size(y,2)));
figure(4);
plot(X(:,1),X(:,2),'kx');
title('未经过EM算法的点');

pi=[1/3,1/3,1/3];%初始化pi的取值
cov(:,:,1)=[1,0;0,1];
cov(:,:,2)=[1,0;0,1];
cov(:,:,3)=[1,0;0,1];
mu_y_init=(max(X(:,1))+min(X(:,1)))/2;
mu_x1_init=max(X(:,2))/4+3*min(X(:,2))/4;
mu_x2_init=2*max(X(:,2))/4+2*min(X(:,2))/4;
mu_x3_init=3*max(X(:,2))/4+min(X(:,2))/4;
gamma=zeros(size(X,1),length(pi));%gamma(i,j),i是样本数量,j是聚类数量,gamma(i,j)表示样本i属于聚类j的概率，最大值表示为i的聚类
mu=[mu_x1_init,mu_y_init;mu_x2_init,mu_y_init;mu_x3_init,mu_y_init];
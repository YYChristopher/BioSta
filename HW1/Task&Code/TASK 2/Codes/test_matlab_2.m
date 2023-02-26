N=500;
pi_real=[3/10,5/10,2/10];%设定真实pi值
mu_real=[-5,4;-2,-5;4,4];%设定真实期望
cov_real(:,:,1)=[1,0;0,1];
cov_real(:,:,2)=[3,1;1,3];
cov_real(:,:,3)=[3,1;1,3];%设定真实方差
X_1=mvnrnd(mu_real(1,:),cov_real(:,:,1),N*pi_real(1));
X_2=mvnrnd(mu_real(2,:),cov_real(:,:,2),N*pi_real(2));
X_3=mvnrnd(mu_real(3,:),cov_real(:,:,3),N*pi_real(3));%生成真实数据点集
X=[X_1;X_2;X_3];
X=X(randperm(size(X,1)),:);
maxiter=100;
res_real=zeros(1,N);
res=zeros(1,N);
res_after=zeros(1,N);
for j=1:floor(N*pi_real(1))
    res_real(1,j)=1;
end
for j=floor(N*pi_real(1))+1:floor(N*(pi_real(1)+pi_real(2)))
    res_real(1,j)=2;
end
for j=floor(N*(pi_real(1)+pi_real(2)))+1:N
    res_real(1,j)=3;
end
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

NMI=zeros(maxiter,1);

%EM算法：
iter=50;
for i=1:iter
    for j=1:length(pi)
        gamma(:,j)=pi(j)*mvnpdf(X,mu(j,:),cov(:,:,j));
    end
    gamma=gamma./repmat(sum(gamma,2),1,size(gamma,2));%建立和gamma一样规模的矩阵，矩阵每一行为gamma(:,1)+gamma(:,2)值
    pi=sum(gamma,1)./size(gamma,1);
    mu=gamma'*X;
    mu=mu./repmat((sum(gamma,1))',1,size(mu,2));
    for j=1:length(pi)
        vari=repmat(gamma(:,j),1,size(X,2)).*(X-repmat(mu(j,:),size(X,1),1));
        cov(:,:,j)=(vari'*vari)/sum(gamma(:,j),1);
    end
    [c estimate]=max(gamma,[],2);

    nor_iter=pi(1)*mvnpdf(mesh,mu(1,:),cov(:,:,1))+pi(2)*mvnpdf(mesh,mu(2,:),cov(:,:,2))+pi(3)*mvnpdf(mesh,mu(3,:),cov(:,:,3));
    one_iter=find(estimate==1);
    two_iter=find(estimate==2);
    three_iter=find(estimate==3);
    for j=1:size(one_iter,1)
        res(1,j)=1;
    end
    for j=size(one_iter,1)+1:size(one_iter,1)+size(two_iter,1)
        res(1,j)=2;
    end
    for j=size(one_iter,1)+size(two_iter,1)+1:N
        res(1,j)=3;
    end
    NMI(i+1,1)=nmi_1(res,res_real);
end

%估计
[c estimate]=max(gamma,[],2);

nor=pi(1)*mvnpdf(mesh,mu(1,:),cov(:,:,1))+pi(2)*mvnpdf(mesh,mu(2,:),cov(:,:,2))+pi(3)*mvnpdf(mesh,mu(3,:),cov(:,:,3));
figure(5);
contour(x,y,reshape(nor,size(x,2),size(y,2)));

one=find(estimate==1);
two=find(estimate==2);
three=find(estimate==3);

for j=1:size(one,1)
    res_after(1,j)=1;
end
for j=size(one,1)+1:size(one,1)+size(two,1)
    res_after(1,j)=2;
end
for j=size(one,1)+size(two,1)+1:N
    res_after(1,j)=3;
end

figure(6);
plot(X(one,1),X(one,2),'rx',X(two,1),X(two,2),'bo',X(three,1),X(three,2),'g<');
title('经过EM算法后的点');
legend('高斯分布1','高斯分布2','高斯分布3');
hold on;
confiEllipse(X_1,0.95);
hold on;
confiEllipse(X_2,0.95);
hold on;
confiEllipse(X_3,0.95);
hold on;
load old_faithful;
figure(1);
plot(old_faithful(:,1),old_faithful(:,2),'k.');
title('需要分类的数据集(喷泉位置)');
N=size(old_faithful,1);
x=1.5:0.5:5.5;
y=40:10:100;
pi=[1/2 1/2];%设置初始的权重
cov(:,:,1)=[1,0;0,1];
cov(:,:,2)=[1,0;0,1];%设置初始方差
X=old_faithful;
mu_x1_init=max(X(:,1))/4+3*min(X(:,1))/4;
mu_x2_init=3*max(X(:,1))/4+min(X(:,1))/4;
mu_y_init=(max(X(:,2))+min(X(:,2)))/2;
gamma=zeros(size(X,1),length(pi));%gamma(i,j),i是样本数量,j是聚类数量,gamma(i,j)表示样本i属于聚类j的概率，最大值表示为i的聚类
mu=[mu_x1_init,mu_y_init;mu_x2_init,mu_y_init];
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
end

[c estimate]=max(gamma,[],2);

one=find(estimate==1);
two=find(estimate==2);
X_1=[X(one,1) X(one,2)];
X_2=[X(two,1) X(two,2)];

figure(2);
plot(X(one,1),X(one,2),'rx',X(two,1),X(two,2),'bo');
axis([1 5.5 40 100]);
title('经过EM算法后的喷泉位置分类');
hold on;
confiEllipse(X_1,0.95);
hold on;
confiEllipse(X_2,0.95);
hold on;




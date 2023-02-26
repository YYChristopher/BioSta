%EM算法：
iter=40;
for i=1:iter
    for j=1:length(pi)
        gamma(:,j)=pi(j)*mvnpdf(X,mu(j,:),cov(:,:,j));
    end
    gamma=gamma./repmat(sum(gamma,2),1,size(gamma,2));%建立和gamma一样规模的矩阵，矩阵每一行为gamma(:,1)+gamma(:,2)值
    pi=sum(gamma,1)./size(gamma,1);%迭代pi值
    mu=gamma'*X;
    mu=mu./repmat((sum(gamma,1))',1,size(mu,2));%迭代期望
    for j=1:length(pi)
        vari=repmat(gamma(:,j),1,size(X,2)).*(X-repmat(mu(j,:),size(X,1),1));
        cov(:,:,j)=(vari'*vari)/sum(gamma(:,j),1);%迭代方差
    end
end

%估计
[c estimate]=max(gamma,[],2);

nor=pi(1)*mvnpdf(mesh,mu(1,:),cov(:,:,1))+pi(2)*mvnpdf(mesh,mu(2,:),cov(:,:,2))+pi(3)*mvnpdf(mesh,mu(3,:),cov(:,:,3));%生成迭代后的混合高斯分布
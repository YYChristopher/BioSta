%EM�㷨��
iter=40;
for i=1:iter
    for j=1:length(pi)
        gamma(:,j)=pi(j)*mvnpdf(X,mu(j,:),cov(:,:,j));
    end
    gamma=gamma./repmat(sum(gamma,2),1,size(gamma,2));%������gammaһ����ģ�ľ��󣬾���ÿһ��Ϊgamma(:,1)+gamma(:,2)ֵ
    pi=sum(gamma,1)./size(gamma,1);%����piֵ
    mu=gamma'*X;
    mu=mu./repmat((sum(gamma,1))',1,size(mu,2));%��������
    for j=1:length(pi)
        vari=repmat(gamma(:,j),1,size(X,2)).*(X-repmat(mu(j,:),size(X,1),1));
        cov(:,:,j)=(vari'*vari)/sum(gamma(:,j),1);%��������
    end
end

%����
[c estimate]=max(gamma,[],2);

nor=pi(1)*mvnpdf(mesh,mu(1,:),cov(:,:,1))+pi(2)*mvnpdf(mesh,mu(2,:),cov(:,:,2))+pi(3)*mvnpdf(mesh,mu(3,:),cov(:,:,3));%���ɵ�����Ļ�ϸ�˹�ֲ�
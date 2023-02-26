
function confiEllipse(datamatrix,p)
%print 2-demension confidence ellipse
%In:n��2 matrix, confidence probability p
 
data = datamatrix;
covariance = cov(data);
[eigenvec,eigenval] = eig(covariance);
 
[sortEigenval,index] = sort(diag(eigenval),'descend');
sortEigenvec = eigenvec(:,index);
 
largestEigenval = sortEigenval(1);
smallestEigenval = sortEigenval(end); %����С����ֵ
largestEigenvec = sortEigenvec(:,1); %�����������������Բ�ĳ��ᣩ
 
angle = atan2(largestEigenvec(2), largestEigenvec(1)); %����x��������������֮��ĽǶ�, [-pi,pi]
 
if(angle < 0) 
    angle = angle + 2*pi;
end
 
avg = mean(data); %�����������ݵľ�ֵ
 
%����������Բ�Ĳ���
chisquareVal = sqrt(chi2inv(p,2)); %����ֵ
thetaGrid = linspace(0,2*pi); 
phi = angle; %��ת�Ƕ�
X0=avg(1);
Y0=avg(2); 
a=chisquareVal*sqrt(largestEigenval); %�־� ����
b=chisquareVal*sqrt(smallestEigenval);
 
ellipseXR = a*cos( thetaGrid ); %������ֱ������ϵ
ellipseYR = b*sin( thetaGrid );
 
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ]; %��ת����
 
rEllipse = [ellipseXR;ellipseYR]' * R; %��ת
plot(rEllipse(:,1) + X0,rEllipse(:,2) + Y0,'k-','linewidth',1.5);
axis square
end


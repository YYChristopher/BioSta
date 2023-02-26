
function confiEllipse(datamatrix,p)
%print 2-demension confidence ellipse
%In:n×2 matrix, confidence probability p
 
data = datamatrix;
covariance = cov(data);
[eigenvec,eigenval] = eig(covariance);
 
[sortEigenval,index] = sort(diag(eigenval),'descend');
sortEigenvec = eigenvec(:,index);
 
largestEigenval = sortEigenval(1);
smallestEigenval = sortEigenval(end); %求最小特征值
largestEigenvec = sortEigenvec(:,1); %求最大特征向量（椭圆的长轴）
 
angle = atan2(largestEigenvec(2), largestEigenvec(1)); %计算x轴和最大特征向量之间的角度, [-pi,pi]
 
if(angle < 0) 
    angle = angle + 2*pi;
end
 
avg = mean(data); %计算两列数据的均值
 
%计算置信椭圆的参数
chisquareVal = sqrt(chi2inv(p,2)); %卡方值
thetaGrid = linspace(0,2*pi); 
phi = angle; %旋转角度
X0=avg(1);
Y0=avg(2); 
a=chisquareVal*sqrt(largestEigenval); %轮距 长度
b=chisquareVal*sqrt(smallestEigenval);
 
ellipseXR = a*cos( thetaGrid ); %作用于直角坐标系
ellipseYR = b*sin( thetaGrid );
 
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ]; %旋转矩阵
 
rEllipse = [ellipseXR;ellipseYR]' * R; %旋转
plot(rEllipse(:,1) + X0,rEllipse(:,2) + Y0,'k-','linewidth',1.5);
axis square
end


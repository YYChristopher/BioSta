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
end%µü´úEMËã·¨
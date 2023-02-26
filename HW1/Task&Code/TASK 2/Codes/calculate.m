load old_faithful;
mu_x=mean(old_faithful(:,1));
mu_y=mean(old_faithful(:,2));
vari=cov(old_faithful);
mu_x
mu_y
vari
figure(1);
plot(old_faithful(:,1),old_faithful(:,2),'k.');

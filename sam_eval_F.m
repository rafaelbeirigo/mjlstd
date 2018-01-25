sam_constants;
sam_parameters;

R=100;
T=100;

thetas=zeros([T,R]);
xs=zeros([size(x_0),T,R]);

Fs=F_subopt;

test_stability(lambda,Fs,As,Bs,P,N)

for r=1:R
  [thetas(:,r) xs(:,:,:,r)]=sam_run_episode(T,theta_0,x_0,P,As,Bs,Fs,r);
end

xs_avg=reshape(mean(xs(2,:,:,:),4),[T,1]);
xs_std=reshape(std(xs(2,:,:,:),0,4),[T,1]);

csvwrite('xs.csv',[[1:size(xs,3)]',xs_avg,xs_std]);

function [thetas xs]=sam_run_episode(T,theta_0,x_0,P,As,Bs,Fs,seed)
  rand('seed',seed);

  thetas=zeros([T,1]);
  xs=zeros([size(x_0),T]);

  thetas(1)=theta_0;
  xs(:,:,1)=x_0;
  for t=2:T
    thetas(t)=find(cumsum(P(thetas(t-1),:))>rand(),1);

    A=As(:,:,thetas(t-1));
    B=Bs(:,:,thetas(t-1));
    F=Fs(:,:,thetas(t-1));
    xs(:,:,t)=(A+B*F)*xs(:,:,t-1);
  end
end

close all;
clear all;

constants;
parameters;

rand('seed',rand_seed);

Theta=zeros(N,K);
for i=1:N
  Theta(i,1)=i;
end

Pc=cumsum(P,2);

for i=1:N
  Upsilons(:,:,i)=eye(size(As,1));
  Sum(:,:,i)=zeros(size(As,1));
  Ys(:,:,i)=zeros(size(As,1));
end

Fs=F_opt;

for i=1:N
  for k=1:K-1
    Theta(i,k+1)=find(Pc(Theta(i,k),:)>=rand())(1);

    U=Upsilons(:,:,i);
    A=As(:,:,Theta(i,k+1));
    B=Bs(:,:,Theta(i,k+1));
    C=Cs(:,:,Theta(i,k+1));
    D=Ds(:,:,Theta(i,k+1));
    F=Fs(:,:,Theta(i,k+1));
    Y1=Ys(:,:,Theta(i,k));
    Y2=Ys(:,:,Theta(i,k+1));

    Bcal(:,:,i)=U'*(C'*C+F'*D'*D*F)*U;
    Ccal(:,:,i)=U'*((A+B*F)'*Y2*(A+B*F)-Y1)*U;
    Dcal(:,:,i)=Bcal(:,:,i)+Ccal(:,:,i);

    Sum_old=Sum(:,:,i);
    Sum(:,:,i)=Sum(:,:,i)+power(lambda,k)*Dcal(:,:,i);

    maxDiff=max(max(abs(Sum(:,:,i)-Sum_old)));
    maxDiff
    if maxDiff<epsilon
      break
    end

    Upsilons(:,:,i)=(A+B*F)*Upsilons(:,:,i);
  end
end

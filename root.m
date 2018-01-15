close all;
clear all;

constants;
parameters;

rand('seed',rand_seed);

Theta=zeros(N,T);
for i=1:N
  Theta(i,1)=i;
end

Pc=cumsum(P,2);

for i=1:N
  Upsilon(:,:,i)=eye(size(A,1));
end

for k=1:T-1
  for i=1:N
    Theta(i,k+1)=find(Pc(Theta(i,k),:)>=rand())(1);
  end

  for i=1:N
    Upsilon(:,:,i)=A(:,:,Theta(i,k+1))+B(:,:,Theta(i,k+1))*F_opt(:,:,Theta(i,k+1))*Upsilon(:,:,i)
  end
end

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
  Upsilon(:,:,i)=eye(N);
end

for k=1:T-1
  for i=1:N
    Theta(i,k+1)=find(Pc(Theta(i,k),:)>=rand())(1);
  end
end

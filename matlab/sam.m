close all;
clear all;

graphics_toolkit('gnuplot');

sam_constants;
sam_parameters;

Fe=zeros(N,J,R);
Fee=zeros(J*T,[3,1],R);
Ys=zeros(2,2,3,10,25);
Ys=zeros(2,2,3,10,25);
Fs=F_opt;

test_stability(lambda,Fs,As,Bs,P,N);

for r=1:R
  [Fss{end+1} Yss{end+1} Fe(:,:,r) Fee(:,:,r)]=mjlstd(lambda,J,T,K,epsilon,N,P,As,Bs,Cs,Ds,F_opt,Fs,X_opt,r,c,eta);
end

%% Plot Y error history
num_lin=size(Yss{1},2);
num_col=size(Yss{1}{1},1)*size(Yss{1}{1},2)*size(Yss{1}{1},3);
Yss_=zeros([num_lin,num_col]);
for i=1:size(Yss{1},2)
  Yss_(i,:,:,:,:,:,:)=reshape(Yss{1}{i}-X_opt,[1,num_col]);
end
figure(1)
plot(Yss_)

%% Plot gain error history
num_lin=size(Fss{1},2);
num_col=size(Fss{1}{1},1)*size(Fss{1}{1},2)*size(Fss{1}{1},3);
Fss_=zeros([num_lin,num_col]);
for i=1:size(Fss{1},2)
  Fss_(i,:,:,:,:,:,:)=reshape(Fss{1}{i}-F_opt,[1,num_col]);
end
figure(2)
plot(Fss_)

Fe_avg=mean(Fe,3);
Fe_std=std(Fe,0,3);
Fee_avg=mean(Fee,3);
Fee_std=std(Fee,0,3);

csvwrite('Fe_avg_std.csv',[[1:size(Fe_avg,1)]',Fe_avg,Fe_std]);
csvwrite('Fee_avg_std.csv',[[1:size(Fee_avg,1)]',Fee_avg,Fee_std]);

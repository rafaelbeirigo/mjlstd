close all;
clear all;

sam_constants;
sam_parameters;

Fe=zeros(N,J,R);
Fee=zeros(J*T,[3,1],R);
Ys=zeros(2,2,3,10,25);
Ys=zeros(2,2,3,10,25);
Fs=F_opt;

maxVal=0;
for i=1:N
  A=As(:,:,i);
  B=Bs(:,:,i);
  F=Fs(:,:,i);
  Gamma(:,:,i)=A+B*F;

  curVal=power(norm(kron(Gamma(:,:,i),Gamma(:,:,i))),-1);
  if curVal>maxVal
    maxVal=curVal;
  end
end
if !(power(lambda,2)<maxVal)
  disp('WARNING: Lemma 2''s sufficient condition for lambda is not satisfied');
end

for r=1:R
  [Fe(:,:,r) Fee(:,:,r) Yse(:,:,:,:,r)]=mjlstd(lambda,J,T,K,epsilon,N,P,As,Bs,Cs,Ds,F_opt,Fs,X_opt,r);
end

Fe_avg=mean(Fe,3);
Fe_std=std(Fe,0,3);
Fee_avg=mean(Fee,3);
Fee_std=std(Fee,0,3);
Yse_avg=mean(Yse,5);
Yse_std=std(Yse,0,5);

csvwrite('Fe_avg_std.csv',[[1:size(Fe_avg,1)]',Fe_avg,Fe_std]);
csvwrite('Fee_avg_std.csv',[[1:size(Fee_avg,1)]',Fee_avg,Fee_std]);
%% csvwrite('Yse_avg_std.csv',[[1:size(Yse_avg,1)]',Yse_avg,Yse_std]);

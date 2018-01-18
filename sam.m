close all;
clear all;

sam_constants;
sam_parameters;

Fe=zeros(N,J,R);
Fee=zeros(J*T,[3,1],R);

for r=1:R
  [Fe(:,:,r) Fee(:,:,r)] = mjlstd(lambda,J,T,K,epsilon,N,P,As,Bs,Cs,Ds,F_opt,X_opt,r);
end

Fe_avg = mean(Fe,3);
Fe_std = std(Fe,0,3);
Fee_avg = mean(Fee,3);
Fee_std = std(Fee,0,3);

csvwrite('Fe_avg_std.csv',[[1:size(Fe_avg,1)]',Fe_avg,Fe_std]);
csvwrite('Fee_avg_std.csv',[[1:size(Fee_avg,1)]',Fee_avg,Fee_std]);

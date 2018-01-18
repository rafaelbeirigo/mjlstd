function [Fe Fee] = mjlstd(lambda,J,T,K,epsilon,N,P,As,Bs,Cs,Ds,F_opt,X_opt)
  Theta=zeros(N,K);
  for i=1:N
    Theta(i,1)=i;
  end
  Pc=cumsum(P,2);
  Fe=zeros(N,J);
  Fee=zeros(J*T,[3,1]);

  rand('seed',seed);
  for i=1:N
    Ys(:,:,i)=zeros(size(As,1));
  end
  Yconverged=zeros(N);
  Fs=F_opt*0.75;
  for j=1:J
    for t=1:T
      for l=1:N
        A=As(:,:,l);
        B=Bs(:,:,l);
        D=Ds(:,:,l);
        S=Ys(:,:,l);

        FsAux(:,:,l)=-inv(B'*S*B+D'*D)*B'*S*A;
        Fee((j-1)*T+t,l,r)=Fee((j-1)*T+t,l,r)+max(abs(F_opt(:,:,l)-FsAux(:,:,l)));
      end

      for i=1:N
        if Yconverged(i)
          continue
        end
        Upsilons(:,:,i)=eye(size(As,1));
        Sum(:,:,i)=zeros(size(As,1));
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

          SumAux=Sum(:,:,i);
          Sum(:,:,i)=Sum(:,:,i)+power(lambda,k)*Dcal(:,:,i);

          maxDiff=max(max(abs(Sum(:,:,i)-SumAux)));
          if maxDiff<epsilon
            break
          end

          Upsilons(:,:,i)=(A+B*F)*Upsilons(:,:,i);
        end
      end

      gamma=1/t;
      for i=1:N
        if Yconverged(i)
          continue
        end
        Yaux=Ys(:,:,i);
        Ys(:,:,i)=Ys(:,:,i)+gamma*Sum(:,:,i);

        maxDiff=max(max(abs(Ys(:,:,i)-Yaux)));
        if maxDiff<epsilon
          Yconverged(i)=1;
        end
      end
    end

    for l=1:N
      A=As(:,:,l);
      B=Bs(:,:,l);
      D=Ds(:,:,l);
      S=Ys(:,:,l);

      Fs(:,:,l)=-inv(B'*S*B+D'*D)*B'*S*A;
      Fe(l,j,r)=max(max(abs(F_opt(:,:,l)-Fs(:,:,l))));
    end

    if j>1
      maxDiff=abs(Fe(:,j,r)-Fe(:,j-1,r));
      if maxDiff<epsilon
        break
      end
    end
  end

  Fe_avg = mean(Fe,3);
  Fe_std = std(Fe,0,3);
  Fee_avg = mean(Fee,3);
  Fee_std = std(Fee,0,3);

  csvwrite('Fe_avg_std.csv',[[1:size(Fe_avg,1)]',Fe_avg,Fe_std]);
  csvwrite('Fee_avg_std.csv',[[1:size(Fee_avg,1)]',Fee_avg,Fee_std]);
end

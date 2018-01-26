function [Fss Fe Fee Yse] = mjlstd(lambda,J,T,K,epsilon,N,P,As,Bs,Cs,Ds,F_opt,Fs,X_opt,seed,c,eta)
  Theta=zeros(N,K);
  for i=1:N
    Theta(i,1)=i;
  end
  Pc=cumsum(P,2);
  Fe=zeros(N,J);
  Fee=zeros(J*T,[3,1]);
  for i=1:N
    Ys(:,:,i)=zeros(size(As,1));
  end
  Yconverged=zeros(N);
  rand('seed',seed);
  for j=1:J
    for t=1:T
      for i=1:N
        A=As(:,:,i);
        B=Bs(:,:,i);
        D=Ds(:,:,i);
        S=Ys(:,:,i);

        FsAux(:,:,i)=-inv(B'*S*B+D'*D)*B'*S*A;
        Fee((j-1)*T+t,i)=Fee((j-1)*T+t,i)+max(abs(F_opt(:,:,i)-FsAux(:,:,i)));
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

      gamma=c*power(t,-eta);

      for i=1:N
        if Yconverged(i)
          continue
        end
        Yaux(:,:,i)=Ys(:,:,i);
        Ys(:,:,i)=Ys(:,:,i)+gamma*Sum(:,:,i);
        diff=abs(Ys(:,:,i)-Yaux(:,:,i));
        maxDiff=max(max(diff));
        if maxDiff<epsilon
          Yconverged(i)=1;
        end
      end
      Yse(:,:,:,(j-1)*T+t)=abs(Ys-Yaux);
    end

    for i=1:N
      A=As(:,:,i);
      B=Bs(:,:,i);
      D=Ds(:,:,i);
      S=Ys(:,:,i);

      Fs(:,:,i)=-inv(B'*S*B+D'*D)*B'*S*A;
      Fe(i,j)=max(max(abs(F_opt(:,:,i)-Fs(:,:,i))));
    end

    if j>1
      maxDiff=abs(Fe(:,j)-Fe(:,j-1));
      if maxDiff<epsilon
        break
      end
    end
  end
end

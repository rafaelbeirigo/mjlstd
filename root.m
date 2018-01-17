close all;
clear all;

constants;
parameters;

Theta=zeros(N,K);
for i=1:N
  Theta(i,1)=i;
end
Pc=cumsum(P,2);
Fe=zeros(N,J);

Fee=zeros(J*T,[3,1]);

for r=1:R
  rand('seed',r);
  for i=1:N
    Ys(:,:,i)=zeros(size(As,1));
  end
  Fs=F_opt/2;
  for j=1:J
    for i=1:N
      for t=1:T
        for l=1:N
          A=As(:,:,l);
          B=Bs(:,:,l);
          D=Ds(:,:,l);
          S=Ys(:,:,l);

          FsAux(:,:,l)=-inv(B'*S*B+D'*D)*B'*S*A;
          Fee((j-1)*T+t,l)=Fee((j-1)*T+t,l)+max(abs(F_opt(:,:,l)-FsAux(:,:,l)));
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

        gamma=1/t;
        Yaux=Ys(:,:,i);
        Ys(:,:,i)=Ys(:,:,i)+gamma*Sum(:,:,i);

        maxDiff=max(max(abs(Ys(:,:,i)-Yaux)));
        if maxDiff<epsilon
          break
        end
      end
    end

    for l=1:N
      A=As(:,:,l);
      B=Bs(:,:,l);
      D=Ds(:,:,l);
      S=Ys(:,:,l);

      Fs(:,:,l)=-inv(B'*S*B+D'*D)*B'*S*A;
      Fe(l,j)=max(max(abs(F_opt(:,:,l)-Fs(:,:,l))));
    end

    if j>1
      maxDiff=abs(Fe(:,j)-Fe(:,j-1));
      if maxDiff<epsilon
        break
      end
    end
  end
end
csvwrite('Fe.csv',[[1:J]',Fe'/R]);
csvwrite('Fee.csv',[[1:size(Fee,1)]',Fee/R]);

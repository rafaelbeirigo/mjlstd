function [r] = stabilizes(Fs,As,Bs,P,N)
  %% Tests if F stabilizes (A,B) according to Lemma 3 of (1)
  %%
  %% (1) O. L. V. Costa and J. C. C. Aya, "Monte Carlo
  %% TD(\lambda)-methods for the optimal control of
  %% discrete-time Markovian jump linear systems,"
  %% Automatica, vol. 38, pp. 217–225, 2002.

  for j=1:N
    I(:,:,j)=zeros([size(As(:,:,1))]);
    for i=1:N
      F=Fs(:,:,i);
      A=As(:,:,i);
      B=Bs(:,:,i);
      G=A+B*F;

      I(:,:,j)=I(:,:,j)+P(i,j)*G*G';
    end
    r(j)=max(abs(eig(I(:,:,j))));
  end
  r
  r_sigma=max(r);
  if (r_sigma<1)
    disp "F stabilizes (A,B)"
  else
    disp "WARNING: F does not stabilize (A,B)"
  end
end

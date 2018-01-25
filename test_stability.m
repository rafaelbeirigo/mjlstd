function [r] = test_stability(Fs,As,Bs,P,N)
  %% Tests stability of (A,B) according to
  %%
  %% (1) O. L. V. Costa and J. C. C. Aya, "Monte Carlo
  %% TD(\lambda)-methods for the optimal control of
  %% discrete-time Markovian jump linear systems,"
  %% Automatica, vol. 38, pp. 217–225, 2002.

  I=[];
  for i=1:N
    F=Fs(:,:,i);
    A=As(:,:,i);
    B=Bs(:,:,i);
    G=A+B*F;

    I=blkdiag(I,kron(G,G))
  end
  r_sigma=max(abs(eig(I)));
  if (r_sigma<1)
    disp "F stabilizes (A,B)"
  else
    disp "WARNING: F does not stabilize (A,B)"
  end
end

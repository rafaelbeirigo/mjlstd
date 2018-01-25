function [r] = test_stability(lambda,Fs,As,Bs,P,N)
  %% Tests stability of (A,B) according to
  %%
  %% (1) O. L. V. Costa and J. C. C. Aya, "Monte Carlo
  %% TD(\lambda)-methods for the optimal control of
  %% discrete-time Markovian jump linear systems,"
  %% Automatica, vol. 38, pp. 217–225, 2002.

  I=[];
  v_max=-Inf;
  for i=1:N
    A=As(:,:,i);
    B=Bs(:,:,i);
    F=Fs(:,:,i);
    G=A+B*F;

    kr=kron(G,G);
    %% Fs test
    I=blkdiag(I,kr);

    %% lambda test
    v=power(norm(kr),-1);
    if v>v_max
      v_max=v;
    end
  end

  if !(max(abs(eig(I)))<1)
    disp('WARNING: F does not satisfy Lemma 3 (stabilizability)');
  end
  if !(power(lambda,2)<v_max)
    disp('WARNING: lambda does not satisfy Lemma 2 (convergence)');
  end
end

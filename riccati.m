function [F] = riccati(T, N, A, B, C, D, R)
  nx = size(A,1);

  for i=1:N
    X(:,:,i) = zeros(nx);
  end

  for t=1:T
    %% Estimate X (EX)
    for i=1:N
      EX(:,:,i) = zeros(nx);
      for j=1:N
        EX(:,:,i) = EX(:,:,i) + R(i,j)*X(:,:,j);
      end
    end

    %% Calculate F and X based on the estimate of X
    for i=1:N
      F(:,:,i) = -(D(:,:,i)'*D(:,:,i) + B(:,:,i)'*EX(:,:,i)*B(:,:,i))\B(:,:,i)'*EX(:,:,i)*A(:,:,i);
      X(:,:,i) = A(:,:,i)'*EX(:,:,i)*A(:,:,i) + A(:,:,i)'*EX(:,:,i)*B(:,:,i)*F(:,:,i) + C(:,:,i)'*C(:,:,i);
    end
  end
end

close all;
clear all;

constants;
parameters;

Theta = zeros(N, T);
for i=1:N
  Theta(i, 1) = i;
end

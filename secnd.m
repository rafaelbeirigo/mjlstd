close all
clear all

lmb = 0.1;    %lambda
epY = 1e-4;  %tolerance for convergence of Y
epF = 1e-4;  %tolerance for convergence of F
c = 1;       %constant that multiplies gamma
elmax = 5;  %maximum number of policy improvements
kmax = 50;   %lookahead horizon
tmax = 1000; %maximum iterates allowed to Y(t)

N = 3;
P = [0.67 0.17 0.16;
     0.30 0.47 0.23;
     0.26 0.10 0.64];

A(:,:,1) = [0 1; -2.5 3.2];
A(:,:,2) = [0 1; -4.3 4.5];
A(:,:,3) = [0 1; 5.3 -5.2];

nx = 2;
nu = 1;

for i=1:N
   B(:,:,i) = [0; 1];
   F(:,:,i) = -A(2,:,i); %initial, stabilizing controller
end

%a suboptimal controller
%F(:,:,1) = [2.4854  -2.2734]; F(:,:,2) = [4.3655  -4.4401]; F(:,:,3) = [-5.2799   6.0555];

%actually optimal controller
%F(:,:,1) = [2.6819 -2.6721]; F(:,:,2) = [4.4230 -3.9265]; F(:,:,3) = [-5.4260 6.0883]; 

Q(:,:,1) = [3.6 -3.8; -3.8 4.87];
Q(:,:,2) = [10 -3; -3 8];
Q(:,:,3) = [5 -4.5; -4.5 4.5];

R(:,:,1) = 2.6;
R(:,:,2) = 1.165;
R(:,:,3) = 1.111;

%%%%%%%%%%%%
% MSS test %
%%%%%%%%%%%%
%Acal = [];
%for i=1:N
%   Acal = blkdiag(Acal, kron(A(:,:,i)+B(:,:,i)*F(:,:,i),A(:,:,i)+B(:,:,i)*F(:,:,i)));
%end
%Acal = kron(P',eye(nx^2))*Acal;
%disp(max(abs(eig(Acal))));


Y = zeros(nx,nx,N);
Yh = [];
Fh = [];

el = 0;
erF = 2*epF;
while (el<elmax && erF>epF)
   el = el+1;
   newY = 0*Y;
   t = 0;
   erY = 2*epY;
   while (t<=tmax && erY>epY)
      t = t+1;
      gam = c/t;
      erY = 0;
      for i=1:N
         AA = A(:,:,i) + B(:,:,i)*F(:,:,i);
         QQ = Q(:,:,i) + F(:,:,i)'*R(:,:,i)*F(:,:,i);
         th = hmmgenerate(kmax+2, [[zeros(1,i),1,zeros(1,N-i)];[zeros(N,1),P]],eye(N+1))-1;
         Dsum = 0;
         Up = eye(nx);
         for k=0:kmax
            Dk = lmb^k*Up'*(QQ + AA'*Y(:,:,th(k+2))*AA - Y(:,:,th(k+1)))*Up;
            Dsum = Dsum + Dk;
            Up = AA*Up;
         end
         newY(:,:,i) = Y(:,:,i) + gam*Dsum;
         erY = max(erY,norm(newY(:,:,i) - Y(:,:,i)));
      end
      Y = newY;
      Yh(:,:,:,end+1) = Y;
   end

   erF = 0;
   for i=1:N
      aux = -B(:,:,i)'*Y(:,:,i)*B(:,:,i) + R(:,:,i);
      newF(:,:,i) = aux\B(:,:,i)'*Y(:,:,i)*A(:,:,i);
      erF = max(erF,norm(newF(:,:,i) - F(:,:,i)));
   end
   F = newF;
   Fh(:,:,:,end+1) = F;
end

fprintf('erF = %f\n',erF);

%Riccati verification
erY = 2*epY;
X = Y;
Xnew = Y;
while (erY>epY)
   erY = 0;
   for i=1:N
      E = 0;
      for j=1:N
         E = E + P(i,j)*X(:,:,j);
      end
      Yric(:,:,i) = E;
      AA = A(:,:,i) + B(:,:,i)*F(:,:,i);
      QQ = Q(:,:,i) + F(:,:,i)'*R(:,:,i)*F(:,:,i);
      newX(:,:,i) = AA'*E*AA + QQ;
      erY = max(erY,norm(newX(:,:,i) - X(:,:,i)));
   end
   X = newX;
end
for i=1:N
   aux = -B(:,:,i)'*Yric(:,:,i)*B(:,:,i) + R(:,:,i);
   FF(:,:,i) = aux\B(:,:,i)'*Yric(:,:,i)*A(:,:,i);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
clf
suptitle('Entries of F at each el-step (blue) vs true optimal gain (red)');

subplot(321)
rw = 1; cl = 1; i = 1;
plot(1:size(Fh,4), [reshape(Fh(rw,cl,i,:),1,size(Fh,4));...
      FF(rw,cl,i,:)*ones(1,size(Fh,4))])
set(gca,'XLim',[1 size(Fh,4)])
ylabel(sprintf('F_%d(%d,%d)',i,rw,cl));
xlabel('el-step');

subplot(322)
rw = 1; cl = 2; i = 1;
plot(1:size(Fh,4), [reshape(Fh(rw,cl,i,:),1,size(Fh,4));...
      FF(rw,cl,i,:)*ones(1,size(Fh,4))])
set(gca,'XLim',[1 size(Fh,4)])
ylabel(sprintf('F_%d(%d,%d)',i,rw,cl));
xlabel('el-step');

subplot(323)
rw = 1; cl = 1; i = 2;
plot(1:size(Fh,4), [reshape(Fh(rw,cl,i,:),1,size(Fh,4));...
      FF(rw,cl,i,:)*ones(1,size(Fh,4))])
set(gca,'XLim',[1 size(Fh,4)])
ylabel(sprintf('F_%d(%d,%d)',i,rw,cl));
xlabel('el-step');

subplot(324)
rw = 1; cl = 2; i = 2;
plot(1:size(Fh,4), [reshape(Fh(rw,cl,i,:),1,size(Fh,4));...
      FF(rw,cl,i,:)*ones(1,size(Fh,4))])
set(gca,'XLim',[1 size(Fh,4)])
ylabel(sprintf('F_%d(%d,%d)',i,rw,cl));
xlabel('el-step');

subplot(325)
rw = 1; cl = 1; i = 3;
plot(1:size(Fh,4), [reshape(Fh(rw,cl,i,:),1,size(Fh,4));...
      FF(rw,cl,i,:)*ones(1,size(Fh,4))])
set(gca,'XLim',[1 size(Fh,4)])
ylabel(sprintf('F_%d(%d,%d)',i,rw,cl));
xlabel('el-step');

subplot(326)
rw = 1; cl = 2; i = 3;
plot(1:size(Fh,4), [reshape(Fh(rw,cl,i,:),1,size(Fh,4));...
      FF(rw,cl,i,:)*ones(1,size(Fh,4))])
set(gca,'XLim',[1 size(Fh,4)])
ylabel(sprintf('F_%d(%d,%d)',i,rw,cl));
xlabel('el-step');

%----------------------------------------------------------%
figure(2)
clf
suptitle('Entries of Y at each t,el-step (blue) vs true Riccati solution (red)');

pltnum = 330;
for i=1:N
   for rw=1:nx
      for cl=rw:nx
         pltnum = pltnum + 1;
         subplot(pltnum);
         plot(1:size(Yh,4), [reshape(Yh(rw,cl,i,:),1,size(Yh,4));Yric(rw,cl,i)*ones(1,size(Yh,4))])
         set(gca,'XLim',[1 size(Yh,4)]);
         ylabel(sprintf('Y_%d(%d,%d)',i,rw,cl));
         xlabel('(t,el)-step');
      end
   end
end

%%print 'fig2' -dpdf
%%!ps2pdf -dEPSCrop fig2.eps fig2.pdf


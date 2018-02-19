clear all

lmb = 0.1;   %lambda
epY = 1e-2;  %tolerance for convergence of Y
epF = 1e-3;  %tolerance for convergence of F
c = 1;       %constant that multiplies gamma
elmax = 20;   %maximum number of policy improvements
kmax = 10;   %lookahead horizon
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
Fop(:,:,1) = [2.6819 -2.6721]; Fop(:,:,2) = [4.4230 -3.9265]; Fop(:,:,3) = [-5.4260 6.0883]; 

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
Yh = Y;
Fh = F;
erFh = [];

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
         th = hmmgenerate(kmax+2, [[zeros(1,i),1,zeros(1,N-i)];[zeros(N,1),P]],eye(N+1))-1;
         Dsum = 0;
         Up = eye(nx);
         for k=0:kmax
            AA = A(:,:,th(k+2)) + B(:,:,th(k+2))*F(:,:,th(k+2));
            QQ = Q(:,:,th(k+2)) + F(:,:,th(k+2))'*R(:,:,th(k+2))*F(:,:,th(k+2));
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
   if (isempty(erFh))
      erFh = [erF erF];
   else
      erFh(end+1) = erF;
   end
end


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
ftsiz = 6;

figure(1)
clf

pltnum = 320;
for i=1:N
   for rw=1:1
      for cl=1:nx
         pltnum = pltnum + 1;
         subplot(pltnum);
         set(gca,'FontSize',ftsiz)
         hold on
         stairs(0:size(Fh,4)-1, reshape(Fh(rw,cl,i,:),1,size(Fh,4)),'b')
         plot(0:size(Fh,4)-1, FF(rw,cl,i,:)*ones(1,size(Fh,4)),'r--')
         hold off
         set(gca,'XLim',[0 size(Fh,4)-1])
         ylabel(sprintf('F_%d(%d,%d)',i,rw,cl));
         xlabel('el-step');
      end
   end
end
suptitle('Entries of F at each el-step (blue) vs true optimal gain (red)');

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
         set(gca,'FontSize',ftsiz)
         hold on
         plot(1:size(Yh,4), reshape(Yh(rw,cl,i,:),1,size(Yh,4)),'b')
         plot(1:size(Yh,4), Yric(rw,cl,i)*ones(1,size(Yh,4)),'r--')
         hold off
         set(gca,'XLim',[1 size(Yh,4)]);
         ylabel(sprintf('Y_%d(%d,%d)',i,rw,cl));
         xlabel('(t,el)-step');
      end
   end
end

%%print 'fig2' -dpdf
%%!ps2pdf -dEPSCrop fig2.eps fig2.pdf

%----------------------------------------------------------%
figure(3)
clf
suptitle('Entries of \Delta (see the 2002 Automatica paper, at the end of Sec. 5)')

erFop = zeros(1,size(Fh,4));
Delt = [];
for t=1:size(Fh,4)
   aux = 0;
   for i=1:N
      aux = max(aux, norm(Fh(:,:,i,t)-Fop(:,:,i)));
   end
   erFop(t) = aux;

   Delt(:,:,:,end+1) = zeros(size(F));
   for i=1:N
      Delt(:,:,i,end) = 100*abs( (Fh(:,:,i,t)-Fop(:,:,i))./Fop(:,:,i) );
   end
end

pltnum = 320;
for i=1:N
   for rw=1:1
      for cl=1:nx
         pltnum = pltnum + 1;
         subplot(pltnum);
         set(gca,'FontSize',ftsiz)
         hold on
         stairs(0:size(Delt,4)-1, reshape(Delt(rw,cl,i,:),1,size(Delt,4)),'b')
         plot(0:size(Delt,4)-1, ones(1,size(Delt,4)),'k:')
         hold off
         set(gca,'XLim',[0 size(Delt,4)-1])
         ylabel(sprintf('\\Delta_%d(%d,%d)',i,rw,cl));
         xlabel('el-step');
         if (isempty(intersect(get(gca,'YTick'),1)))
            set(gca,'YTick',sort([1 get(gca,('YTick'))]))
         end
      end
   end
end

for i=1:gcf
   figure(i)
   aux = sprintf('fig%d',i);
   print('-depsc2',aux)
   eval(sprintf('!ps2pdf -dEPSCrop %s.eps %s.pdf\n',aux,aux))
end


close all
clear 
clc

global Tm TM  pC w g kn pI kR RunT  ny   sigmaa sigmag  k_eta C
%%
disp('Initializing...') 

%%%%% Intermittency %%%%%% 
Tm = 0.1;
TM = 0.2;  

% landmarks
ny = 0.1;
sigmaa = 0.01;
sigmag = 0.01;
w = 1;
g = 9.81*[0,0,1]';

n = 10; 
pI  = [2*randn(2,n);zeros(1,n)];
% pI = [ 0.6739   -0.5004    1.6567    1.6212   -0.3342    3.0828    0.1430   -1.6905;
%    -0.0837   -0.8073    0.4521    0.0959   -1.0428   -0.1323    2.3234   -0.9586;
%          0         0         0         0         0         0         0         0];

% kn   = 1./(diag((pI-mean(pI,2))'*(pI-mean(pI,2))));
% kn   = 1./(diag(pI'*pI)); 
% kn   = diag(kn); 
kn   = eye(n)/n;
pC   = sum(pI*kn,2)./sum(diag(kn)); 
M    = (pI-pC)*kn*(pI-pC)';
Mbar = (trace(M)*eye(3)-M)/2;
Mbar_M = max(eig(Mbar));
Mbar_m = min(eig(Mbar));

Mline = trace(Mbar^2)*eye(3) - 2*Mbar^2;
MLine_M = max(eig(Mline));

%%%%%% zero-hold
% mu2 = norm(Mbar,'fro');%max(MLine_M/(4*Mbar_m),norm(Mline));%2*Mbar_M - Mbar_m;% 
% mu1 = mysolve(2*Mbar_m,-0.5*TM*sqrt(MLine_M),mu2*TM^2,0.1); 
% 
% P1 = [2*Mbar_m -0.5*TM*sqrt(MLine_M);
%       -0.5*TM*sqrt(MLine_M) (mu1+mu2*TM^2)/2];
% epsion = min(1,min(eig(P1))/2)
% 
% kRTemp(1) = 1/(sqrt(mu1*mu2));
% kRTemp(2) = 2*Tm/(mu2*Tm^2+mu1);
% kRTemp(3) = 2*TM/(mu2*TM^2+mu1);
% kR = 0.9*min(kRTemp)
% k_eta = 0;

 
% ee    = 0.9* Mbar_m/Mbar_M % 
% minP1 = ee* 2*Mbar_M; %/Mbar_M
% mu    =  mysolve2(2*Mbar_m,-0.5*TM*sqrt(MLine_M),minP1);
% kR = 0.3/Mbar_m %15*TM/(mu*exp(Tm))  %0.9*TM/(mu*exp(Tm))
% k_eta = 0;
% P1 = [2*Mbar_m -0.5*TM*sqrt(MLine_M);
%       -0.5*TM*sqrt(MLine_M) mu];
% epsion = min(eig(P1)/(2*Mbar_M))


norm_M = norm(M,'fro');
kR     = 0.8*2/(TM*norm_M) 
k_eta  = 0;

n     = size(pI,2);
C     = zeros(3*n,9);
for ii=1:n
    C(3*ii-2:3*ii,:) = [Skew(pI(:,ii)) zeros(3) -eye(3)];
end

% Run time
N = 200;
RunT = Tm + (TM-Tm).*rand(N,1);

% initial conditions
R    = eye(3);
Q    = rotm2quat(R)';
% p    = 10*[1 0 1]'; % circle
% v    = 10*w*[0 1 0]';
p    = 10*[0 0 1]'; % 8-shpae
v    = 10*w*[1 1 0]';


u     =  randn(3,1); %[0.0900 0.6676 -0.7391]';  %
u     = u/norm(u);
Su    = Skew(u);
Rhat  = expm(0.1*pi*Su);
% vRhat = reshape(Rhat,9,1);
Qhat  = rotm2quat(Rhat)';
phat  = 0*p;%zeros(3,1);
vhat  = 0*v;%zeros(3,1);
eta   = zeros(3,1);

sqrt(trace(eye(3)-R*Rhat'))/2

xIn   = [Q;p;v;Qhat;phat;vhat;eta];
 
tic
%%%%%%%%%% 
[tout1,error1,phatout1]=HINO1_F(xIn); 
% disp('HINO1_F completed.') 
toc
tic
%%%%%%%%%% 
[tout2,error2,phatout2]=HINO1_V(xIn);
% disp('HINO1_V completed.')
toc
tic
%%%%%%%%% 
[tout3,error3,phatout3]=HINO2_F(xIn);
% disp('HINO2_F completed.')
tic
toc
%%%%%%%%% 
[tout4,error4,phatout4]=HINO2_V(xIn);
% disp('HINO2_V completed.')
tic
toc
%%%%%%%%%%%%
[tout5,error5,phatout5]=IEKF(xIn);
% disp('IEKF completed.')
toc


% save('output.mat','tout1','tout2','tout3','tout4','tout5',...
%      'error1','error2','error3','error4','error5','pI','w',...
%      'phatout1','phatout2','phatout3','phatout4','phatout5')
 
%%

figure

% p  = 10*[sin(w*t) sin(w*t)*cos(w*t) 1]' % 8-shape
plot3(10*sin(w*tout1),10*sin(w*tout1).*cos(w*tout1),10*ones(size(tout1)),'k--'), hold on
plot3(pI(1,:),pI(2,:),pI(3,:),'k*')
% p  = 10*[cos(w*t) sin(w*t)-1 1]' % circle
% plot3(10*cos(w*tout1),10*sin(w*tout1),10*ones(size(tout1)),'k--')
plot3(phatout1(:,1),phatout1(:,2),phatout1(:,3),'-','linewidth',1) 
plot3(phatout2(:,1),phatout2(:,2),phatout2(:,3),'-','linewidth',1)
plot3(phatout3(:,1),phatout3(:,2),phatout3(:,3),'-','linewidth',1)
plot3(phatout4(:,1),phatout4(:,2),phatout4(:,3),'-','linewidth',1)
plot3(phatout5(:,1),phatout5(:,2),phatout5(:,3),'-.','linewidth',1)
grid on 
set(gca,'GridLineStyle',':','GridAlpha',0.8)
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
% zlim([0 12])
legend('True trajectory','Landmarks','HINO1-F','HINO1-V','HINO2-F','HINO2-V','IEKF')

% set(gcf, 'Renderer', 'Painters');
% print('-depsc','E:\Dropbox (Personal)\Research Note\6-INS\HINO with ILM\Trajectories.eps')
 
%%
figure 
subplot(2,2,1)
plot(tout1,error1(:,1),'-','linewidth',1), hold on
plot(tout2,error2(:,1),'-','linewidth',1),
plot(tout3,error3(:,1),'-','linewidth',1),
plot(tout4,error4(:,1),'-','linewidth',1),  
plot(tout5,error5(:,1),'-.','linewidth',1),  
xlabel('time(s)')
ylabel('$|\tilde{R}|_I$','FontSize',12,'interpreter','latex')
legend('HINO1-F','HINO1-V','HINO2-F','HINO2-V','IEKF')
grid on 
set(gca,'GridLineStyle',':','GridAlpha',0.8)

subplot(2,2,2)
plot(tout1,error1(:,2),'-','linewidth',1), hold on
plot(tout2,error2(:,2),'-','linewidth',1), hold on
plot(tout3,error3(:,2),'-','linewidth',1),
plot(tout4,error4(:,2),'-','linewidth',1) 
plot(tout5,error5(:,2),'-.','linewidth',1), 
xlabel('time(s)')
ylabel('$\|p-\hat{p}\|$','FontSize',12,'interpreter','latex')
% legend('HINO1-F','HINO1-V','HINO2-F','HINO2-V')
grid on 
set(gca,'GridLineStyle',':','GridAlpha',0.8)

subplot(2,2,3)
plot(tout1,error1(:,3),'-','linewidth',1), hold on
plot(tout2,error2(:,3),'-','linewidth',1),
plot(tout3,error3(:,3),'-','linewidth',1),
plot(tout4,error4(:,3),'-','linewidth',1),  
plot(tout5,error5(:,3),'-.','linewidth',1),
xlabel('time(s)')
ylabel('$\|v-\hat{v}\|$','FontSize',12,'interpreter','latex')
% legend('HINO1-F','HINO1-V','HINO2-F','HINO2-V')
grid on 
set(gca,'GridLineStyle',':','GridAlpha',0.8)

subplot(2,2,4)
plot(tout1,error1(:,4),'-',tout2,error2(:,4),'-','linewidth',1), hold on
xlabel('time(s)')
ylabel('$\|\mathrm{g}-\hat{\mathrm{g}}\|$','FontSize',12,'interpreter','latex')
% legend('HINO1-F','HINO1-V','HINO2-F','HINO2-V')
grid on 
set(gca,'GridLineStyle',':','GridAlpha',0.8)


% set(gcf, 'Renderer', 'Painters');
% print('-depsc','E:\Dropbox (Personal)\Research Note\6-INS\HINO with ILM\SimulationError.eps')
 
function [tout,error,phatout]=HINO2_F(xIn) 
global  RunT

ts       = 0;
error    = []; 
tout     = []; 
phatout  = [];
for jj=1:length(RunT)
    
%     tspan = [ts ts+RunT(jj)];
%     
%     % flows
%     options = odeset('RelTol',1e-3,'MaxStep',.1);
%     [t,x] = ode45(@flow,tspan,xIn,options);
    tspan = ts :0.01: ts+RunT(jj);
    x = ode4(@flow,tspan,xIn);
    t = tspan';
    % jumps
    xIn = jump(x(end,:)');
    
    x(end,:) = xIn';
    ts = t(end);
    
    % plots 
    
    errorR = zeros(size(t));
    errorp = zeros(size(t));
    errorv = zeros(size(t));
    for ll=1:length(t)
        ii    = 0;
        Q     = x(ll,ii+1:ii+4);
        p     = x(ll,ii+5:ii+7);
        v     = x(ll,ii+8:ii+10);
        ii    = ii+ 10;
        Qhat = x(ll,ii+1:ii+4);
        phat  = x(ll,ii+5:ii+7);
        vhat  = x(ll,ii+8:ii+10);
%         ii    = ii+ 10;
%         eta   = x(ll,ii+1:ii+3)';

        
        R     = quat2rotm(Q);
        Rhat  = quat2rotm(Qhat);

        errorR(ll,1) =  sqrt(trace(eye(3)-R*Rhat'))/2;
        errorp(ll,1) = norm(p-phat);
        errorv(ll,1) = norm(v-vhat);

    end 
    error = [error; errorR errorp errorv];
    tout  = [tout;t]; 
     phatout  = [phatout;x(:,15:17)]; 
end
end


function dx=flow(t,x)
global pC   g k_eta

ii    = 0;
Q    = x(ii+1:ii+4,1);
p     = x(ii+5:ii+7,1);
v     = x(ii+8:ii+10,1);
ii    = ii+ 10;
Qhat = x(ii+1:ii+4,1);
phat  = x(ii+5:ii+7,1);
vhat  = x(ii+8:ii+10,1);
ii    = ii+ 10;
eta   = x(ii+1:ii+3,1); 

R     = quat2rotm(Q');
Rhat  = quat2rotm(Qhat'); 
 
%IMU  
[a,omega,a_y,omega_y] = IMU(R,t);
 
% dR    = R*Skew(omega);
Qw    = [0;omega];
dQ    = 0.5*quatmultiply(Q',Qw')';
dp    = v;
dv    = g + R*a; 

% dRhat  = Rhat*Skew(omega+Rhat'*eta);
Qw    = [0;omega_y+Rhat'*eta];
dQhat    = 0.5*quatmultiply(Qhat',Qw')';
dphat  = Skew(eta)*(phat-pC) + vhat;
dvhat  = Skew(eta)* vhat + g + Rhat*a_y;
deta   = -k_eta*eta;%zeros(3,1); 
dx=[dQ;dp;dv;dQhat;dphat;dvhat;deta];
end

function xplus = jump(x)
global kn pC  pI kR ny
ii    = 0;
Q    = x(ii+1:ii+4,1);
p     = x(ii+5:ii+7,1);
v     = x(ii+8:ii+10,1);
ii    = ii+ 10;
Qhat = x(ii+1:ii+4,1);
phat  = x(ii+5:ii+7,1);
vhat  = x(ii+8:ii+10,1);
% ii    = ii+ 10;
% eta   = x(ii+1:ii+3,1); 

R     = quat2rotm(Q');
Rhat  = quat2rotm(Qhat'); 

% gains
kp = 0.5;
kv = 1.0; 

%landmark measurements
y = R'*(pI-p) + ny*randn(size(pI));

% sigmaR = zeros(3,1);
% sigmap = zeros(3,1);
% for ii=1:length(pI)
%     temp = pI(:,ii)-phat-Rhat*y(:,ii);
%     sigmaR = sigmaR + 0.5*kR* kn(ii,ii)*Skew(pI(:,ii)-pC)*temp;
%     sigmap = sigmap + kn(ii,ii)*temp;
% end 
Temp = pI-phat-Rhat*y;
sigmap = sum(Temp*kn,2);
Temp = (Temp*kn*(pI-pC)');
sigmaR = 0.5*kR*[Temp(3,2)-Temp(2,3),Temp(1,3)-Temp(3,1),Temp(2,1)-Temp(1,2)]'; 


phat = phat+ kp*sigmap/(sum(diag(kn)));
vhat = vhat+ kv*sigmap/(sum(diag(kn)));
eta  = sigmaR;

xplus = [Q;p;v;Qhat;phat;vhat;eta];

end










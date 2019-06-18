function [tout,error,phatout]=IEKF(xIn)

global  RunT

xIn = xIn(1:end-3,1);
P = eye(9);
vP = reshapeT(P);
xIn = [xIn;vP];

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
global  g sigmaa sigmag
ii    = 0;
Q    = x(ii+1:ii+4,1);
p     = x(ii+5:ii+7,1);
v     = x(ii+8:ii+10,1);
ii    = ii+ 10;
Qhat = x(ii+1:ii+4,1);
phat  = x(ii+5:ii+7,1);
vhat  = x(ii+8:ii+10,1); 
vP    = x(ii+11:end,1); 

R     = quat2rotm(Q');
Rhat  = quat2rotm(Qhat');  
P     = reshapeT(vP);


%IMU  
[a,omega]=IMU(R,t);
 
% dR    = R*Skew(omega);
Qw    = [0;omega];
dQ    = 0.5*quatmultiply(Q',Qw')';
dp    = v;
dv    = g + R*a; 


A     = [zeros(3,9);Skew(g) zeros(3,6);zeros(3) eye(3) zeros(3)];
CovW    = diag([sigmag.^2*ones(1,3), sigmaa.^2*ones(1,3) zeros(1,3)]);
Temp    = [Rhat zeros(3,6);Skew(vhat)*Rhat Rhat zeros(3);
           Skew(phat)*Rhat zeros(3) Rhat];
Vt      = Temp*CovW*Temp'; %0.01*eye(9);%

% dRhat  = Rhat*Skew(omega);
Qw    = [0;omega];
dQhat    = 0.5*quatmultiply(Qhat',Qw')';
dphat  = vhat;
dvhat  = g + Rhat*a;
dP     = A*P + P*A' + Vt; 
dvP    = reshapeT(dP);

dx=[dQ;dp;dv;dQhat;dphat;dvhat;dvP];
end

function xplus = jump(x)
global  pI  ny 
ii    = 0;
Q    = x(ii+1:ii+4,1);
p     = x(ii+5:ii+7,1);
v     = x(ii+8:ii+10,1);
ii    = ii+ 10;
Qhat = x(ii+1:ii+4,1);
phat  = x(ii+5:ii+7,1);
vhat  = x(ii+8:ii+10,1); 
vP    = x(ii+11:end,1); 

R     = quat2rotm(Q');
Rhat  = quat2rotm(Qhat'); 
P     = reshapeT(vP);



%landmark measurements
y = R'*(pI-p) + ny*randn(size(pI));


n     = size(pI,2);
C     = zeros(3*n,9);
for ii=1:n
    C(3*ii-2:3*ii,:) = [Skew(pI(:,ii)) zeros(3) -eye(3)];
end
Qt     = ny^2*eye(3*n);
Ln     = P*C'/(C*P*C'+Qt);
P      = (eye(size(P))-Ln*C)*P;

Xhat   = [Rhat vhat phat;zeros(2,3) eye(2)];
r      = [pI;zeros(1,n);ones(1,n)];
b      = [y;zeros(1,n);ones(1,n)];
II     = [eye(3) zeros(3,2)];
Temp   = II*(Xhat*b-r);
Temp    = Ln*reshape(Temp,3*n,1);
Delta  = [Skew(Temp(1:3,1)) Temp(4:6,1) Temp(7:9,1);zeros(2,5)];
Xhat   = expm(Delta)*Xhat;

% vRhat  = reshape(Xhat(1:3,1:3),9,1);
vhat   = Xhat(1:3,4);
phat   = Xhat(1:3,5);
vP     = reshapeT(P); 

Qhat   = rotm2quat(Xhat(1:3,1:3))';

xplus = [Q;p;v;Qhat;phat;vhat;vP];

end

function [a,omega,a_y,omega_y] = IMU(R,t)
global w g sigmaa sigmag

omega = 1*[sin(0.3*pi*t) 0.1 cos(0.3*pi*t)]' ;
omega_y = omega + 0*sigmag*randn(3,1);
%%% - circle
% p  = 10*[cos(w*t) sin(w*t)-1 1]'
% dp = 10*w *[-sin(w*t) cos(w*t) 0];
% dv    = -10*w^2*[cos(w*t) sin(w*t) 0]'+ 0*sigmaa*randn(3,1);

%%% - 8 shape
% p  = 10*[sin(w*t) sin(w*t)*cos(w*t) 1]'
% dp = 10*w*[cos(w*t) (cos(w*t))^2-(sin(w*t))^2 0]'
dv    = -10*w^2*[sin(w*t) 4*sin(w*t)*cos(w*t) 0]' ;


a     = R'*(dv-g); 
a_y   = a+  0*sigmaa*randn(3,1);

end
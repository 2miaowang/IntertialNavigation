clear all
clc

%%%%% Intermittency %%%%%%
Tm = .1;
TM = .2;

% %%%%% Gains %%%%%%
kp = 0.5;
kv = 1.0;
kg = 0.6; 


% %%%%% Plant1 %%%%%%
% Af = [zeros(3) eye(3) zeros(3); zeros(3) zeros(3) eye(3); zeros(3,9)];
% C = [eye(3) zeros(3,6)];
% K = [kp*eye(3);kv*eye(3);kg*eye(3)]; 
% % % 


%%%%% Gains %%%%%%
% kp = 0.6;
% kv = 1.0; 
%%%% Plant2 %%%%%%
Af = [zeros(3) eye(3); zeros(3) zeros(3)];
C = [eye(3) zeros(3,3)];
K = [kp*eye(3);kv*eye(3)]; 


 
Ag = eye(size(Af)) - K*C;
  
SolverP(Af,Ag,Tm,TM)
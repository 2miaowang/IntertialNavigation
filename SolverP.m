function SolverP(Af,Ag,Tm,TM)

n = length(Af);
ba = 1;
cvx_begin sdp
    variable II(n,n) symmetric  
    (-Af + ba*eye(n))'*II + II*(-Af + ba*eye(n))>= eye(n);
    II >= eye(n);
cvx_end
gamma = sqrt(max(eig(II))/min(eig(II)));

eA  = gamma*exp(ba*TM);

%%%%%%%%% SOLVE FOR P %%%%%%%%%%%
 
mu = 1; 
T = [Tm TM];
while (1)
 
    [pM,P]=mOptfun(Af,Ag,T,mu);
    dT = mu/(2*eA*pM*norm(Af)*norm(Ag));
    Td = Tm:2*dT:TM; 
    if isnan(pM)
%         clc
        disp('Status: Infeasible')
        break;
    end
    MinE = zeros(size(Td));
    for i=1:length(Td)
        H = expm(Af*Td(i))*Ag; 
        MinE(i) = max(eig(H'*P*H - P));  
    end
    if (max(MinE)<-mu) 
%         clc
        disp('Status: Solved')
        disp(['max(H^T*P*H - P)=' num2str(max(MinE))]) 
        break; 
    elseif (max(MinE)>=-mu)&&(length(T)<length(Td)) 
        [~,index] = max(MinE)
        T = [T; Td(index)];
    else
%         clc
        disp('Status: Infeasible')
        break;
    end
     
end

% plot(Td,MinE) 
end

% solve P and pM
function [pM,P]=mOptfun(Af,Ag,T,mu)

n = length(Af);

cvx_begin sdp
    variable P(n,n) symmetric
    variable pM
    minimize (pM)    
    for i=1:length(T)
        (expm(Af*T(i))*Ag)'*P*expm(Af*T(i))*Ag - P <= -2*mu*eye(n);
    end 
    P >= eye(n);
    P <= pM*eye(n);
cvx_end

end


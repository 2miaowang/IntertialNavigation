function Y=reshapeT(A)
    [c,r]=size(A);
%     square to vector
    if c==r
        Y  = zeros(c*(c+1)/2,1);
        id = 0;
        for i=1:c 
            Y(id+1:id+c-i+1,1) = A(i:c,i);
            id = id + (c-i+1);
        end
%     vector to square    
    elseif r==1
        n  = (-1+sqrt(1+8*c))/2;
        Y  = zeros(n);
        id = 0;
        for i=1:n
            Y(i:n,i) = A(id+1:id+n-i+1,1); 
            id = id + (n-i+1);
        end
        Y = (Y + Y')-diag(diag(Y));        
    else        
        disp('Dimensional eror of input')
    end

end
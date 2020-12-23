function [x,y,n] = StartingPoints(N,a,b)
     
    if mod(N,2) == 0
        N = N + 1;
    end
    
    if N == 1 
         x = 0;  y = 0;    
         n = 1;
        return;
    end
    
    M = ceil(2*(N-1)/sqrt(3)); 
    if mod(M,2) == 0
        M = M + 1;
    end
    
    x = zeros(1,N*M);
    y = zeros(1,N*M);
    
    centerline = (M+1)/2;
    count = 0;
    for j = 0:(M-1)
        for i = 0:(N-1)
            X = - a + 2*i*a/(N-1);
            Y = - (M-1)*sqrt(3)*a/(2*(N-1)) + j*sqrt(3)*a/(N-1);
                      
            if mod(j - centerline,2) == 0
                    X = X + a/(N-1);
            end
            
            if sqrt(X^2+Y^2) - a <= 10e-6                
                count = count + 1;
                x(count) = X;
                y(count) = b*Y/a;
            end     
        end
    end
    n = count;
    x = x(1:n);
    y = y(1:n); 
end    

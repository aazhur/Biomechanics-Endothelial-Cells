function phi = angle_to_neighbor(xc,yc,nc,R)
    X = ones(nc,nc)*diag(xc)-diag(xc)*ones(nc,nc);
    Y = ones(nc,nc)*diag(yc)-diag(yc)*ones(nc,nc); 
    
    Dist = sqrt(X.^2+Y.^2);
    Dist(Dist == 0) = R*2;

    [~,i] = min(Dist,[],2);
    i = i' + linspace(0,nc*(nc-1),nc);
  
    phi = mod(atan2d(Y(i),X(i))+360,360);
    disp(phi);
    
end
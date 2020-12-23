function [prot_len,es,ae] = prot_lengths(es,xc,yc,np,nd,nc,dl,prot_len,begx,begy,endx,endy,ca,ba,dla,R,effective_r,close_r,effective_n,close_n) 
 
    es = es.*[zeros(nc,np),ones(nc,nd)];
    ae = zeros(nc,np);
    [Clust_far,Clust_close,Box_clust_far,Box_clust_close,Intersec] = cluster(xc,yc,nc,effective_r,close_r,R); 
        
    if sum(Box_clust_far) ~= 0
        idx = find(Box_clust_far);
        for i = 1:numel(idx)
            number = sum(ba(idx(i),1:np) > 0)+sum((ca(idx(i),1:np) > 0).*ismember(ca(idx(i),1:np),Clust_far(idx(i),:)));
            if number < effective_n
                cells_dx = xc(idx(i)); 
                cells_dy = yc(idx(i));

                prot_dx = (begx(idx(i),1:np)-xc(idx(i))).*(prot_len(idx(i),1:np)==0)+(endx(idx(i),1:np)-begx(idx(i),1:np)).*(prot_len(idx(i),1:np)>0);
                prot_dy = (begy(idx(i),1:np)-yc(idx(i))).*(prot_len(idx(i),1:np)==0)+(endy(idx(i),1:np)-begy(idx(i),1:np)).*(prot_len(idx(i),1:np)>0); 

                angles = atan2d(cells_dx*prot_dy-cells_dy*prot_dx,cells_dx*prot_dx+cells_dy*prot_dy);
                att = ((-180/np) < angles & angles < (180/np)).*(ba(idx(i),1:np) == 0); 
                if sum(att) > (effective_n-number)
                    ind = find(att);
                    att(ind(effective_n+1-number:end)) = 0;
                end   
                es(idx(i),1:np) = es(idx(i),1:np) + att;
            end
        end
    end      
    
    if sum(Box_clust_close) ~= 0
        idx = find(Box_clust_close);
        for i = 1:numel(idx)
            number = sum((ba(idx(i),1:np)+ca(idx(i),1:np)) > 0);
            if number < close_n
                cells_dx = xc(idx(i)); 
                cells_dy = yc(idx(i));

                prot_dx = (begx(idx(i),1:np)-xc(idx(i))).*(prot_len(idx(i),1:np)==0)+(endx(idx(i),1:np)-begx(idx(i),1:np)).*(prot_len(idx(i),1:np)>0);
                prot_dy = (begy(idx(i),1:np)-yc(idx(i))).*(prot_len(idx(i),1:np)==0)+(endy(idx(i),1:np)-begy(idx(i),1:np)).*(prot_len(idx(i),1:np)>0);

                angles = atan2d(cells_dx*prot_dy-cells_dy*prot_dx,cells_dx*prot_dx+cells_dy*prot_dy);
                att = ((-180/np) < angles & angles < (180/np)).*(ba(idx(i),1:np) == 0); 
                if sum(att) > (close_n-number)
                    ind = find(att);
                    att(ind(close_n+1-number:end)) = 0;
                end   
                es(idx(i),1:np) = es(idx(i),1:np) + att;
            end
        end
    end
    
    for i=1:nc        
        number = sum((ba(i,1:np)+ca(i,1:np)) > 0);
        if number < close_n
        if sum(Clust_close(i,:)) > 0
            Idx = find(Clust_close(i,:));
            [~,~,att_cells] = find(ca(i,:));
            [~,~,idx] = find(Idx.*(1-ismember(Idx,att_cells)));
            if ~isempty(idx)
                [~,I] = sort(Intersec(i,idx));
                idx = idx(I);
                
                cells_dx = diag(xc(idx)-xc(i))*ones(numel(idx),np); 
                cells_dy = diag(yc(idx)-yc(i))*ones(numel(idx),np);

                prot_dx = ones(numel(idx),np)*diag((begx(i,1:np)-xc(i)).*(prot_len(i,1:np)==0)+(endx(i,1:np)-begx(i,1:np)).*(prot_len(i,1:np)>0));
                prot_dy = ones(numel(idx),np)*diag((begy(i,1:np)-yc(i)).*(prot_len(i,1:np)==0)+(endy(i,1:np)-begy(i,1:np)).*(prot_len(i,1:np)>0));

                cprot_dx = begx(idx,1:np)-xc(i);
                cprot_dy = begy(idx,1:np)-yc(i);

                A = atan2d(cells_dx.*cprot_dy-cells_dy.*cprot_dx,cells_dx.*cprot_dx+cells_dy.*cprot_dy); 
                A1 = diag(min(A.*(A <= 0),[],2))*ones(numel(idx),np);
                A2 = diag(max(A.*(A >= 0),[],2))*ones(numel(idx),np);  
                angles = atan2d(cells_dx.*prot_dy-cells_dy.*prot_dx,cells_dx.*prot_dx+cells_dy.*prot_dy); 

                Dist = (cells_dx.^2+cells_dy.^2)-ones(numel(idx),np)*diag(prot_len(i,1:np).^2); 
                att = (angles.*(A1 <= angles & angles <= A2).*(Dist > 0)).*(ones(numel(idx),np)*diag((ca(i,1:np) == 0)));
                [row,col] = find(att);
                [col,gamma,~] = unique(col,'stable');
                row = row(gamma);
                if numel(col) > (close_n-number)
                    col = col(1:(close_n-number));
                    row = row(1:(close_n-number));
                end
                att = zeros(1,np); where = zeros(1,np);
                att(col) = 1; where(col) = idx(row);
                ae(i,:) = ae(i,:) + where;
                es(i,1:np) = es(i,1:np)+att;
            end
        end
        end
                
        number = sum(ba(i,1:np) > 0)+sum((ca(i,1:np) > 0).*ismember(ca(i,1:np),Clust_far(i,:)));
        if number < effective_n
        if sum(Clust_far(i,:)) > 0
            Idx = find(Clust_far(i,:));
            [~,~,att_cells] = find(ca(i,:));
            [~,~,idx] = find(Idx.*(1-ismember(Idx,att_cells)));
            if ~isempty(idx)
                if numel(idx) > (effective_n-number)
                    [~,I] = sort(Intersec(i,idx));
                    idx = idx(I(1:(effective_n-number)));
                end    
                cells_dx = diag(xc(idx)-xc(i))*ones(numel(idx),np); 
                cells_dy = diag(yc(idx)-yc(i))*ones(numel(idx),np);

                prot_dx = ones(numel(idx),np)*diag((begx(i,1:np)-xc(i)).*(prot_len(i,1:np)==0)+(endx(i,1:np)-begx(i,1:np)).*(prot_len(i,1:np)>0));
                prot_dy = ones(numel(idx),np)*diag((begy(i,1:np)-yc(i)).*(prot_len(i,1:np)==0)+(endy(i,1:np)-begy(i,1:np)).*(prot_len(i,1:np)>0));

                cprot_dx = begx(idx,1:np)-xc(i);
                cprot_dy = begy(idx,1:np)-yc(i);

                A = atan2d(cells_dx.*cprot_dy-cells_dy.*cprot_dx,cells_dx.*cprot_dx+cells_dy.*cprot_dy); 
                A1 = diag(min(A.*(A <= 0),[],2))*ones(numel(idx),np);
                A2 = diag(max(A.*(A >= 0),[],2))*ones(numel(idx),np);  
                angles = atan2d(cells_dx.*prot_dy-cells_dy.*prot_dx,cells_dx.*prot_dx+cells_dy.*prot_dy); 

                Dist = (cells_dx.^2+cells_dy.^2)-ones(numel(idx),np)*diag(prot_len(i,1:np).^2); 
                A = (A1 <= angles & angles <= A2);
                angles = (angles.*A.*(Dist > 0)).*(ones(numel(idx),np)*diag((ca(i,1:np) == 0)));
                angles(A == 0) = NaN;
                medians = median(angles,2,'omitnan');
                medians(isnan(medians)) = 0;
                M = angles.*(angles >= diag(medians)*ones(numel(idx),np));
                M(M == 0) = NaN;
                M = min(M,[],2);
                medians = diag(M)*ones(numel(idx),np);
                att = (angles == medians);   
                [row,col] = find(att);
                [col,gamma,~] = unique(col,'stable');
                row = row(gamma);
                where = zeros(1,np);     
                
                if numel(idx) > 1
                    att = sum(att);
                end 
                where(col) = idx(row);
                ae(i,:) = ae(i,:).*(where == 0) + where;
                es(i,1:np) = es(i,1:np)+(att > 0);
            end
        end
        end
    end
    
    es = es > 0;
    am = (ca ~= 0 | ba ~= 0);
    lp = prot_len;
    
    %springs = ((endx(:,1:np)-spx(:,1:np)).^2+(endy(:,1:np)-spy(:,1:np)).^2);
    %lp = lp+dla*[(springs <= 4),zeros(nc,nd)].*am.*(2*es-1)+dl*(1-am).*(2*es-1);  
    lp = lp+dla*[ones(nc,np),zeros(nc,nd)].*am.*(2*es-1)+dl*[ones(nc,np),0.21*ones(nc,nd)].*(1-am).*(2*es-1); 
    prot_len = (lp+abs(lp))/2;
    
end   

function [Clust_far,Clust_close,Box_clust_far,Box_clust_close,Intersec] = cluster(xc,yc,nc,r,close,R)      
    Intersec = R - sqrt(xc.^2+yc.^2);
    Box_clust_far = (Intersec > close & Intersec <= r);
    Box_clust_close = (Intersec <= close);  
    
    X = (diag(xc)*ones(nc,nc)-ones(nc,nc)*diag(xc)).^2;
    Y = (diag(yc)*ones(nc,nc)-ones(nc,nc)*diag(yc)).^2; 
    Intersec = sqrt(X+Y);
    Clust_far = (Intersec > close & Intersec <= r).*(ones(nc,nc)*diag(linspace(1,nc,nc)));
    Clust_far(logical(eye(size(Clust_far)))) = 0;
    
    Clust_close = (Intersec <= close).*(ones(nc,nc)*diag(linspace(1,nc,nc)));
    Clust_close(logical(eye(size(Clust_close)))) = 0;   
end
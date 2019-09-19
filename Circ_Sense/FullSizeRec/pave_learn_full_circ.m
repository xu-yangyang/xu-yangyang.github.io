function U = pave_learn_full_circ(D,N1,N2,n1,n2,K)
nr = floor(N1/n1); nc = floor(N2/n2);
resnr = N1-nr*n1; resnc = N2-nc*n2;
N = nr*nc+sign(resnr)*nc+sign(resnc)*nr+sign(resnr*resnc);
Y = zeros(N1*N2,N*K);

k = 0;
for j = 1:nc
    for i = 1:nr
        y = zeros(N1,N2);
        start_row = (i-1)*n1+1; start_col = (j-1)*n2+1;
        end_row = i*n1; end_col = j*n2;
        for dnum = 1:K
            k = k+1;
            y(start_row:end_row,start_col:end_col) = ...
                reshape(D(:,dnum),n1,n2);
            Y(:,k) = reshape(ifft2(y),N1*N2,1);
        end
    end
end

if resnr~=0
    for j = 1:nc
        start_row = nr*n1+1;
        start_col = (j-1)*n2+1; end_col = j*n2;
        y = zeros(N1,N2);
        for dnum = 1:K
            k = k+1;
            xi = reshape(D(:,dnum),n1,n2);
            y(start_row:end,start_col:end_col) = xi(end-resnr+1:end,:);
            Y(:,k) = reshape(ifft2(y),N1*N2,1);
        end
    end
end
if resnc~=0
    for i = 1:nr
        start_col = nc*n2+1;
        start_row = (i-1)*n1+1; end_row = i*n1;
        y = zeros(N1,N2);
        for dnum = 1:K
            k = k+1;
            xi = reshape(D(:,dnum),n1,n2);
            y(start_row:end_row,start_col:end) = xi(:,end-resnc+1:end);
            Y(:,k) = reshape(ifft2(y),N1*N2,1);
        end
    end
end
if resnr~=0 && resnc~=0
    start_row = nr*n1+1; start_col = nc*n2+1;
    y = zeros(N1,N2);
        for dnum = 1:K
            k = k+1;
            xi = reshape(D(:,dnum),n1,n2);
            y(start_row:end,start_col:end) =...
                xi(end-resnr+1:end,end-resnc+1:end);
            Y(:,k) = reshape(ifft2(y),N1*N2,1);
        end
end
Y = Y*sqrt(N1*N2);
Ys = Y*Y'; Ys = (Ys+Ys')/2;
H = conj(Ys).*Ys; f = diag(Ys);
quadopts = optimset('Display','off');
u = quadprog(H,-f,[],[],[],[],zeros(N1*N2,1),[],[],quadopts);
U = reshape(u,N1,N2);
function seq = zadoff(root,N) 

if mod(N,2)==1
for n=0:(N-1)
    seq(n+1)=exp(-j*(pi*root*n*(n+1))/N);
end
elseif mod(N,2)==0
    for n=0:(N-1)
        seq(n+1)=exp(-j*(pi*root*n*(n))/N);
    end
end



% nn=16;max(max(abs(conj(dftmtx(nn))*diag(zadoff(1,nn))*dftmtx(nn)*dctmtx(nn))/sqrt(nn)))
% nn=16;max(max(abs(conj(dftmtx(nn))*diag(zadoff(1,nn))*dftmtx(nn)/sqrt(nn)))
% nn=16;max(max(abs(conj(dftmtx(nn))*diag(zadoff(1,nn))*dftmtx(nn)/sqrt(nn)
% )
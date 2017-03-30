% Gold code generation, period 31
% order-5 primitive polynomial 1 + x^3 + x^5 (L = 5)
% the preferred pair is obtained by decimating the m sequence
% decimation rule: r = 2^[(L+2)/2] + 1 = 9;
% extended to period 32, Mar. 2009
function CS_Matrix = Goldcode32x1024

clear;
L = 5;
N = 2^L-1;
r = 2^floor((L+2)/2) + 1;
mseqa=zeros(N,1);
mseqb=zeros(N,1);  %preferred pair
u=zeros(N,1);
v=zeros(N,1);
v_shift=zeros(N,1);
v_extend=zeros(N,1);
v_extend_shift=zeros(N,1);
SM=[];

x=ones(L,1);
for i=1:N
   xx=xor(x(5),x(3)); 
   for k=L:-1:2
     x(k)=x(k-1);
   end
   x(1)=xx;
   mseqa(i)=(-1)^xx;
end

for i=1:N-1
  k=rem(r*i,N);
  mseqb(i)=mseqa(k);
end
mseqb(N)=mseqa(N);

SM=[mseqa mseqb];

for shift = 0:N-1
  for k=1:N
     u(k)=mseqa(k);
     if k + shift <= N  % get a code sequence v
        v(k)=mseqb(k+shift) * mseqa(k);
     else
        v(k)=mseqb(rem(k+shift,N+1)+1) * mseqa(k);
     end
  end
  SM = [SM v];
end

CS_Matrix = [];
for n=1:N+2
  v_extend = SM(:,n);
  CS_Matrix = [CS_Matrix v_extend];
  for kk=1:N-1      % shift a code sequence v
      for k=1:N-1
          v_extend_shift(k+1)=v_extend(k);
      end
      v_extend_shift(1)=v_extend(N);
      v_extend=v_extend_shift;
      CS_Matrix = [CS_Matrix v_extend];
  end
end

CS_Matrix = [ones(N,1) CS_Matrix];
CS_Matrix = [ones(1,(N+1)*(N+1)); CS_Matrix];
% CS_Matrix = [prod(CS_Matrix); CS_Matrix];

% u=SM(:,3);
% v=SM(:,12);
% u=CS_Matrix(:,10);
% v=CS_Matrix(:,456);
% fa=fft(u); fb=fft(v);
% ccf=ifft(fa.*conj(fb));
% eccf=abs(ccf);
%  max(eccf)


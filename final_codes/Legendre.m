function Lege = Legendre(p,N); 

%p=7; % p=3 mod 4
%N=7;

Lege = zeros(1,N);
for k = 1:N
Lege(k) = mod((k-1)^((p-1)/2), p);
if Lege(k)== p-1;
    Lege(k)= -1;
end

end
Lege(1)=1;
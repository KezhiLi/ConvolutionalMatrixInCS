function Lege = Legendre1(p,N); 

%p=7; % p=3 mod 4
%N=7;
% 
% Lege = -ones(1,N);
% for k = 0:N-1
%     if ~(isempty(find(mod(k , p)==[1,4,9,16,25,36,49,64,81,100])))
%     Lege(k+1) = 1;
%     elseif ~(isempty(find(mod(k , p)==0)))
%     Lege(k+1) = 0;
%     end
% end
sqr=[1:1024].^2;
Lege = -ones(1,N);
for k = 0:N-1
    test=round(max(sqr)/p);
    candi= k+[1:test].*p;
    for ii=1:length(candi);
        if ~(isempty(find(candi(ii)==sqr)))
            Lege(k+1) = 1;
        end
    end
end

%Lege(1)=1;
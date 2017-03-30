function otf_conjs = conjugate_symmetrize(otf)
% CONJUGATE_SYMMETRIZE  generate a conjugate symmetric OTF
%   function otf_conjs = conjugate_symmetrize(otf)
%  
%   From an arbitrary 2D OTF, it generates a conjuage symmetric OTF.
%   Author: Wotao Yin, wotao.yin@rice.edu. Mar 01, 2010.

cnj = @(N,HF) mod(N-HF+1,N)+1;

[m,n] = size(otf);

if mod(m,2)==0
    mhalf1 = (1:m/2); mhalf2 = (m/2+1:m);
else
    mhalf1 = (1:(m+1)/2); mhalf2 = ((m+1)/2+1:m);
end

if mod(n,2)==0
    nhalf = (1:n/2+1); 
else
    nhalf = (1:(n+1)/2);
end

% method = 2;
% switch method
%     case 1
%         otf_conjs = complex(zeros(m,n),zeros(m,n));
%
%         otf_conjs(mhalf1,nhalf) = otf(mhalf1,nhalf) + conj(otf(cnj(m,mhalf1),cnj(n,nhalf)));
%         otf_conjs(cnj(m,mhalf1),cnj(n,nhalf)) = conj(otf_conjs(mhalf1,nhalf));
%         
%         otf_conjs(mhalf2,nhalf) = otf(mhalf2,nhalf) + conj(otf(cnj(m,mhalf2),cnj(n,nhalf)));
%         otf_conjs(cnj(m,mhalf2),cnj(n,nhalf)) = conj(otf_conjs(mhalf2,nhalf));
%     case 2
        otf_conjs = otf;
        m_real = ((1:m)==cnj(m,(1:m)));
        n_real = ((1:n)==cnj(n,(1:n)));
        if any(m_real) && any(n_real)
            otf_conjs(m_real,n_real) = abs(otf_conjs(m_real,n_real));
        end
        otf_conjs(cnj(m,mhalf1),cnj(n,nhalf)) = conj(otf_conjs(mhalf1,nhalf));
        otf_conjs(cnj(m,mhalf2),cnj(n,nhalf)) = conj(otf_conjs(mhalf2,nhalf));
% end

% --- test ---
%     if ~all(abs(otf_conjs(:)))
%         fprintf('error\n');
%         keyboard;
%     end
%     for ii = 1:m
%         for jj = 1:n
%             if abs(otf_conjs(ii,jj)-conj(otf_conjs(cnj(m,ii),cnj(n,jj))))>eps
%                 fprintf('error\n');
%                 keyboard;
%             end
%         end
%     end
%     tmp = randn(m,n); tmp = ifft2(otf_conjs.*fft2(tmp));
%     if ~isreal(tmp)
%         fprintf('error\n');
%         keyboard;
%     end

end


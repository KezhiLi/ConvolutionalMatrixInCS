function [U,Out] = RecPF_Circ(m,n,aTV,aL1,picks,Mask,B,TVtype,opts,PsiT,Psi,URange,uOrg)
% [U,Out] = RecPF_Circ(m,n,aTV,aL1,picks,B,TVtype,opts,PsiT,Psi,URange,uOrg)
%
% RecPF solves the TVL1-L2 model:
%
%   min aTV*TV(u) + aL1*|PsiT*U|_1 + 0.5|P*C*U - B|_2^2
%
% where P is a selection operator and C=F*diag(d)F is a circulant matrix.
%
% Inputs:
%
%  aTV, aL1 -- regularization parameters in the model (positive)
%           IMPORTANT: see Lines 60, 66, 67 for data/parameter normalization
%  picks    -- sample positions corresponding to P, logical array
%  Mask     -- the kernel matrix, or PSF, that determines C=circulate(Mask)
%  B        -- measurment vector
%  TVtype   -- 2 (isotropic) or 1 (anisotropic), future versions allow various neighborhood types for anisotropic TV
%  opts     -- parameters for algorithm, see below
%  PsiT     -- sparsifying basis of U, 
%              x = PsiT*U is the sparse representation of U under Psi
%  Psi      -- inverse of PsiT, Psi*x = U reconstructs the image
%  URange   -- grayscale range of U, e.g., 1, 255, 65535
%  uOrg     -- (optional) true image
%
% Outputs: 
%          
%  U      -- reconstructed image
%  Out    -- iteration information, e.g., iter number, relative errors,
%            function values, etc.
%
% Options:
%
%  [to be completed]
%
% This code is written by Wotao Yin following previous code of RecPF 
% and Junfeng Yang's partical convolution code.
%
% Wotao Yin, CAAM, Rice University
% Junfeng Yang, Math, Nanjing University
%
% Feburary 2010.

%% set options
    % --- getting from opts ---
    maxItr = opts.maxItr;
    gamma = opts.gamma;
    beta1 = opts.beta1;
    beta2 = opts.beta2;
    beta3 = opts.beta3;
    bSymm = opts.bsymm;
    relchg_tol = opts.relchg_tol;

    bPrint  = false;    % screen display switch; turning it on slows the code down
    if isfield(opts,'bPrint'); bPrint = opts.bPrint; end
        
    bComplex = true;    % use complex intermediat computation for U?
    if isfield(opts,'bComplex'); bComplex = opts.bComplex; end
    

%% initial vars U, etc.
% If changing U to anything nonzero, you must change the
% initialization of Ux and Uy below
if bComplex; 
    U = complex(zeros(m,n),zeros(m,n));
    g = complex(zeros(m,n),zeros(m,n));
    CircU = complex(zeros(m,n),zeros(m,n));
else
    U = zeros(m,n);
    g = zeros(m,n);
    CircU = zeros(m,n);
end

%% test input
if (aTV < 0); error('aTV must be non-negative'); end
if (aL1 < 0); error('aL1 must be non-negative'); end
if (aL1+aTV <= 0); error('One of aTV and aL1 must be strictly positive'); end
if exist('uOrg','var'); snr(U,uOrg); end % initial SNR computation

%% normalize parameters and data
if opts.normalize
    if ~isreal(URange)||~isscalar(URange)||URange<eps; error('URange must be postive and real.'); end
    fctr = 1/URange;
    B = fctr*B;
    if exist('uOrg','var'); uOrg = fctr*uOrg; snr(U,uOrg); end

    aTV = nnz(picks)/sqrt(m*n)*aTV;
    aL1 = nnz(picks)/sqrt(m*n)*aL1;
end

%% initialize constant parts of numinator and denominator (in order to save computation time)
Numer1 = zeros(m,n); Numer1(picks) = B/beta3;
D = psf2otf(Mask);
Ds = conj(D);
Denom1 = beta3*real(D.*Ds); % Denom1
if aTV > 0
    prd = sqrt(aTV*beta1);
    Denom2 = (abs(psf2otf([prd,-prd],[m,n])).^2 + abs(psf2otf([prd;-prd],[m,n])).^2);
end
if aL1 > 0 
    if aTV > 0; Denom2 = Denom2 + aL1*beta2; else Denom2 = aL1*beta2; end
end
Denom = Denom1 + Denom2;

%% initialize constants
if aTV > 0
    if bComplex
        Ux = complex(zeros(m,n),zeros(m,n)); Uy = complex(zeros(m,n),zeros(m,n));
        bx = complex(zeros(m,n),zeros(m,n)); by = complex(zeros(m,n),zeros(m,n));
    else
        Ux = zeros(m,n); Uy = zeros(m,n);
        bx = zeros(m,n); by = zeros(m,n);
    end
end
if aL1 > 0
    if bComplex
        PsiTU = complex(zeros(m,n),zeros(m,n)); % equal to PsiT(U) as U=0;
        z = complex(zeros(m,n),zeros(m,n));
        d = complex(zeros(m,n),zeros(m,n));
    else
        PsiTU = zeros(m,n); % equal to PsiT(U) as U=0;
        z = zeros(m,n);
        d = zeros(m,n);
    end
end

last_ii = 0;

%% Main loop
for ii = 1:maxItr
        
    % ================================
    %  Begin Alternating Minimization
    % ----------------
    %   W-subprolem
    % ----------------
    if aTV > 0
        switch TVtype
            case 1;   % anisotropic TV
                Ux = Ux + bx; Uy = Uy + by;      % latest Ux and Uy are already calculated
                Wx = sign(Ux).* max(abs(Ux)-1/beta1,0);
                Wy = sign(Uy).* max(abs(Uy)-1/beta1,0);
            case 2;   % isotropic TV
                [Wx, Wy] = Compute_Wx_Wy(Ux,Uy,bx,by,1/beta1);
            otherwise; 
                error('TVtype must be 1 or 2. Future versions will support more neighborhood types for anisotropic TV.');
        end
    end

    % ----------------
    %   Z-subprolem
    % ----------------
    if aL1 > 0;
        PsiTU = PsiTU + d;
        Z = sign(PsiTU).*max(abs(PsiTU)-1/beta2,0);
    end;
    
    % ----------------
    %   V-subprolem
    % ----------------
    
    V = Numer1 + CircU - g;
    V(picks) = V(picks)/(1+1/beta3);
    
    % ----------------
    %   U-subprolem
    % ----------------
    Uprev = U;
    
    rhs = Ds.*fft2(beta3*(V+g));
    if aTV > 0;
        rhs1 = Compute_rhs_DxtU_DytU(Wx,Wy,bx,by,(aTV*beta1)); % = (aTV*beta1)*(DxtU(Wx-bx)+DytU(Wy-by));
    end
    if (aL1 > 0); 
        if aTV > 0; rhs1 = rhs1 + (aL1*beta2)*Psi(Z-d); else rhs1 = (aL1*beta2)*Psi(Z-d); end
    end

    U = (fft2(rhs1)+rhs)./Denom;    % intermediate U
    if bSymm
        CircU = ifft2(D.*U,'symmetric');
        U = ifft2(U,'symmetric');
    else
        CircU = ifft2(D.*U);
        U = ifft2(U);
    end
    
    if ~bComplex; U=real(U); end
    
    %
    %  End Alternating Minimization
    % ================================

    % -------------------------------------------------
    % check stopping criterion
    %
    relchg = norm(U-Uprev,'fro')/norm(U,'fro');
    if bPrint && mod(ii,10)==0; 
        fprintf('itr=%d relchg=%4.1e', ii, relchg);
        fprintf(' diff_CU=%4.1e', norm(CircU-V,'fro'));
        if exist('uOrg','var'); fprintf(' snr=%4.1f',snr(U)); end
        fprintf(' ||PV-B||=%4.1e', norm(V(picks)-B));
        fprintf('\n');
    end

    if relchg < relchg_tol
        break;
    end
    
    % ------------------------------------------
    % Update quantities of U
    if aTV > 0; [Ux,Uy]=Compute_Ux_Uy(U); end
    if aL1 > 0; PsiTU = PsiT(U); end
    
    % ------------------------------------------
    % ADM update
    %
    if aTV > 0;
        bx = bx + gamma*(Ux - Wx);
        by = by + gamma*(Uy - Wy);
    end
    if (aL1 > 0); d = d + gamma*(PsiTU - z); end
    g = g + gamma*(V - CircU);

    if false && ii-last_ii >= 10 && relchg<1e-2
        last_ii = ii;
        tmp = 1.2;
        Numer1 = Numer1/tmp;
        Denom = tmp*Denom1 + Denom2;
        beta3 = tmp*beta3;
        g = g/tmp;
        fprintf('beta3=%f ',beta3);
    end

end % outer

Out.iter = ii;

%% reverse normalization?
if opts.normalize; U = U/fctr; end

%% cast to real?
if opts.real_sol; U = real(U); end

end

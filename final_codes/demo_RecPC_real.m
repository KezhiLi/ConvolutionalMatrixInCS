close all; clear; clc;

% --- load image ---
F = phantom(128);
%F=imread('lena.bmp');
%F=double(F)/256;
%F=(double(F)-min(min(double(F))))/(max(max(double(F)))-min(min(double(F))));

[m n] = size(F);

% --- # of subsamples ---
nSamples = round(0.15*m*n); % sample rate 15%

% --- generate subsample mask ---
picks = randsample(m*n,nSamples);
picks = sort(picks);
if (picks(1) ~= 1); picks(1) = 1; end  % make sure 1st sample is taken
pick = false(m,n); pick(picks) = true; % generate logical matrix "pick"

% --- OTF (for real sampling) generation ---
% 1. generate a Romberg-type OTF: random phase, unit magnitude
OTF = exp(rand(m,n)*(2*pi*1i));
% 2. for real samples, make OTF conjugate symmetric
OTF = conjugate_symmetrize(OTF); OTF = OTF./abs(OTF);
% 3. generate a corresponding PSF
PSF = otf2psf(OTF,[m n]); if ~isreal(PSF); error('PSF is not real.'); end
% 4. take circulant measurements and take subsamples
CB = ifft2(OTF.*fft2(F),'symmetric'); CB = CB(pick);
%

%
% --- PF parameters ---
% importnat parameters are .maxItr .beta3 .relchg_tol
% for REAL data, set: opts.bsymm = true, opts.bComplex = true
aTV = 1e-8;
opts = [];
    opts.maxItr = 3000; % max # of iterations
    opts.gamma = 1.618; % noiseless choice = 1.618
    opts.beta1 = 100; % TV ADM parameter
    opts.beta2 = 10; % L1 ADM parameter, not used if aL1=0
    opts.beta3 = 100*aTV; % CU = V ADM parameter
    opts.bsymm = true; % use conj_symm ifft2? If set yes, intermediate U are real
    opts.relchg_tol = 1e-5; % stopping tolerance based on relative change
    opts.real_sol = false; % return solution of real type (i.e., U = real(U))
    opts.bPrint = true; %true; % screen display option
    opts.normalize = false; % New parameter/data normalization was added to make
    opts.bComplex = true;

% --- call PC solver ---
aL1 = 0; WT = []; W = []; % turn off L1-wavelets; TV-minimization is kept
[U,Out] = RecPF_Circ_Eq(m,n,aTV,aL1,pick,PSF,CB,2,opts,WT,W,range(F(:)),F); % this uses equality constraints
% [U,Out] = RecPF_Circ(m,n,aTV,aL1,pick,PSF,CB,2,opts,WT,W,range(F(:)),F); % this applies quadratic penalty

% --- show results ---
U = abs(U); % make U from complex to real
figure(1); clf;
subplot(221); imshow(F,[]); title('original');
subplot(222); imshow(pick,[]); title(sprintf('sample mask, %4.1f%%',100*nnz(pick)/numel(pick)));
subplot(223); imshow(U,[]); title(sprintf('recovery, SNR %4.1f',snr(U)));
subplot(224); imshow(F-U,[]); colorbar; title('error');

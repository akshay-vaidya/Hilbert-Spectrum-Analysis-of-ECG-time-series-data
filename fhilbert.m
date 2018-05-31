function [omega,amp,mag,phase,w] = fhilbert(varargin);
%FHILBERT takes a matrix of IMFs or a single column vector
%         and performs the Hilbert-Huang transform.  The 
%         original input signal can be recovered by 
%         the following:
%
%            x = real((sum(amp.*exp(i*phase),2)));
%
%   Usage:
%    unfiltered
%    [omega,amp,mag,phase,w] = fhilbert(imf,t0,tf,[])
%
%    filtered
%    [omega,amp,mag,phase,w] = fhilbert(imf,t0,tf,H)
%
%   Inputs:
%    imf   = IMFs produced from empirical mode decomposition
%    t0    = signal start time 
%    tf    = signal stop time 
%    H     = frequency response (2 col array [f |H|]) from freqz.m
%
%   Outputs:
%    omega = matrix of instantaneous frequencies in cycles/time
%    amp   = matrix of variable amplitudes
%    mag   = time-frequency matrix of amplitudes (Hilbert spectrum)
%    phase = instantaneous phase angle in radians
%    w     = frequency scale for 'mag'
%
% author: Bradley M. Battista
%   University of South Carolina
%   Department of Geological Sciences
%   701 Sumter Street, EWS 617
%   Columbia, SC. 29208
%
% COPYRIGHT: see the associated COPYRIGHT.txt file, and also
% http://software.seg.org/disclaimer2.txt
% This source code may be found online at:
% http://software.seg.org/2007/0003
%

imf = varargin{1};
t0  = varargin{2};
tf  = varargin{3};
t = 1:size(imf,1);

% Hilbert transform of IMFs
idx = 1:size(imf,1);
h = hilbert([imf-repmat(mean(imf),size(imf,1),1);zeros(size(imf))]);
h = h(idx,:);

% variable amplitudes and instantaneous frequencies
amp   = abs(h);
phase = unwrap(angle(h));

% phase average (or not) for smooth spectogram
flen  = 3;
filtr = ones(1,flen)/flen; % phase averaging filter
if flen == 1
    filtr = [0 1 0]; % no phase averaging filter
end
phase = filtfilt(filtr,1,phase);
phase = filtfilt(filtr,1,phase);
omega = abs(diff(phase))/(2*pi*(tf-t0)/(length(t)-1));
omega = [omega(1,:); omega];

% mandatory frequency cuts
n   = 2; % min samples that define a frequency (n>=5)
lof = 1/(length(t)-1);
%hif = min([(length(t)-1)/(n*(tf-t0)) 100000]);
hif = min([(length(t)-1)/(n*(tf-t0)) 48000]);
lo = find(omega < lof);
hi = find(omega > hif);

omega(lo) = lof;
omega(hi) = hif;
amp(lo)   = 0;
amp(hi)   = 0;

% Uncomment if you want signed polarity in the t-f domain
% omega = repmat(sign(sum(imf,2)),1,size(omega,2)).*omega;

% optional frequency cuts
if ~isempty(varargin{4})
    F = varargin{4};
    [om,i,j] = unique(omega);
    H = interp1(F(:,1),F(:,2),om);plot(om,H);
    H = reshape(H(j),size(amp,1),size(amp,2));
    amp = amp.*H;
end

% create a frequency scale and index it
nW       = 1000; % floor(length(t)*hif); for max freq resolution
omegaMin = min(omega(:));
omegaMax = max(omega(:));
widx     = round( (nW-1)*(omega-omegaMin)/(omegaMax-omegaMin) )+1;
w        = linspace(omegaMin,omegaMax,nW);

% output the hilbert spectrum matrix
[b,k,j] = unique(repmat(shiftdim(t),size(widx,2),1));
mag     = sparse(reshape(widx,size(widx,1)*size(widx,2),1),j,...
    reshape(amp,size(amp,1)*size(amp,2),1),nW,length(t));

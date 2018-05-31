function [indmin, indmax, indzer] = extr(x,t);

% [indmin, indmax, indzer] = EXTR(x,t) finds extrema and zero-crossings
%
% inputs : - x : analyzed signal
%          - t (optional) : sampling times, default 1:length(x)
%
% outputs : - indmin = indices of minima
%           - indmax = indices of maxima
%           - indzer = indices of zero-crossings
%
%  originally done by:
%   Rilling, G., Flandrin, P., and Goncalves, P. (2002)
%       http://perso.ens-lyon.fr/patrick.flandrin/emd.html
%  commented by Bradley Matthew Battista
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

if(nargin==1)
  t=1:length(x);
end

m = length(x);
x1=x(1:m-1);
x2=x(2:m);
indzer = find(x1.*x2<0); % find negative values

if any(x == 0) % find existing zeros if any
  iz = find( x==0 );
  indz = [];
  if any(diff(iz)==1) % find nonconsecutive zeros
    zer = x == 0;
    dz = diff([0 zer 0]);
    debz = find(dz == 1); % find peak-side of a zero
    finz = find(dz == -1)-1; % find trough-side of a zero
    indz = round((debz+finz)/2); % avg peaks and troughs to get zero crossings
  else
    indz = iz;
  end
  indzer = sort([indzer indz]); % organize zero crossings
end
  
% take 2nd order derivative of x to find extrema
d = diff(x);
n = length(d);
d1 = d(1:n-1);
d2 = d(2:n);
indmin = find(d1.*d2<0 & d1<0)+1;
indmax = find(d1.*d2<0 & d1>0)+1;

% repeat for extrema not associated with a zero crossing
if any(d==0)
  
  imax = [];
  imin = [];
  
  bad = (d==0);
  dd = diff([0 bad 0]);
  debs = find(dd == 1);
  fins = find(dd == -1);
  if debs(1) == 1
    if length(debs) > 1
      debs = debs(2:end);
      fins = fins(2:end);
    else
      debs = [];
      fins = [];
    end
  end
  if length(debs) > 0
    if fins(end) == m
      if length(debs) > 1
        debs = debs(1:(end-1));
        fins = fins(1:(end-1));

      else
        debs = [];
        fins = [];
      end      
    end
  end
  lc = length(debs);
  if lc > 0
    for k = 1:lc
      if d(debs(k)-1) > 0
        if d(fins(k)) < 0
          imax = [imax round((fins(k)+debs(k))/2)];
        end
      else
        if d(fins(k)) > 0
          imin = [imin round((fins(k)+debs(k))/2)];
        end
      end
    end
  end
  
  if length(imax) > 0
    indmax = sort([indmax imax]);
  end

  if length(imin) > 0
    indmin = sort([indmin imin]);
  end
  
end  

function [indmin,indmax,indzer] = extr2(x,num_sd);
% [indmin,indmax,indzer] = EXTR2(x,num_sd) improves the 
%     function EXTR(x).  The problem is that some extrema 
%     may be missed where strong trends occur, resulting in a 
%     divergence of a signal from the envelope defined by its 
%     initial extrema. This function examines the temporal separation 
%     of extrema and identifies locations in the signal containing 
%     such trends.  The trends are fitted to a 3rd order and removed 
%     before determining extrema in these locations.  The trends are
%     not, however, removed from the data.
% 
%     The trends are fitted to 3rd order to respect the cubic spline
%     fits that are associated with the extrema during the EMD.  Further,
%     the location of areas demanding these fits is determined between any
%     two extrema whose temporal separation is greater than 'num_sd' standard
%     deviations from the overall extrema separation.
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

warning off

if(nargin==1)
  num_sd = 2;
end

t = [1:length(x)];
x = x(:)';

% find normal extrema
[indmin,indmax,indzer] = extr(x);
idx = sort([indmin indmax]);

% find gaps where no extrema were found
% perhaps because of strong trends
if std(diff(idx)) > 0
    idx2 = find(abs(diff(idx)-mean(diff(idx)))>num_sd*std(diff(idx)))';
else
    idx2 = [];
end
if ~isempty(idx2)
    idx2(:,2) = idx2+1;
    
    % loop through pieces and determine additional extrema
    for n = 1:size(idx2,1)
        tbuf = t(idx(idx2(n,1)):idx(idx2(n,2)));
        xbuf = x(idx(idx2(n,1)):idx(idx2(n,2)));
        
        % determine coefficients of 3rd order fit to piece
        p = polyfit(tbuf,xbuf,3);
        
        % remove 3rd order fit and find more extrema
        [indmx,indmn,jnk] = extr([xbuf-polyval(p,tbuf)]);
        
        indmin = [indmin indmn+idx(idx2(n))];
        indmax = [indmax indmx+idx(idx2(n))];
    end
    
    indmin = unique(indmin);
    indmax = unique(indmax);
end

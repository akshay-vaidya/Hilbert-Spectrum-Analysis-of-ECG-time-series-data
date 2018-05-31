function ppout = splaid(pp,breaks)
%SPLAID Insert additional break points into a piecewise polynomial.  This
%       is used to adjust the upper and lower extrema splines to
%       have common break points such that their polynomial coefficients
%       can be averaged, thus producing the polynomial for the mean
%       spline and saving the emd a little time while improving accuracy.
%
%   PPOUT = SPLAID(PP,ADDBREAKS)
%
% Bradley M. Battista- simplified from the Spline Toolbox
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

breaks = sort(breaks(:).'); lb = length(breaks);

b = pp.breaks;
c = pp.coefs;
l = pp.pieces;
k = pp.order;
d = pp.dim;

index0 = find(breaks<b(1)); l0 = length(index0);
% any of these become left-end points for new pieces to the left of the first
% piece, with their coefs all computed from the first piece, i.e., jl(j) = 1.

index2 = find(breaks>b(l+1)); l2 = length(index2);

% now look at the entries of BREAKS in [B(1) .. B(L+1)]. Any of these which
% are not equal to some B(j) become new left-end points, with the coefs
% computed from the relevant piece.
index1 = (l0+1):(lb-l2);
if isempty(index1)
   index = index1;
   jl = [ones(1,l0),l*ones(1,l2)];
else
   pointer = sorted(b(1:l+1),breaks(index1));
   % find any BREAKS(j) not in B.
   % For them, the relevant left-end point is B(POINTER(INDEX)).
   index = find(b(pointer)~=breaks(index1));
   jl = [ones(1,l0),pointer(index),l*ones(1,l2)];
end
ljl = length(jl);
     % If all entries of BREAKS are already in B, then just return the input.
if ljl==0, ppout = pp; return, end

% if there are any BREAKS to the right of B(L+1), then B(L+1) and all but the
% rightmost of these must become left-end points, with coefs computed from
% the last piece, i.e., JL(j) = L  for these, and the rightmost BREAKS
% becomes the new right endpoint of the basic interval.
if l2>0
   tmp = breaks(lb);
   breaks(lb:-1:(lb-l2+1)) = [breaks(lb-1:-1:(lb-l2+1)),b(l+1)];
   b(l+1) = tmp;
end

% These are all the additional left-end points:
addbreaks = breaks([index0,index1(index),index2]);
% Now compute the new coefficients in lockstep:
x = addbreaks - b(jl);
if d>1 % repeat each point D times if necessary
   x = x(ones(d,1),:);
   omd = (1-d:0).'; jl = d*jl(ones(d,1),:)+omd(:,ones(1,ljl));
end
a = c(jl,:); x = x(:);
for ii=k:-1:2
   for i=2:ii
      a(:,i) = x.*a(:,i-1)+a(:,i);
   end
end

% Now, all that's left is to insert the coefficients appropriately.
% First, get the enlarged breaks sequence:
newbreaks = sort([b, addbreaks]);
% This should be of length  L + length(JL)  +  1, requiring
newc = zeros(d*(length(newbreaks)-1),k);
if d>1
   temp = d*sorted(newbreaks,b(1:l));
   newc(temp(ones(d,1),:)+omd(:,ones(1,l)),:) = c;
   temp = d*sorted(newbreaks,addbreaks);
   newc(temp(ones(d,1),:)+omd(:,ones(1,ljl)),:) = a;
else
   newc(sorted(newbreaks,b(1:l)),:) = c;
   newc(sorted(newbreaks,addbreaks),:) = a;
end

ppout = mkpp(newbreaks,newc,d);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pointer = sorted(meshsites, sites)

[ignored,index] = sort([meshsites(:).' sites(:).']);
pointer = find(index>length(meshsites))-[1:length(sites)];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

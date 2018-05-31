function [imf,nbits] = semd(varargin);
%EMD is the empirical mode decomposition.  A signal is reduced to its
%lowest frequency trend by repetitively removing its mean spline.  The
%resulting intrinsic mode functions may then be fed to the hilbert
%transform for analysis.  
%
% imf = semd(x);          for the general application
% imf = semd(x,n);        to limit number of IMFs to 'n'
% [imf,nbits] = semd(...);
%
%   where nbits is the number of sifting iterations performed for each IMF
%
%  originally done by:
%   Rilling, G., Flandrin, P., and Goncalves, P. (2002)
%       http://perso.ens-lyon.fr/patrick.flandrin/emd.html
% Optimized by Bradley M. Battista
%   University of South Carolina
%   Department of Geological Sciences
%   701 Sumter Street, EWS 617
%   Columbia, SC. 29208
%
%  The spline function was modified to reduce the calculation of the mean
%  spline to only one interpolation per sift instead of two.  Some time is 
%  saved but the major effect is a more accurate mean spline.  The improved
%  accuracy is a result of adjusting the upper and lower extr2ema piecewise-
%  polynomials to have the same break points and then averaging the
%  coefficients to create a new polynomial describing the mean spline instead of
%  averaging the interpolated results.  
%
% COPYRIGHT: see the associated COPYRIGHT.txt file, and also
% http://software.seg.org/disclaimer2.txt
% This source code may be found online at:
% http://software.seg.org/2007/0003
% 

warning off;
clipcnt = 1;
x = varargin{1};
t = 1:length(x);


if size(x,1) > 1
    x = x';
end

sd = min(abs(x))*1e-3; % intitial starting condition

% maximum number of IMFs
if nargin >= 2
    nm_imf = varargin{2};
else
    nm_imf = 5;
end

% maximum number of iterations
if nargin == 3
    MAXITERATIONS = varargin{3};
else
    MAXITERATIONS = inf;
end

% maximum number of symmetrized points for interpolations
NBSYM = 2;

lx = length(x);

% maximum number of extr2ema and zero-crossings in residual
ner = lx;
nzr = lx;

r = x;
imf = [];
k = 1;

% iterations counter for extraction of 1 mode
nbit=0;

% total iterations counter
NbIt=0;

while ner > 2
    
    % current mode
    m = r;
    
    % mode at previous iteration
    mp = m;
    
   sx = sd+1; % initial starting condition
    
    % tests if enough extrema to proceed
    test = 0;
    
    [indmin,indmax,indzer] = extr(m);
    if length(indmin) + length(indmax) <= NBSYM
        imf(k,:) = detrend(m,'constant');
        imf(k+1,:) = m - imf(1,:);
        brkflag = 1;
    else
        lm=length(indmin);
        lM=length(indmax);
        nem=lm + lM;
        nzm=length(indzer);
        brkflag = 0;
        
        j=1;
        col = rand(1,3);
        
        % sifting loop
        while ( nbit <= MAXITERATIONS ) & test ~= 1
            ext = abs(nem-nzm); % num extrema-num zero crossings

            [indmin,indmax,indzer] = extr2(m',2);

            if sx <= sd & ...
                    ( ext <= 1 ) | ...
                    ( nem <= 3 )
                test=1;
%                 disp(['convergence achieved in ',num2str(nbit),' sifts'])
            end
            
            % boundary conditions for interpolations
            % extend signal head and tail to reduce edge effects of splines
            if indmax(1) < indmin(1)
                if m(1) > m(indmin(1))
                    lmax = fliplr(indmax(2:min(end,NBSYM+1)));
                    lmin = fliplr(indmin(1:min(end,NBSYM)));
                    lsym = indmax(1);
                else
                    lmax = fliplr(indmax(1:min(end,NBSYM)));
                    lmin = [fliplr(indmin(1:min(end,NBSYM-1))),1];
                    lsym = 1;
                end
            else
                if m(1) < m(indmax(1))
                    lmax = fliplr(indmax(1:min(end,NBSYM)));
                    lmin = fliplr(indmin(2:min(end,NBSYM+1)));
                    lsym = indmin(1);
                else
                    lmax = [fliplr(indmax(1:min(end,NBSYM-1))),1];
                    lmin = fliplr(indmin(1:min(end,NBSYM)));
                    lsym = 1;
                end
            end
            
            if indmax(end) < indmin(end)
                if m(end) < m(indmax(end))
                    rmax = fliplr(indmax(max(end-NBSYM+1,1):end));
                    rmin = fliplr(indmin(max(end-NBSYM,1):end-1));
                    rsym = indmin(end);
                else
                    rmax = [lx,fliplr(indmax(max(end-NBSYM+2,1):end))];
                    rmin = fliplr(indmin(max(end-NBSYM+1,1):end));
                    rsym = lx;
                end
            else
                if m(end) > m(indmin(end))
                    rmax = fliplr(indmax(max(end-NBSYM,1):end-1));
                    rmin = fliplr(indmin(max(end-NBSYM+1,1):end));
                    rsym = indmax(end);
                else
                    rmax = fliplr(indmax(max(end-NBSYM+1,1):end));
                    rmin = [lx,fliplr(indmin(max(end-NBSYM+2,1):end))];
                    rsym = lx;
                end
            end
            
            tlmin = 2*t(lsym)-t(lmin);
            tlmax = 2*t(lsym)-t(lmax);
            trmin = 2*t(rsym)-t(rmin);
            trmax = 2*t(rsym)-t(rmax);
            
            
            % in case symmetrized parts do not extend enough
            try
                if tlmin(1) > t(1) | tlmax(1) > t(1)
                    if lsym == indmax(1)
                        lmax = fliplr(indmax(1:min(end,NBSYM)));
                    else
                        lmin = fliplr(indmin(1:min(end,NBSYM)));
                    end
                    if lsym == 1
                        error('bug')
                    end
                    lsym = 1;
                    tlmin = 2*t(lsym)-t(lmin);
                    tlmax = 2*t(lsym)-t(lmax);
                end
            catch
%                 disp('left sym error')            
            end
            
            try
                if trmin(end) < t(lx) | trmax(end) < t(lx)
                    if rsym == indmax(end)
                        rmax = fliplr(indmax(max(end-NBSYM+1,1):end));
                    else
                        rmin = fliplr(indmin(max(end-NBSYM+1,1):end));
                    end
                    if rsym == lx
                        error('bug')
                    end
                    rsym = lx;
                    trmin = 2*t(rsym)-t(rmin);
                    trmax = 2*t(rsym)-t(rmax);
                end
            catch
%                 disp('right sym error')
            end
            mlmax =m(lmax); 
            mlmin =m(lmin);
            mrmax =m(rmax); 
            mrmin =m(rmin);
            
            % combine extended signal with extrema
            t_mx = [tlmax t(indmax) trmax];
            t_mn = [tlmin t(indmin) trmin];
            
            m_mx = [mlmax m(indmax) mrmax];
            m_mn = [mlmin m(indmin) mrmin];
            
            % ensure no duplicate extrema exist
            [t_mx,idxmx,jnk] = unique(t_mx);
            [t_mn,idxmn,jnk] = unique(t_mn);
            
            % compute coefficients for each spline
            envmx = spline(t_mx,m_mx(idxmx));
            envmn = spline(t_mn,m_mn(idxmn));
            
            % add spline knots such that upper spline and lower spline
            % contain the same points and can be averaged in spline space
            envmn2 = splaid(envmn,envmx.breaks);
            envmx2 = splaid(envmx,envmn.breaks);
            if ~isequal(size(envmx2.coefs),size(envmn2.coefs))
                % if upper and lower polynomials are not of the
                % same order we have to pad the lower order one
                % with zero coefs to make it higher order or matlab will
                % crash.  It's faster to do this than to evaluate the
                % splines and then average them.
                whichpoly = envmn2.order - envmx2.order;
                if whichpoly < 0
                    envmn2.order = envmx2.order;
                    envmn2.coefs = [zeros(envmn2.pieces,1) envmn2.coefs];
                elseif whichpoly > 0
                    envmx2.order = envmn2.order;
                    envmx2.coefs = [zeros(envmx2.pieces,1) envmx2.coefs];
                end
            end
            % generate mean-spline coefficients
            envmean = mkpp(envmx2.breaks,(envmx2.coefs+envmn2.coefs)/2);
            % evaluate to get mean spline as a time series
            envmoy = ppval(envmean,t);
            
%             plot(t,[m;envmoy]);drawnow
            mp = m;
            m = m - envmoy;

            % generate stopping criterion statistic
            sx = norm(mp)/length(mp) - norm(m)/length(m);
            
            % determine if enough extrema to continue sifting
            [indmin,indmax,indzer] = extr(m);
            lm=length(indmin);
            lM=length(indmax);
            nem = lm + lM;
            nzm = length(indzer);
                        
            % update counters
            nbit=nbit+1;
            NbIt=NbIt+1;
            
        end
%         if test ~= 1
%             disp(char({'convergence not met...';...
%                     ['approximating IMF ',num2str(k)]}));
%         else
%             disp(['found IMF ',num2str(k)]);
%         end
        
        imf(k,:) = m;
        nbits(k) = nbit;
        k = k+1;
        r = r - m;
        
        [indmin,indmax,indzer] = extr(r);
        ner = length(indmin) + length(indmax);
        nzr = length(indzer);
        nbit=1;
        
        if k > nm_imf
            imf(k,:) = r;
            break
        end
    end
    if brkflag == 1
        break
    end
end
imf = imf';
nbits = nbits(:)-1;



function imf = emd2(x);
tic;

c = x(:)'; % copy of the input signal (as a row vector)
N = length(x);

%-------------------------------------------------------------------------
% loop to decompose the input signal into successive IMF

imf = []; % Matrix which will contain the successive IMF, and the residue

while (1) % the stop criterion is tested at the end of the loop
   
   %-------------------------------------------------------------------------
   % inner loop to find each imf
   
   h = c; % at the beginning of the sifting process, h is the signal
   SD = 1; % Standard deviation which will be used to stop the sifting process
   
   while SD > 0.3
      % while the standard deviation is higher than 0.3 (typical value)
      
      % find local max/min points
      d = diff(h); % approximate derivative
      maxmin = []; % to store the optima (min and max without distinction so far)
      for i=1:N-2
         if d(i)==0                        % we are on a zero
            maxmin = [maxmin, i];
         elseif sign(d(i))~=sign(d(i+1))   % we are straddling a zero so
            maxmin = [maxmin, i+1];        % define zero as at i+1 (not i)
         end
      end
      
      if size(maxmin,2) < 2 % then it is the residue
         break
      end
      
      % divide maxmin into maxes and mins
      if maxmin(1)>maxmin(2)              % first one is a max not a min
         maxes = maxmin(1:2:length(maxmin));
         mins  = maxmin(2:2:length(maxmin));
      else                                % is the other way around
         maxes = maxmin(2:2:length(maxmin));
         mins  = maxmin(1:2:length(maxmin));
      end
      
      % make endpoints both maxes and mins
      maxes = [1 maxes N];
      mins  = [1 mins  N];
      
      
      %-------------------------------------------------------------------------
      % spline interpolate to get max and min envelopes; form imf
      maxenv = spline(maxes,h(maxes),1:N);
      minenv = spline(mins, h(mins),1:N);
      
      m = (maxenv + minenv)/2; % mean of max and min enveloppes
      prevh = h; % copy of the previous value of h before modifying it
      h = h - m; % substract mean to h
      
      % calculate standard deviation
      eps = 0.0000001; % to avoid zero values
      SD = sum ( ((prevh - h).^2) ./ (prevh.^2 + eps ) );
      
   end
   
   imf = [imf; h]; % store the extracted IMF in the matrix imf
   % if size(maxmin,2)<2, then h is the residue
   
   % stop criterion of the algo.
   if size(maxmin,2) < 2
      break
   end
   
   c = c - h; % substract the extracted IMF from the signal
   
end

% a=1; 
% b=1; 
% c=size(imf); 
% while b<=c(1)  
%    if a==1 
%         k=2; 
%         figure(a); 
%         subplot(5,1,1) 
%         plot(x); 
%         title('INTRINSIC MODE FUNCTIONS'); 
%     while k<=5  
%         if b>c(1) 
%             break; 
%         end 
%         subplot(5,1,k) 
%         plot(imf(b,:)); 
%         k=k+1; 
%         b=b+1; 
%     end 
%    else 
%         k=1; 
%     while k<=5  
%         if b>c(1) 
%             break; 
%         end 
%         figure(a); 
%         subplot(5,1,k) 
%         plot(imf(b,:)); 
%         k=k+1; 
%         b=b+1; 
%     end 
%    end 
%     a=a+1; 
% end
% figure(90)
% plot(imf(1,:))

% % acharya=imf(1,:);
% % dini=imf(2,:);
% % dinesh=norm(acharya);
% % shiva=norm(dini);
% % 
% % love=acos((((acharya)')*(dini))/((dinesh)*(shiva)));
% % % disp(love)
% % breakup = convang(love(:,:),'rad','deg');
% % disp(breakup)
% % uthra = 89:0.01:90;
% % blast = breakup;
% % figure(500)
% % hist(blast,uthra)
% % % %Calculate number of elements in each bin
% % % 
% % % n_elements = histc(blast,uthra);
% % % size(n_elements)
% % % %Calculate the cumulative sum of these elements using cumsum
% % % shock=0;
% % % for i=1:21
% % %     for j=1:1280
% % %         shock=shock+n_elements(i,j);
% % %         
% % %     
% % %     end
% % % end
% % %  for i=1:21
% % %     for j=1:1280
% % %         pdf=n_elements(i,j)/shock;
% % %         end
% % % end       
% % % %Plot the cumulative distribution like a histogram using bar:
% % 
% % % figure(200)
% % % bar(uthra,pdf,'BarWidth',1)
% % % u=length(blast)^2;
% % % figure(800)
% % % [n,x]=hist(blast,uthra);
% % % % a = randn(1,10000);
% % % % [n,x] = hist(a,50);
% % % plot(x,n/u)
% % % hold on
% % % plot(x,normpdf(x,0,1),'r')
% % acharya=imf(2,:);
% % dini=imf(3,:);
% % dinesh=norm(acharya);
% % shiva=norm(dini);
% % 
% % love=acos((((acharya)')*(dini))/((dinesh)*(shiva)));
% % % disp(love)
% % breakup = convang(love(:,:),'rad','deg');
% % disp(breakup)
% % uthra = 89:0.01:90;
% % blast = breakup;
% % figure(600)
% % hist(blast,uthra)
% % Y=abs(fft(imf(1,:)));
% % figure(778)
% % plot(Y)
% % size(Y)
% m=(imf(1,:))
% figure(1256)
% plot(m)
% 
% Fs =0.005; % Sampling frequency 
% t = 0:1/Fs:1200; % Time vector of 1 second
% % x = chirp(t,0,1,Fs/6);
% nfft = 1024; % Length of FFT
% X = fft(x,nfft);
% X = X(1:nfft/2);
% mx= abs(X);
% f = (0:nfft/2-1)*Fs/nfft; 
% % Y=abs(fft(imf(1,:)));
% figure(207); plot(f,mx); title('Power Spectrum of  Signal'); xlabel('Frequency (Hz)'); ylabel('Power'); 


% summ=0
% for i=1:1280
% entropy=summ+((Y(i))*log(Y(i)));
% end
% figure(999)
% plot(entropy)
% figure(98)
% 
% hist(Y)
% for i=0:1:50
% histc(Y)

toc;
warning('off','all');
warning;
return

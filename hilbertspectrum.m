 Octave = exist('OCTAVE_VERSION');
  load('04015m.mat');
  fid = fopen('04015m.info', 'rt');
  fgetl(fid);
  fgetl(fid);
  fgetl(fid);
  [freqint] = sscanf(fgetl(fid),'Sampling frequency: %f Hz  Sampling interval: %f sec');
  interval = freqint(2);
  fgetl(fid);
  if(Octave)
    for i=1:size(val,1)
      R=split(fgetl(fid),char(9));
      signal{i} = R(2,:);
      gain(i) = str2num(R(3,:));
      base(i) = str2num(R(4,:));
      units{i} = R(5,:);
    end
  else
    for i=1:size(val,1)
      [row(i), signal(i), gain(i), base(i), units(i)]=strread(fgetl(fid),'%d%s%f%f%s','delimiter','\t');
    end
  end
  fclose(fid);
  val(val==-32768) = NaN;
  for i=1:size(val,1)
    val(i,:) = (val(i,:)-base(i))/gain(i);
  end
  t = (1:size(val,2))*interval;
  plot(t',val');
  for i=1:length(signal),labels{i}=strcat(signal{i},' (',units{i},')'); end
  legend(labels);
  xlabel('Time (sec)');
p=val';
p=detrend(p);
figure(143)
plot(t',p);
title('(detrend) base line shift');
xlabel('Time (sec)');
ylabel('amplitude(mV)');
imf=semd(p(1:30000));
t=t(1:30000);
tn=length(t);

t1=min(t(1),t(tn));
t2=max(t(1),t(tn));
imfStart = 1;
imfStop = size(imf,2);
imfLast = imfStop+1;

for j=imfStart:imfLast;
    subplot(imfLast-imfStart+1,1,j+1-imfStart);
    if j==imfStart
        z=sum(imf(:,imfStart:imfStop),2);
        plot(t,z);
        yMax = max(z);
        yMin = min(z);
    else
        plot(t,imf(:,j-1));
        yMax = max(imf(:,j-1));
        yMin = min(imf(:,j-1));
    end
    if( abs(yMax-yMin) > 1e-10 )
        axis([t1 t2 yMin yMax]);
    else
        axis([t1 t2 -1 1]);
    end
   s1='c';
   s2=[s1 num2str(j-1)];
   ylabel(s2);
   if( j == 1 )
       title( 'Time Series and IMF Components' );
   end
   if( j < imfLast );
      set(gca,'xticklabel','');
   end;
end;
xlabel('Time');


%%%%%%%%%%%%%%%%%%% Plot the Hilbert Spectrum %%%%%%%%%%%%%%%%%%%

% Perform the Hilbert transform and 
% get instantaneous attributes
[omega,amp,mag,phase,w] = fhilbert(imf,t1,t2,[]);

% Plot the Hilbert Spectrum
[m,n] = size(mag);
h     = full(sum(mag,2)/n);

% Need time,frequency and amplitude to be column vectors
T = repmat(t(:),size(imf,2),1);
O = omega(:);
A = amp(:);

% Sum amp for identical freqs in different IMFs
stoa = sortrows([T O A],[2 1]);
T = stoa(:,1);
O = stoa(:,2);
A = stoa(:,3);
myeps = max(max(abs(T)),max(abs(O)))*eps^(1/3);
ind = [0; ((abs(diff(O)) < myeps) & (abs(diff(T)) < myeps)); 0];
if sum(ind) > 0
  disp(['Summing amplitudes at identical frequencies.']);
  fs = find(ind(1:end-1) == 0 & ind(2:end) == 1);
  fe = find(ind(1:end-1) == 1 & ind(2:end) == 0);
  for n = 1 : length(fs)
    % summing z values
    A(fe(n)) = sum(A(fs(n):fe(n)));
  end
  T = T(~ind(2:end));
  O = O(~ind(2:end));
  A = A(~ind(2:end));
end

figure;
ax1 = subplot(3,1,1);plot(t,sum(imf,2));
ylabel('Amplitude');set(gca,'XTickLabel','','xlim',[min(t) max(t)]);
grid on;box on;

ax2 = subplot(3,1,2);scatter(T,O,10,20*log10(A),'o','filled');
set(ax2,'ydir','normal','xlim',[min(t) max(t)]);

xlabel('Time');ylabel('Frequency');
%set(gca,'yscale','log');
grid on;box on

ax4 = colorbar('horiz');axes(ax4);
xlabel('Power (dB)')
grid on;box on;

ax3 = subplot(3,1,3);plot(h,w);
xlabel('Amp');set(gca,'YTickLabel','');
grid on;box on

set(ax1,'Position',[0.1300 0.5822 0.6770 0.3222]);
set(ax2,'Position',[0.1300 0.2100 0.6770 0.3222]);
set(ax3,'Position',[0.8480 0.2100 0.0770 0.3222]);
set(ax4,'Position',[0.1300 0.1000 0.6770 0.0100]);

% Plot the regenerated signal against the original signal
sig = real((sum(amp.*exp(i*phase),2)));%+imf(:,end);
figure;plot(t,[sum(imf,2) sig]);
xlabel('Time');
legend('\Sigma IMFs','\Sigma a_je^i^\Theta^t');

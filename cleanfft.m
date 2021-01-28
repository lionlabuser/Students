function [sfft,f,bigtaillefft] =cleanfft(s,fs,dessin,col)

% function [sfft,f,bigtaillefft] =cleanfft(s,fs,[dessin,col])
% Computes the Fast Fourier Transform (fft) of signal s (hanning window)
% parameters are the sampling rate fs, and the picture type 'dessin' (if any) 
% dessin is a 2 letters string, the first one setting the frequency axis 
% and the second one the amplitude ('l' linear, 'd' dB).
% col specifies the color - useful when overlaying plots
% Returns the fft (_complex_), the frequency axis and the fft size.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fenetrage
if size(s,2)~=1;
   s = s';
end

s = s.*hanning(length(s));
% Padding avedc zeros pour avoir puissance de 2
bigtaillefft=2^nextpow2(length(s));
% bigtaillefft = 16384;
sfenetre = zeros(1, bigtaillefft);
sfenetre(1:length(s)) = s;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
   col = 'b';
end

bigsfft = fft(sfenetre,bigtaillefft);

%on  reduit la fft de moitie car on a un signal reel
taillefft = bigtaillefft/2 +1;

%Attention: taillefft comprends le 1er terme
sfft = bigsfft(1:taillefft);

%on a jete la moitie, donc
sfft = sfft*2; 

%pour ne pas compter le composant dc *2
sfft(1) = sfft(1)/2;

%de meme pour nyquist qui existe car bigtaille est pair  
% sfft((taillefft-1)/2) = sfft((taillefft-1)/2)/2; 

%on normalise par rapport a la taille
sfft=sfft/length(s);

%transpose
sfft = sfft';

%On construit l'echelle des frequences
f=(0:taillefft-1)*fs/bigtaillefft;
f=f';
  
if nargin > 2
  if dessin == 'ld'
    plot(f,20*log10(abs(sfft)),'color',col,'linewidth',1)
    xlabel('Frequency (Hz, linear)')
    ylabel('Amplitude (dB)')
  elseif dessin == 'll'
    plot(f,abs(sfft),'color',col)
    xlabel('Frequency (Hz, linear)')
    ylabel('Amplitude (linear)')
  elseif dessin == 'dl'   
    semilogx(f,abs(sfft),'color',col)
    xlabel('Frequency (Hz, log)')
    ylabel('Amplitude (linear)')
  elseif dessin == 'dd'   
    semilogx(f,20*log10(abs(sfft)),'color',col)
    xlabel('Frequency (Hz, log)')
    ylabel('Amplitude (dB)')
  end
  zoom on
end

clc
clear all
close all

sampf = 100;
t = [1:1/sampf:301]';
s = 1*sin(2*pi/10*t);

%[specp, specf, spect]=pspectrum(s, sampf,'spectrogram');
%[psp, psf]=pspectrum(s, sampf, 'FrequencyResolution',0.01);
[psp, psf, pst]=pspectrum(s, sampf, 'FrequencyResolution',0.02,'spectrogram','FrequencyLimits',[0 10],'OverlapPercent',50,'Leakage',0.85);

[fH, fA] = fftBasic(s, sampf);


%%
% figure(1)
% plot(t, s)

% figure(2)
% hold on
% plot(psf, psp)
% plot(fH, fA.*fA/2)
% xlim([0, 1])


figure(2)
contourf(pst,psf,psp,'LineStyle','none')
colorbar
ylim([0, 0.5])

figure(11)
stft(t,s,50,0.5)
ylim([0 0.5])


%%
function stft(t, data, windowSec, overlapRatio)
    fs = 1/(t(2) - t(1));% sampling frequency
    
    colormap(jet)
    windowsize = ceil(windowSec*fs)% size of the "viewing window" to see/analyze a chunk/portion of the time series
    window = hanning(windowsize);% creating a HANN window of size="windowsize" to prevent "edg-artifacts"
    nfft = windowsize;% no. of FFT points
    noverlap = windowsize*overlapRatio;% no. of points to overlap between adjacent chunks of the time series to be FFT'd
    [S,F,T]=spectrogram(data,window,noverlap,nfft,fs);% command to plot the spectrogram with the given inputs
    %imagesc(T,F,abs(S))% display the image with scaled colours / utilize the full range of the colormap
    %pbaspect([1 1 1])% similar to "set size square" of GNUPLOT
    contourf(T,F,abs(S),'LineStyle','none')
    T
    
    
    %axis([5 45 0.1 10])
    h=colorbar;
    %caxis([0 4])
    set(gca,'YDir','Normal')% plot-tool
    xlabel('Time (secs)', 'Fontsize',20)
    ylabel('Freq (Hz)', 'Fontsize',20)
    set(gca, 'FontSize',15)    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%  Q3  %%%%%%%%%%%%%%%%%%%%%
clc;clear all ; close all ;
load('Electrodes.mat')
locsX = Electrodes.X;
locsY = Electrodes.Y;
labels =  Electrodes.labels;
fs = 250;
%% Part 1.A
load('NewData1.mat')
EEG_plotter(EEG_Sig,1);
title("EEG signals")
%% 1.C
[F,W,~]=COM2R(EEG_Sig,21);
S = W*EEG_Sig;
%% 1.D.1 :time domain
EEG_plotter(S,2);
title('Extracted sources for the 1st signal');
%%
t = (1/fs)*(0:1:(length(S)-1));
for j = 1:7
    figure(j)
    for i = 1:3
        subplot(3,1,i);
        plot(t,S(3*(j-1)+i,:));
        title("time domain of source"+num2str(3*(j-1)+i))
        xlabel("time(s)"); ylabel("amplitute")
    end
end
%% 1.D.2 frequency domain
figure()
for j =1:7
    figure(j)
    for i = 1:3
        subplot(3,1,i);
        [fft_pw,freq] = pwelch(S(3*(j-1)+i,:),[],[],[],fs);
        plot(freq,fft_pw);
        xlabel('frequency(Hz)');ylabel("Amplitute")
        title("frequency response of source"+num2str(3*(j-1)+i))
    end
end
%% 1.D.3:spatial features
for j=1:6
figure()
    for i = 1:4
        if 4*(j-1)+i<22
            subplot(2,2,i)
            plottopomap(locsX,locsY,labels,F(:,4*(j-1)+i));
            title("spatial plot of source"+num2str(4*(j-1)+i))
        end
    end
end
%%selected vectors
SelectedSources = [ 2 3 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21];
X_den = F(:,SelectedSources)*S(SelectedSources,:);
EEG_plotter(X_den,1);
title('Denoised signal 1');
%%
%% Part 2.A
load('NewData2.mat')
EEG_plotter(EEG_Sig,1);
title("EEG signals")
%%2.C
[F,W,~]=COM2R(EEG_Sig,21);
S = W*EEG_Sig;
%%2.D.1 :time domain
EEG_plotter(S,2);
title('Extracted sources for the 1st signal');
%%%
t = (1/fs)*(0:1:(length(S)-1));
for j = 1:7
    figure(j)
    for i = 1:3
        subplot(3,1,i);
        plot(t,S(3*(j-1)+i,:));
        title("time domain of source"+num2str(3*(j-1)+i))
        xlabel("time(s)"); ylabel("amplitute")
    end
end
%%2.D.2 frequency domain
figure()
for j =1:7
    figure(j)
    for i = 1:3
        subplot(3,1,i);
        [fft_pw,freq] = pwelch(S(3*(j-1)+i,:),[],[],[],fs);
        plot(freq,fft_pw);
        xlabel('frequency(Hz)');ylabel("Amplitute")
        title("frequency response of source"+num2str(3*(j-1)+i))
    end
end
%%2.D.3: spatial features
for j=1:6
figure()
    for i = 1:4
        if 4*(j-1)+i<22
            subplot(2,2,i)
            plottopomap(locsX,locsY,labels,F(:,4*(j-1)+i));
            title("spatial plot of source"+num2str(4*(j-1)+i))
        end
    end
end
%%2.E:selected vectors
SelectedSources = [1,2,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21];
X_den = F(:,SelectedSources)*S(SelectedSources,:);
EEG_plotter(X_den,1);
title('Denoised signal 2');
%% Part 3.A
load('NewData3.mat')
EEG_plotter(EEG_Sig,1);
title("EEG signals")
%% 3.C
[F,W,~]=COM2R(EEG_Sig,21);
S = W*EEG_Sig;
%%3.D.1 :time domain
EEG_plotter(S,2);
title('Extracted sources for the 1st signal');
%% 
t = (1/fs)*(0:1:(length(S)-1));
for j = 1:7
    figure(j)
    for i = 1:3
        subplot(3,1,i);
        plot(t,S(3*(j-1)+i,:));
        title("time domain of source"+num2str(3*(j-1)+i))
        xlabel("time(s)"); ylabel("amplitute")
    end
end
%% 3.D.2 frequency domain
figure()
for j =1:7
    figure(j)
    for i = 1:3
        subplot(3,1,i);
        [fft_pw,freq] = pwelch(S(3*(j-1)+i,:),[],[],[],fs);
        plot(freq,fft_pw);
        xlabel('frequency(Hz)');ylabel("Amplitute")
        title("frequency response of source"+num2str(3*(j-1)+i))
    end
end
%% 3.D.3:spatial features
for j=1:6
figure()
    for i = 1:4
        if 4*(j-1)+i<22
            subplot(2,2,i)
            plottopomap(locsX,locsY,labels,F(:,4*(j-1)+i));
            title("spatial plot of source"+num2str(4*(j-1)+i))
        end
    end
end
%% 3.E selected vectors
SelectedSources = [2, 3, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21];
X_den = F(:,SelectedSources)*S(SelectedSources,:);
EEG_plotter(X_den,1);
title('Denoised signal 3');
%% Part 4.A
load('NewData4.mat')
EEG_plotter(EEG_Sig,1);
title("EEG signals")
%% 4.C
[F,W,~]=COM2R(EEG_Sig,21);
S = W*EEG_Sig;
%% 4.D.1 :time domain
EEG_plotter(S,2);
title('Extracted sources for the 1st signal');
%% 
t = (1/fs)*(0:1:(length(S)-1));
for j = 1:7
    figure(j)
    for i = 1:3
        subplot(3,1,i);
        plot(t,S(3*(j-1)+i,:));
        title("time domain of source"+num2str(3*(j-1)+i))
        xlabel("time(s)"); ylabel("amplitute")
    end
end
%% 4.D.2 frequency domain
figure()
for j =1:7
    figure(j)
    for i = 1:3
        subplot(3,1,i);
        [fft_pw,freq] = pwelch(S(3*(j-1)+i,:),[],[],[],fs);
        plot(freq,fft_pw);
        xlabel('frequency(Hz)');ylabel("Amplitute")
        title("frequency response of source"+num2str(3*(j-1)+i))
    end
end
%% 4.D.3: spatial features
for j=1:5
figure()
    for i = 1:4
        if 5*(j-1)+i<22
            subplot(2,2,i)
            plottopomap(locsX,locsY,labels,F(:,5*(j-1)+i));
            title("spatial plot of source"+num2str(5*(j-1)+i))
        end
    end
end
%% 4.E:selected vectors
SelectedSources = [2,4,5,6,7,9,10,11,12,13,14,15];
X_den = F(:,SelectedSources)*S(SelectedSources,:);
EEG_plotter(X_den,1);
title('Denoised signal 4');
function EEG_plotter(X,mode)
%X : 32-channel data
load('Electrodes') ;
% Plot Data
% Use function disp_eeg(X,offset,feq,ElecName)
offset = max(abs(X(:))) ;
feq = 250 ;
ElecName = {};
if mode ==1
    ElecName = Electrodes.labels ;
else
    for i =1:21
        str = num2str(i);
        ElecName{i} ="s"+str; 
    end
end
disp_eeg(X,offset,feq,ElecName);
end
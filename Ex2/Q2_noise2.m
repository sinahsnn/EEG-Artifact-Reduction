%%%%%%%%%%%%%%%%%%%%  Q2  %%%%%%%%%%%%%%%%%%%%%%%%%
%% part 0 : first initializations
clc; clear all ; close all ;
fs = 250;
load('Ex2.mat')
X_org = X_org;
%%% selecting noise 
X_noise_select = X_noise_2;
%% part A
%%%% SNR = -10
SNR = -10;
P_signal = sum(sum(X_org.^2));
P_noise = sum(sum(X_noise_select.^2));
sigma = sqrt((P_signal/P_noise)*10^(-1*(SNR)/10));
X_SNR10 = zeros(size(X_org,1),size(X_org,2)); 
for i = 1:size(X_org,1)
   X_SNR10(i,:) = X_org(i,:) + sigma*(X_noise_select(i,:));
end
EEG_plotter(X_org,1)
xlabel('time'); ylabel('channels');
title("original EEG signals ")
EEG_plotter(X_SNR10,1);
xlabel('time'); ylabel('channels');
title("Noisy EEG signals with SNR =-10")
%%%% SNR = -20
SNR = -20;
P_signal = sum(sum(X_org.^2));
P_noise = sum(sum(X_noise_select.^2));
sigma = sqrt((P_signal/P_noise)*10^(-1*(SNR)/10));
X_SNR20 = zeros(size(X_org,1),size(X_org,2)); 
for i = 1:size(X_org,1)
   X_SNR20(i,:) = X_org(i,:) + sigma*(X_noise_select(i,:));
end
EEG_plotter(X_SNR20,1);
xlabel('time'); ylabel('channels');
title("Noisy EEG signals with SNR =-20")
%% part B.1 :PCA method
%%%% SNR = -10
source_pca_X_SNR10 = pca(X_SNR10')'*X_SNR10;
EEG_plotter(source_pca_X_SNR10,0);
title('Output with SNR = -10 (using PCA)');
%%%% SNR = -20
source_pca_X_SNR20 = pca(X_SNR20')'*X_SNR20;
EEG_plotter(source_pca_X_SNR20,0);
title('Output with SNR = -20 (using PCA)');
%%
t = (1/fs)*(0:1:(size(source_pca_X_SNR20(1,:),2)-1));
for i = 1:32
    subplot(8,4,i);
    plot(t,source_pca_X_SNR20(i,:));
    title(['channel = s',num2str(i)]);
end
%% Part B.2 :ICA method 
%%%% SNR = -10
[den_ica_X_SNR10,ICA_X_SNR10,~]=COM2R(X_SNR10,32);
source_ICA_X_SNR10 = ICA_X_SNR10*X_SNR10;
EEG_plotter(source_ICA_X_SNR10,0);
title('Output with SNR = -10 (using ICA)');
%%%% SNR = -20
[den_ica_X_SNR20,ICA_X_SNR20,~]=COM2R(X_SNR20,32);
source_ICA_X_SNR20 = ICA_X_SNR20*X_SNR20;
EEG_plotter(source_ICA_X_SNR20,0);
title('Output with SNR = -20 (using ICA)');
%%
t = (1/fs)*(0:1:(size(source_ICA_X_SNR20(1,:),2)-1));
for i = 1:32
    subplot(8,4,i);
    plot(t,source_ICA_X_SNR20(i,:));
    title(['channel = s',num2str(i)]);
end
%% Part D.1:PCA method
chosen_channel_PCA_X_SNR10 = [1,2,3,4];
chosen_channel_PCA_X_SNR20 = [4,17];
pca_X_SNR10 = pca(X_SNR10');
pca_X_SNR20 = pca(X_SNR20');
den_pca_X_SNR10 = inv(pca_X_SNR10);
den_pca_X_SNR20 = inv(pca_X_SNR20);
X_den_pca_X_SNR10 = den_pca_X_SNR10(:,chosen_channel_PCA_X_SNR10)*source_pca_X_SNR10(chosen_channel_PCA_X_SNR10,:);
X_den_pca_X_SNR20 = den_pca_X_SNR20(:,chosen_channel_PCA_X_SNR20)*source_pca_X_SNR20(chosen_channel_PCA_X_SNR20,:);
EEG_plotter(X_den_pca_X_SNR10,1);
title('denoised signal using pca(SNR-10)');
EEG_plotter(X_den_pca_X_SNR20,1);
title('denoised signal using pca(SNR-20)')
%% part D.2:ICA method
chosen_channel_ICA_X_SNR10 = [6,12,25];
chosen_channel_ICA_X_SNR20 = 21;
X_den_ica_X_SNR10 = den_ica_X_SNR10(:,chosen_channel_ICA_X_SNR10)*source_ICA_X_SNR10(chosen_channel_ICA_X_SNR10,:);
X_den_ica_X_SNR20 = den_ica_X_SNR20(:,chosen_channel_ICA_X_SNR20)*source_ICA_X_SNR20(chosen_channel_ICA_X_SNR20,:);
EEG_plotter(X_den_ica_X_SNR10,1);
title('denoised signal using Ica(SNR-10)');
EEG_plotter(X_den_ica_X_SNR20,1);
title('denoised signal using Ica(SNR-20)');
%% E
%%% channel =13,24
t = (1/fs)*(0:1:(size(X_den_ica_X_SNR10(1,:),2)-1));
channel_num =24;
figure(10)
subplot(4,1,1)
plot(t,X_org(channel_num,:))
title("X original in time domain")
subplot(4,1,2)
plot(t,X_SNR10(channel_num,:))
title("X noisy SNR=-10  in time domain")
subplot(4,1,3)
plot(t,X_den_pca_X_SNR10(channel_num,:))
title("X denoised (using pca) in time domain")
subplot(4,1,4)
plot(t,X_den_ica_X_SNR10(channel_num,:))
title("X denoised (using ica) in time domain")
%%% channel =13,24
t = (1/fs)*(0:1:(size(X_den_ica_X_SNR20(1,:),2)-1));
figure(20)
subplot(4,1,1)
plot(t,X_org(channel_num,:))
title("X original in time domain")
subplot(4,1,2)
plot(t,X_SNR20(channel_num,:))
title("X noisy SNR=-20  in time domain")
subplot(4,1,3)
plot(t,X_den_pca_X_SNR20(channel_num,:))
title("X denoised (using pca) in time domain")
subplot(4,1,4)
plot(t,X_den_ica_X_SNR20(channel_num,:))
title("X denoised (using ica) in time domain")
%% PART F 
RRMSE_PCA10 = RRMSE(X_org,X_den_pca_X_SNR10)
RRMSE_ICA10 = RRMSE(X_org,X_den_ica_X_SNR10)
RRMSE_PCA20=RRMSE(X_org,X_den_pca_X_SNR20)
RRMSE_ICA20=RRMSE(X_org,X_den_ica_X_SNR20)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%
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
    for i =1:32
        str = num2str(i);
        ElecName{i} ="s"+str; 
    end
end
disp_eeg(X,offset,feq,ElecName);
end
function y = RRMSE(X_org,X_den)
a = sum(sum((X_org - X_den) .^ 2),2);
b = sum(sum(X_org .^ 2,2));
y = sqrt(a / b);
end

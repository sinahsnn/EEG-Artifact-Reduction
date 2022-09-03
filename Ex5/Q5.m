clc ; clear all ; close all ; 
%% Part A 
fs = 200;           
load('contaminated.mat');load('pure.mat')
disp_eeg(pure,[],200,[],' pure signal in time domain');
disp_eeg(contaminated,[],200,[],'contaminated signal in time domain');
disp_eeg(contaminated - pure,[],200,[],'Predicted EOG signal in time domain');
T_on = zeros(1,length(pure));           
t = 1/fs:1/fs:length(pure)/fs;
idx = find ((t>1.45 & t<3.25) | (t>15.8 & t<16.67) | (t>19.95 & t<23.1));
T_on(1,idx) =1 ;
figure()
plot(t,T_on,'b',"linewidth",1.5); xlabel("Time[s]"); title("Predicted T_{on}");grid on
%% GEVD
C_x = 0;
for i = 1:length(contaminated)
    C_x = C_x + contaminated(:,i)*transpose(contaminated(:,i));
end
C_x = (1/length(contaminated))*C_x;
P_X = 0;
for i = 1:length(contaminated)
    if T_on(i) == 1
        P_X = P_X + contaminated(:,i)*transpose( contaminated(:,i));
    end
end
P_X = (1/sum(T_on)) * P_X;
[V,A] = eig(P_X,C_x);
[~, index] = sort(diag(A), 'descend');
V = V(:, index);
s_all = transpose(V)*contaminated;
disp_eeg(contaminated,[],200,[],'predicted sources using GEVD');
s_selected = s_all;
s_selected(1:2,:) = 0;
X_den_GEVD = mldivide(transpose(V),s_selected);
disp_eeg(X_den_GEVD,[],200,[],'denoised signal using GEVD');
%% DSS 
contaminated_m = zeros(19,length(contaminated));
contaminated_m = contaminated-mean(contaminated, 2);
C = cov(contaminated_m');
[V,D] = eig(C);
D_wh = D^(-1/2);
contaminated_w = D_wh*V'*contaminated_m;
pca_mat =  inv(D_wh*V');
W = zeros(19,19);
for i = 1:19
    w = rand(19,1);
    for j = 1:2000
        r = w'*contaminated_w;
        r = r.*T_on;
        wp = 0;
        for k = 1:length(contaminated)
            wp = contaminated_w(:,k)*r(k) + wp;
        end
        if i>1
            A = W(:,1:i-1);
            wp = (eye(length( A*(A')))-A*(A'))*wp;
        end
        w_new = wp/norm(wp);   
        if immse(w,w_new) < 1e-12
            break;
        end           
        w = w_new;
        
    end
    W(:,i) = w;
end
sources_all_dss = (transpose(W)*contaminated_w);
disp_eeg(sources_all_dss,[],200,[],'predicted sources using DSS');
sources_selected_dss = sources_all_dss;
sources_selected_dss(1:2,:) = 0;
X_den_dss = pca_mat * mldivide(transpose(W),sources_selected_dss);
disp_eeg(X_den_dss,[],200,[],'denoised signal using DSS');
%% part V 
RRMSE_GEVD_method = sqrt(sum(sum((pure - X_den_GEVD) .^ 2))) / sqrt(sum(sum(pure .^ 2)))
RRMSE_DSS_method  = sqrt(sum(sum((pure - X_den_dss ) .^ 2))) / sqrt(sum(sum(pure .^ 2)))



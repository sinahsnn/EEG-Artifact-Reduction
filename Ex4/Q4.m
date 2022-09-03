clc;clear all;close all;
%%
load('Ex1_data.mat');
fs = 100;
C_x = 0;
for i = 1:length(X_org)
    C_x = C_x + X_org(:,i)*transpose(X_org(:,i));
end
C_x = (1/1e4)*C_x;
fs = 100;
t = 1/fs:1/fs:length(X_org)*1/fs; 
len = size(X_org,2);
%% ---------------------------Part A -------------------------------
%% GEVD
T = 4;
period = fs*T; 
P_x = 0;
for i = 1:length(X_org)-period
    P_x = P_x + X_org(:,i)*X_org(:,i+period)';
end
P_x = (1 / (1e4 - period))*P_x;
P_tilda_x = (1/2)*(P_x + transpose(P_x));
[V,D] = eig(P_tilda_x,C_x);
s1_all = V'*X_org;
t = 1/fs:1/fs:length(s1_all(1,:))*1/fs;
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,s1_all(i,:));
    title("Predicted source number"+num2str(i));
end
s1_selected = s1_all;
s1_selected(1:7,:) = 0;
x1_hat = mldivide(V',s1_selected);
t = 1/fs:1/fs:length(x1_hat)*1/fs;
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,x1_hat(i,:),'k');
    hold on
    plot(t,X1(i,:),'r');title("Predicted x1 channel number"+num2str(i));
    legend('Predicted signal','original signal');
end
%% DSS 
X_org_m = zeros(8,10000);
X_org_m = X_org-mean(X_org, 2);
C = cov(X_org_m');
[V,D] = eig(C);
D_wh = D^(-1/2);
X_org_w = D_wh*V'*X_org_m;
pca_mat =  inv(D_wh*V');
W = zeros(8,8);
for i = 1:8
    w = rand(8,1);
    for j = 1:2000
        r = w'*X_org_w;
        n_period = length(r)/period ;
        r_tilda = zeros(1,period);
        for k =1:n_period
            r_tilda = r_tilda + r(1+(k-1)*400:400*k);
        end
        r_tilda = r_tilda/n_period;
        r = repmat(r_tilda,[1 n_period]);
        wp = 0;
        for k = 1:length(X_org)
            wp = X_org_w(:,k)*r(k) + wp;
        end
        if i>1
            A = W(:,1:i-1);
            temp = length( A*(A'));
            
            wp = (eye(8)-A*(A'))*wp;
        end
        w_new = wp/norm(wp);   
        if immse(w,w_new) < 1e-12
            break;
        end           
        w = w_new; 
    end
    W(:,i) = w;
end

s1_all = (W'*X_org_w);
t = 1/fs:1/fs:(1/fs)*length(s1_all(1,:));
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,s1_all(i,:));
    title("predicted source number"+num2str(i));
end
s1_selected = s1_all;
s1_selected(2:8,:) = 0;
x1_hat_dss = pca_mat*inv(W')*s1_selected;
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,x1_hat_dss(i,:),'k');
    hold on
    plot(t,X1(i,:),'r');title("Predicted x1 channel number"+num2str(i));
    legend('Predicted signal','original signal');
end
RRMSE_GEVD = sqrt(sum(sum((X1 - x1_hat) .^ 2))) / sqrt(sum(sum(X1 .^ 2)))
RRMSE_DSS = sqrt(sum(sum((X1 - x1_hat_dss) .^ 2))) / sqrt(sum(sum(X1 .^ 2)))
%% ---------------------------Part B -------------------------------
%% GEVD
T_inter = 3:7;
lambda = zeros(1,length(T_inter));
counter = 1;
for T = 3:7
    period = fs*T; 
    P_x = 0;
    for i = 1:length(X_org)-period
        P_x = P_x + X_org(:,i)*X_org(:,i+period)';
    end
    P_x = (1 / (1e4 - period))*P_x;
    P_tilda_x = (1/2)*(P_x + transpose(P_x));
    [V,D] = eig(P_tilda_x,C_x);
    lambda(counter) = max(diag(D));
    counter = counter + 1 ;
end
[M,Idx] = max(lambda);
%%%
T = T_inter(Idx);
period = fs*T; 
P_x = 0;
for i = 1:length(X_org)-period
    P_x = P_x + X_org(:,i)*X_org(:,i+period)';
end
P_x = (1 / (1e4 - period))*P_x;
P_tilda_x = (1/2)*(P_x + transpose(P_x));
[V,D] = eig(P_tilda_x,C_x);
s1_all = V'*X_org;
s1_selected = s1_all;
t = 1/fs:1/fs:length(s1_all(1,:))*1/fs;
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,s1_all(i,:));
    title("Predicted source number"+num2str(i));
end
s1_selected(1:7,:) = 0;
x1_hat = mldivide(V',s1_selected);
t = 1/fs:1/fs:length(x1_hat)*1/fs;
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,x1_hat(i,:),'k');
    hold on
    plot(t,X1(i,:),'r');
    title("Predicted x1 channel number"+num2str(i));
    legend('Predicted signal','original signal');
end
RRMSE_GEVD = sqrt(sum(sum((X1 - x1_hat) .^ 2))) / sqrt(sum(sum(X1 .^ 2)))
%% DSS
T_inter = 3:7;
X_org_m = zeros(8,10000);
X_org_m = X_org-mean(X_org, 2);
C = cov(X_org_m');
[V,D] = eig(C);
D_wh = D^(-1/2);
X_org_w = D_wh*V'*X_org_m;
pca_mat =  inv(D_wh*V');
%%%
RRMSE_vec = zeros(1,5);
counter = 1 ;
for T = 3:7
    period = fs*T;  
    W = zeros(8,8);
    for i = 1:8
        w = rand(8,1);
        for j = 1:2000
            r = w'*X_org_w;
            n_period = length(r)/period ;
            r_tilda = zeros(1,period);
            for k =1:floor(n_period)
                r_tilda = r_tilda + r(1+(k-1)*period:period*k);
            end
            r_tilda = r_tilda/floor(n_period);
            r = repmat(r_tilda,[1 floor(n_period)]);
            r = [r , r(1:mod(1e4, period))];
            wp = 0;
            for k = 1:length(X_org)
                wp = X_org_w(:,k)*r(k) + wp;
            end
            if i>1
                A = W(:,1:i-1);
                temp = length( A*(A'));

                wp = (eye(8)-A*(A'))*wp;
            end
            w_new = wp/norm(wp);   
            if immse(w,w_new) < 1e-12
                break;
            end           
            w = w_new; 
        end
        W(:,i) = w;
    end
    s1_all = W'*X_org_w;
    s1_selected = s1_all;
    s1_selected(2:8,:) = 0;
    x1_hat_dss = pca_mat*inv(W')*s1_selected;
    RRMSE_vec(counter) = sqrt(sum(sum((X1 - x1_hat_dss) .^ 2))) / sqrt(sum(sum(X1 .^ 2)));
    counter = counter + 1 ;
end
[val,Idx] = min(RRMSE_vec)
X_org_m = zeros(8,10000);
X_org_m = X_org-mean(X_org, 2);
C = cov(X_org_m');
[V,D] = eig(C);
D_wh = D^(-1/2);
X_org_w = D_wh*V'*X_org_m;
pca_mat =  inv(D_wh*V');
period = T_inter(Idx) * fs ; 
W = zeros(8,8);
for i = 1:8
    w = rand(8,1);
    for j = 1:2000
        r = w'*X_org_w;
        n_period = length(r)/period ;
        r_tilda = zeros(1,period);
        for k =1:floor(n_period)
            r_tilda = r_tilda + r(1+(k-1)*period:period*k);
        end
        r_tilda = r_tilda/floor(n_period);
        r = repmat(r_tilda,[1 floor(n_period)]);
        r = [r , r(1:mod(1e4, period))];
        wp = 0;
        for k = 1:length(X_org)
            wp = X_org_w(:,k)*r(k) + wp;
        end
        if i>1
            A = W(:,1:i-1);
            temp = length( A*(A'));

            wp = (eye(8)-A*(A'))*wp;
        end
        w_new = wp/norm(wp);   
        if immse(w,w_new) < 1e-12
            break;
        end           
        w = w_new; 
    end
    W(:,i) = w;
end
s1_all = W'*X_org_w;
figure();
for i = 1:8
    subplot(4,2,i)
    plot(t,s1_all(i,:));
    title("predicted source number"+num2str(i));
end

s1_selected = s1_all;
s1_selected(2:8,:) = 0;
x1_hat_dss = pca_mat*inv(W')*s1_selected;
figure();
for i = 1:8
    subplot(4,2,i)
    plot(t,x1_hat_dss(i,:),'k');
    hold on
    plot(t,X1(i,:),'r');title("Predicted x1 channel number"+num2str(i));
    legend('Predicted signal','original signal');
end
RRMSE_dss = sqrt(sum(sum((X1 - x1_hat_dss) .^ 2))) / sqrt(sum(sum(X1 .^ 2)))
%% ---------------------------Part c -------------------------------
%% GEVD
C_x_tilda = 0;
for i = 1:length(X_org)
    if T1(i) == 1
        C_x_tilda = C_x_tilda + X_org(:,i)*transpose( X_org(:,i));
    end
end
C_x_tilda = (1/sum(T1))*C_x_tilda;
[V,~] = eig(C_x_tilda,C_x);
s2_all = V'*X_org;
s2_selected = s2_all;
t = 1/fs:1/fs:length(s2_all(1,:))*1/fs;
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,s2_all(i,:));
    title("Predicted source number"+num2str(i));
end
s2_selected(1:7,:) = 0;
x2_hat = mldivide(V',s2_selected);
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,x2_hat(i,:),'k');
    hold on
    plot(t,X2(i,:),'r');title("Predicted x2 channel number"+num2str(i));
    legend('Predicted signal','original signal');
end
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,x2_hat(i,:),'k');title("Predicted x2 channel number"+num2str(i));
end
RRMSE_GEVD = sqrt(sum(sum((X2 - x2_hat) .^ 2))) / sqrt(sum(sum(X2 .^ 2)))
%% DSS
X_org_m = zeros(8,10000);
X_org_m = X_org-mean(X_org, 2);
C = cov(X_org_m');
[V,D] = eig(C);
D_wh = D^(-1/2);
X_org_w = D_wh*V'*X_org_m;
pca_mat =  inv(D_wh*V');
W = zeros(8,8);
for i = 1:8
    w = rand(8,1);
    for j = 1:2000
        r = w'*X_org_w;
        r = r.*T1;
        wp = 0;
        for k = 1:length(X_org)
            wp = X_org_w(:,k)*r(k) + wp;
        end
        if i>1
            A = W(:,1:i-1);
            temp = length( A*(A'));
            
            wp = (eye(temp)-A*(A'))*wp;
        end
        w_new = wp/norm(wp);   
        if immse(w,w_new) < 1e-12
            break;
        end           
        w = w_new;
        
    end
    W(:,i) = w;
end
s2_all = (W'*X_org_w);
t = 1/fs:1/fs:(1/fs)*length(s2_all(1,:));
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,s2_all(i,:));
    title("predicted source number"+num2str(i));
end
s2_selected = s2_all;
s2_selected(2:8,:) = 0;
x2_hat_dss =  mldivide(W',s2_selected);
x2_hat_dss = pca_mat*x2_hat_dss;
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,x2_hat_dss(i,:),'k');
    hold on
    plot(t,X2(i,:),'r');title("Predicted x2 channel number"+num2str(i));
    legend('Predicted signal','original signal');
end
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,x2_hat_dss(i,:),'k');title("Predicted x2 channel number"+num2str(i));
end
RRMSE_DSS = sqrt(sum(sum((X2 - x2_hat_dss) .^ 2))) / sqrt(sum(sum(X2 .^ 2)))
%% ---------------------------Part d -------------------------------
%% GEVD
C_x_tilda = 0;
for i = 1:length(X_org)
    if T2(i) == 1
        C_x_tilda = C_x_tilda + X_org(:,i)*transpose( X_org(:,i));
    end
end
C_x_tilda = (1/sum(T2))*C_x_tilda;
[V,~] = eig(C_x_tilda,C_x);
s2_all = V'*X_org;
T1_estimated = my_macd(s2_all(8, :), 2, 0.5, 2, 100);
C_x_tilda = 0;
for i = 1:length(X_org)
    if T1_estimated(i) == 1
        C_x_tilda = C_x_tilda + X_org(:,i)*transpose( X_org(:,i));
    end
end
C_x_tilda = (1/sum(T1_estimated))*C_x_tilda;
[V,~] = eig(C_x_tilda,C_x);
s2_all = V'*X_org;
%%%
s2_selected = s2_all;
t = 1/fs:1/fs:length(s2_all(1,:))*1/fs;
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,s2_all(i,:));
    title("Predicted source number"+num2str(i));
end
s2_selected(1:7,:) = 0;
x2_hat = mldivide(V',s2_selected);
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,x2_hat(i,:),'k');
    hold on
    plot(t,X2(i,:),'r');title("Predicted x2 channel number"+num2str(i));
    legend('Predicted signal','original signal');
end
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,x2_hat(i,:),'k');title("Predicted x2 channel number"+num2str(i));
end
RRMSE_GEVD = sqrt(sum(sum((X2 - x2_hat) .^ 2))) / sqrt(sum(sum(X2 .^ 2)))
%% DSS
X_org_m = zeros(8,10000);
X_org_m = X_org-mean(X_org, 2);
C = cov(X_org_m');
[V,D] = eig(C);
D_wh = D^(-1/2);
X_org_w = D_wh*V'*X_org_m;
pca_mat =  inv(D_wh*V');
W = zeros(8,8);
T1_estimated = T2;
for k=1:2
    for i = 1:8
        w = rand(8,1);
        for j = 1:40
            r = w'*X_org_w;
            r = r.*T1_estimated;
            wp = 0;
            for k = 1:length(X_org)
                wp = X_org_w(:,k)*r(k) + wp;
            end
            if i>1
                A = W(:,1:i-1);
                temp = length( A*(A'));

                wp = (eye(temp)-A*(A'))*wp;
            end
            w_new = wp/norm(wp);   
            if immse(w,w_new) < 1e-3
                break;
            end           
            w = w_new;   
        end
        W(:,i) = w;
    end
    T1_estimated = my_macd(r, 2, 0.5, 2, 100);
end
%%%
s2_all = (W'*X_org_w);
t = 1/fs:1/fs:(1/fs)*length(s2_all(1,:));
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,s2_all(i,:));
    title("predicted source number"+num2str(i));
end
s2_selected = s2_all;
s2_selected(2:8,:) = 0;
x2_hat_dss =  mldivide(W',s2_selected);
x2_hat_dss = pca_mat*x2_hat_dss;
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,x2_hat_dss(i,:),'k');
    hold on
    plot(t,X2(i,:),'r');title("Predicted x2 channel number"+num2str(i));
    legend('Predicted signal','original signal');
end
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,x2_hat_dss(i,:),'k');title("Predicted x2 channel number"+num2str(i));
end
RRMSE_DSS = sqrt(sum(sum((X2 - x2_hat_dss) .^ 2))) / sqrt(sum(sum(X2 .^ 2)))

%% ---------------------------Part e -------------------------------
%% GEVD
len = size(X_org,2);
f_inter =[10,15];
idxf1 = (len/fs) *f_inter(1) ;
idx_symf1 = len - idxf1;
idxf2 = (len/fs) *f_inter(2) ;
idx_symf2 = len - idxf2;
X_fft = fft(X_org(:,1:len-1)')';
idx = [idxf1:idxf2,idx_symf2:idx_symf1];
S_x = 0;
for i = 1:length(idx)
    S_x = S_x + abs(X_fft(:,idx(i)))*transpose( abs(X_fft(:,idx(i))));
end
S_x = S_x/length(idx);
[V,D] = eig(S_x, C_x);
s3_all = V'*X_org;
f =(-len/2:len/2-1)*fs/len;
figure()
for i = 1:8
    subplot(4,2,i)
    plot(f,abs(fftshift(fft(s3_all(i,:)))));
    title("Predicted source number"+num2str(i));
end
s3_selected = s3_all;
s3_selected(1:7,:) = 0;
x3_hat = mldivide(V',s3_selected);
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,x3_hat(i,:),'k');
    hold on
    plot(t,X3(i,:),'r');title("Predicted x2 channel number"+num2str(i));
    legend('Predicted signal','original signal');
end
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,x3_hat(i,:),'k');title("Predicted x2 channel number"+num2str(i));
end
RRMSE_GEVD = sqrt(sum(sum((X3 - x3_hat) .^ 2))) / sqrt(sum(sum(X3 .^ 2)))
%% DSS
f_inter =[10,15];
X_org_m = zeros(8,10000);
X_org_m = X_org-mean(X_org, 2);
C = cov(X_org_m');
[V,D] = eig(C);
D_wh = D^(-1/2);
X_org_w = D_wh*V'*X_org_m;
pca_mat =  inv(D_wh*V');
W = zeros(8,8);
for i = 1:8
    w = rand(8,1);
    for j = 1:2000
        r = w'*X_org_w;
        r = bandpass(r,f_inter, fs);
        wp = 0;
        for k = 1:length(X_org)
            wp = X_org_w(:,k)*r(k) + wp;
        end
        if i>1
            A = W(:,1:i-1);
            temp = length( A*(A'));
            
            wp = (eye(temp)-A*(A'))*wp;
        end
        w_new = wp/norm(wp);   
        if immse(w,w_new) < 1e-12
            break;
        end           
        w = w_new;       
    end
    W(:,i) = w;
end
s3_all = (W'*X_org_w);
f =(-len/2:len/2-1)*fs/len;
figure()
for i = 1:8
    subplot(4,2,i)
    plot(f,abs(fftshift(fft(s3_all(i,:)))));
    title("Predicted source number"+num2str(i));
end
s3_selected = s3_all;
s3_selected(2:8,:) = 0;
x3_hat_dss =  mldivide(W',s3_selected);
x3_hat_dss = pca_mat*x3_hat_dss;
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,x3_hat_dss(i,:),'k');
    hold on
    plot(t,X3(i,:),'r');title("Predicted x2 channel number"+num2str(i));
    legend('Predicted signal','original signal');
end
RRMSE_DSS = sqrt(sum(sum((X3 - x3_hat_dss) .^ 2))) / sqrt(sum(sum(X3 .^ 2)))
%% ---------------------------Part v -------------------------------
%% GEVD
f_inter = [5,25];
n = 300;                    % Filter order
f = [0 (f_inter(1)-1)/(fs/2) f_inter(1)/(fs/2) f_inter(2)/(fs/2) (f_inter(2)+1)/(fs/2) 1];         % Frequency band edges
a = [0  0  1 1 0 0];           % Amplitudes
Num = firpm(n,f,a);
for count = 1:20
   for i = 1:8
        X_org_filtered(i,:) = filtfilt(Num,1,X_org(i,:));
   end
    S_x = 0;
    for i = 1:length(X_org)
        S_x = S_x + X_org_filtered(:,i)*transpose(X_org_filtered(:,i));
    end
    S_x = S_x/length(X_org);
    [V,D] = eig(S_x,C_x);
    s3_all = V'*X_org;
    s3_selected = s3_all;
    s3_selected(1:7,:) = 0;
    [out,f] = dfft_fs(s3_selected(8,:),length(s3_selected(8,:)),fs);
    index = find( f<25 & f>5);
    out = out(index);
    f = f(index);
    temp = zeros(size(out));
    for i = 1:length(out)
        if (abs(out(i))) > 1.5
            temp(i) = 1;
        end
    end
    index = find(temp ==1);
    f1_new = f(index(1));
    f2_new = f(index(end));
    x3_hat = transpose(W)\s3_selected;
    if f1_new == f_inter(1)
        if f2_new == f2
            break;
        end
    end
    f_inter(1) = f1_new;
    f2 = f2_new;
    n = 300;                    % Filter order
    f = [0 (f_inter(1)-1)/(fs/2) f_inter(1)/(fs/2) f2/(fs/2) (f2+1)/(fs/2) 1];         % Frequency band edges
    a = [0  0  1 1 0 0];           % Amplitudes
    Num = firpm(n,f,a);
end
t = 1/fs:1/fs:length(x3_hat)*1/fs;

figure()
for i = 1:8
    subplot(4,2,i)
    [out,f] = dfft_fs(s3_all(i,:),length(s3_all(i,:)),fs);
    plot(f,abs(out));
    title("frequency reponse of estimated source number"+num2str(i));
end
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,x3_hat(i,:),"k");
    hold on
    plot(t,X3(i,:),"r");
    title("pridicted x3 channel number"+num2str(i));
    legend('pridicted signal','original signal');
end

figure()
for i = 1:8
    subplot(4,2,i)
    plot(x3_hat(i,:));
    title("pridicted x3 channel number"+num2str(i));
end
RRMSE_GEVD = sqrt(sum(sum((X3 - x3_hat) .^ 2))) / sqrt(sum(sum(X3 .^ 2)))
%% DSS
X_org_m = zeros(8,10000);
X_org_m = X_org-mean(X_org, 2);
C = cov(X_org_m');
[V,D] = eig(C);
D_wh = D^(-1/2);
X_org_w = D_wh*V'*X_org_m;
pca_mat =  inv(D_wh*V');for i = 1:8
    X_org_m(i,:) = X_org(i,:) - mean(X_org(i,:));
end
f_inter =[5,25];
for count = 1:20
    W = zeros(8,8);
    for i = 1:8
        w = rand(8,1);
        for j = 1:1000
            r = w'*X_org_w;
            r = bandpass(r,f_inter, fs);
            wp = 0;
            for k = 1:length(X_org)
                wp = X_org_w(:,k)*r(k) + wp;
            end
            if i>1
                A = W(:,1:i-1);
                temp = length( A*(A'));

                wp = (eye(temp)-A*(A'))*wp;
            end
            w_new = wp/norm(wp);   
            if immse(w,w_new) < 1e-12
                break;
            end           
            w = w_new;
        end
        W(:,i) = w;
    end
    s3_all = (transpose(W)*X_org_w);
    s3_chosen = s3_all(8,:);
    s3_selected = s3_all;
    s3_selected(2:8,:) = 0;
    x3_hat_dss = inv(W')*s3_selected; %#ok
    [out,f] = dfft_fs(s3_selected(1,:),length(s3_selected(1,:)),fs);
    index = find( f<25 & f>5);
    out = out(index);
    f = f(index);
    temp = zeros(size(out));
    for i = 1:length(out)
        if (abs(out(i))) > 5
            temp(i) = 1;
        end
    end
    index = find(temp ==1);
    f1_new = f(index(1));
    f2_new = f(index(end));
    if f1_new == f_inter(1)
        if f2_new == f_inter(2)
            break;
        end
    end
    f_inter(1) = f1_new;
    f_inter(2) = f2_new;
end   
t = 1/fs:1/fs:(1/fs)*length(s3_all);

figure()
for i = 1:8
    subplot(4,2,i)
    [out,f] = dfft_fs(s3_all(i,:),length(s3_all(i,:)),fs);
    plot(f,abs(out));
    title("frequency reponse of estimated source number"+num2str(i));
end


x3_hat_dss = pca_mat*x3_hat_dss;
%%
figure()
for i = 1:8
    subplot(4,2,i)
    plot(t,x3_hat_dss(i,:),'k');
    hold on
    plot(t,X3(i,:),'r');
    title("pridicted x3 channel number"+num2str(i));

    legend('original signal','estimated signal');
end

RRMSE_DSS = sqrt(sum(sum((X3 - x3_hat_dss) .^ 2))) / sqrt(sum(sum(X3 .^ 2)))

function [out,f] = dfft_fs(signal,n,fs)
fori=fft(signal,n);
fori = fori/fs;
out = fftshift(fori);
temp2 = linspace(-pi,pi,length(out));
f = temp2*fs/(2*pi);
end
function y = my_macd(x, thr1, thr2, high_f, f)
L = thr2*f;
y = x > thr1;
temp = zeros(size(y));
for i = 1:length(x)
    if i - L / 2 <= 0
	buff =[y(1:i-1),y(i+1:i+L/2)];
    elseif i + L / 2 > length(x)
	buff=[y(i-L/2:i-1),y(i+1:end)];
    else
	buff = [y(i-L/2:i-1),y(i+1:i+L/2)];
    end
    if sum(buff) > high_f
        temp(i) = 1;
    else
        temp(i) = 0;
    end
end
y = temp;
end
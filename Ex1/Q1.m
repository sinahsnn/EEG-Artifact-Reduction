%%%%%%%%%%%%%%%%%%%%  Q1  %%%%%%%%%%%%%%%%%%
clc;clear all;close all;
X_load=load("Ex1.mat");
X_load=X_load.X;
%% part a
X =X_load';
for dim =1:3
    X(:,dim) = X(:,dim) - sum(X(:,dim))/(size(X,1)); 
end
figure()
subplot (3,3,[1:6])
scatter3(X(:,1),X(:,2),X(:,3))
title("Visualized Data in 3D domain")
xlabel('X'); ylabel("Y"); zlabel("Z");
subplot(3,3,7)
scatter(X(:,1),X(:,2))
title("Visualized Data in 2D x-y plane")
xlabel('X'); ylabel("Y");
grid on
subplot(3,3,8)
scatter(X(:,1),X(:,3))
title("Visualized Data in 2D x-z plane")
xlabel('X'); ylabel("Z");
grid on
subplot(3,3,9)
scatter(X(:,2),X(:,3))
title("Visualized Data in 2D y-z plane")
xlabel('Y'); ylabel("Z");
grid on
%% part B.1 
Covariance_mat = cov(X)
[V,D] = eig(Covariance_mat)
D_wh = D^(-1/2);
C_final = D_wh*V'*Covariance_mat*V*D_wh;
kernel=-1:0.001:1;
figure()
subplot (3,3,[1:6])
scatter3(X(:,1),X(:,2),X(:,3));
hold on
plot3(kernel*V(1,2),kernel*V(2,2),kernel*V(3,2),'linewidth',2.5)
hold on
plot3(kernel*V(1,1),kernel*V(2,1),kernel*V(3,1),'linewidth',2.5)
hold on
plot3(40*kernel*V(1,3),40*kernel*V(2,3),40*kernel*V(3,3),'linewidth',2.5)
title("Visualized Data with eigen vectors")
xlabel("X"); ylabel("Y"); zlabel("Z");
subplot(3,3,7)
scatter3(X(:,1),X(:,2),X(:,3));
hold on
plot3(kernel*V(1,2),kernel*V(2,2),kernel*V(3,2),'linewidth',2.5)
hold on
plot3(kernel*V(1,1),kernel*V(2,1),kernel*V(3,1),'linewidth',2.5)
hold on
plot3(40*kernel*V(1,3),40*kernel*V(2,3),40*kernel*V(3,3),'linewidth',2.5)
title("Visualized Data with eigen vectors")
xlabel("X"); ylabel("Y"); zlabel("Z");
subplot(3,3,8)
scatter3(X(:,1),X(:,2),X(:,3));
hold on
plot3(kernel*V(1,2),kernel*V(2,2),kernel*V(3,2),'linewidth',2.5)
hold on
plot3(kernel*V(1,1),kernel*V(2,1),kernel*V(3,1),'linewidth',2.5)
hold on
plot3(40*kernel*V(1,3),40*kernel*V(2,3),40*kernel*V(3,3),'linewidth',2.5)
title("Visualized Data with eigen vectors")
xlabel("X"); ylabel("Y"); zlabel("Z");
subplot(3,3,9)
scatter3(X(:,1),X(:,2),X(:,3));
hold on
plot3(kernel*V(1,2),kernel*V(2,2),kernel*V(3,2),'linewidth',2.5)
hold on
plot3(kernel*V(1,1),kernel*V(2,1),kernel*V(3,1),'linewidth',2.5)
hold on
plot3(40*kernel*V(1,3),40*kernel*V(2,3),40*kernel*V(3,3),'linewidth',2.5)
title("Visualized Data with eigen vectors")
xlabel("X"); ylabel("Y"); zlabel("Z");
%% part B.2
Y = transpose(D_wh*V'*X'); 
Covariance_matY =cov(Y)
figure();
scatter3(Y(:,1),Y(:,2),Y(:,3));
title('Visualization of whitened Data');
xlabel('X');ylabel('Y');zlabel('Z'); 
%% Part C.1
COEF_pca = pca(X);
figure()
scatter3(X(:,1),X(:,2),X(:,3));
hold on
plot3(4*kernel*COEF_pca(1,2),kernel*COEF_pca(2,2),kernel*COEF_pca(3,2),'linewidth',2)
hold on
plot3(4*kernel*COEF_pca(1,1),kernel*COEF_pca(2,1),kernel*COEF_pca(3,1),'linewidth',2)
hold on
plot3(40*kernel*COEF_pca(1,3),40*kernel*COEF_pca(2,3),40*kernel*COEF_pca(3,3),'linewidth',2)
title('Visualized Data with eigen vectors (matlab)');
xlabel('X');ylabel('Y');zlabel('Z');
%% Part C.2
figure()
Y_matlab = COEF_pca'*X';
scatter3(Y_matlab(1,:),Y_matlab(2,:),Y_matlab(3,:));
title('Visuallization of whitened Data (using matlab pca)');
xlabel('X');ylabel('Y');zlabel('Z');
covarianceY_matlab =cov(Y_matlab')
%%
[U_svd,S_svd,V_svd] = svd(X);




 


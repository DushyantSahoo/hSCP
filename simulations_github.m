% This script first simulates correlation matrices and then use it for
% generating components from them. This will give an idea of how to use
% hscp_amsgrad code

% number of samples
N = 300;

% number of nodes
p = 50;

% number of components at level 1 and 2
k1 = 15;
k2 = 6;

% density of components
density1 = 0.8;
density2  = 0.5;

% Generate components and subject specific information
W1 = sprandn(p,k1,density1);
W2 = sprandn(k1,k2,density2);
W1 = full(W1);
W2 = full(W2);
lambda1 = cell(N,1);
lambda2 = cell(N,1);
temp_lamb1 = rand(k1,1);
temp_lamb2 = rand(k2,1);
mean1 = zeros(k1,1);
mean2 = zeros(k2,1);
for i=1:N
    lambda1{i} = temp_lamb1 + 0.1*rand(k1,1);
    lambda2{i} = temp_lamb2 + 0.1*rand(k2,1);
    mean1 = mean1 + lambda1{i};
    mean2 = mean2 + lambda2{i};
end
mean2 = mean2/N;

L = W1*W2*diag(mean2)*transpose(W1*W2);
temp_L = diag(1./sqrt(diag(L)));

new_W1 = temp_L*W1;

% corre stores the correlation matrix of each subject
% notice that it is a cell of length N
corre = cell(N,1);

or i =1:N

    temp = new_W1*W2*diag(lambda2{i}) *transpose(new_W1*W2);
    temp = diag(1./sqrt(diag(temp)))*temp*diag(1./sqrt(diag(temp)));
    temp = (temp + temp')/2;
    corre{i}=temp;

end

% Below is the standard setting for hyperparameters which should work in
% almost all the cases
svd_check = 1;
loop = 6000;
eta = 0.1;
beta1 = 0.99;
beta2 = 0.999;
eps = 10^-8;
tol1 = 10^(-4);

% set the number of components you desire and the sparsity
k = [p,10,5];
alpha = [0.1,0.1];
[W_train, lambda_train, error] = hscp_amsgrad(corre,k,alpha,loop,eta,beta1,beta2,eps,tol1,svd_check);



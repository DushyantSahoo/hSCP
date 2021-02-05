function [W, lambda, error] = hscp_amsgrad(A,k,alpha,loop,eta,beta1,beta2,eps,tole,svd_check)

% A contains the input matrices in in a cell format
% for example- A{1}, A{2}, ....

% k contains the number of components in each hierarchy
% for example k = [120,20,4]
% where 120 is the number of nodes in data or the size of input matrix
% 20 is the number of components at level 1
% 4 is the number of components at level 2
% k can also be [120,10] having only one level or [120,40,20,4] having
% three levels

% alpha is the sparsity level of each components at each level, you will
% have to play with it a little bit to get a balance between sprasity in
% the components and good approximation error
% if k = [120,20,4] then alpha can be [1,1] i.e. sparsity for each level

% loop is the number of iterations of gradient descent

% eta, beta1, beta2, eps are the hyperparameters for amsgrad

% tole is the % change in the error before gradient descent stops

% svd_check is used for initializing the algorithm with SVD, details are
% given in the main paper

% There are three outputs-
% 1) W stores components at different level in cells, each cell of W will
% store components at each level
% 2) lambda stores subject specific information in cell, each cell of
% lambda will store subject and level specific information
% 3) error stores % information captured by the decomposition, ideally it
% should be decreasing with the iterations.

% intialize the components and subject specific information using svd
% initialization
[~,subjects] = size(A);
hierarchy = length(k)-1;
svd_error = 0;
for hi=1:hierarchy
    W{hi} = zeros(k(hi),k(hi+1));
    for sub3=1:subjects
        
        if hi==1
            %U = = randn
            %A{sub3} = vineBeta(n,1);
            [U,S,V]=svd(A{sub3});
            W{hi} = W{hi}+ U(:,1:k(hi+1));
            lambda{hi,sub3} = S(1:k(hi+1),1:k(hi+1));
            if svd_check==0
                lambda{hi,sub3} = diag(diag(rand(k(hi+1),k(hi+1)))); 
                W{hi} = rand(k(hi),k(hi+1));
            end
            
        else
            [U1,S1,~]=svd(lambda{(hi-1),sub3});
            %[U,S,V]=svd(A{sub3});

            lambda{hi,sub3} = S1(1:k(hi+1),1:k(hi+1));
            
            W{hi} = W{hi} + U1(:,1:k(hi+1));
            
            if svd_check==0
                lambda{hi,sub3} = diag(diag(rand(k(hi+1),k(hi+1)))); 
                W{hi} = rand(k(hi),k(hi+1));
            end
            
        end
        past_m_lambda{hi,sub3} = zeros(size(lambda{hi,sub3}));
        past_v_lambda{hi,sub3} = zeros(size(lambda{hi,sub3}));
    end
    
    past_v_W{hi} = zeros(size(W{hi}));
    past_m_W{hi} = zeros(size(W{hi}));
    W{hi} = W{hi}/subjects;
    svd_error=svd_error/(subjects*hierarchy);
end

% run amsgrad
for l = 1:loop
    
    Y{1} = W{1};
    for hi=2:hierarchy
       Y{hi} = Y{hi-1}*W{hi}; 
    end
    
    % compute error
    error(l) = 0;
    for hi=1:hierarchy
        for sub = 1:subjects
            error(l) = error(l) + (norm(A{sub}-Y{hi}*lambda{hi,sub}*Y{hi}','fro')/norm(A{sub},'fro'))^2;
            for new = hi:hierarchy
                if hi-new==0
                    X{hi,new,sub} = lambda{hi,sub};
                else
                    temp = W{new};
                    if (new-hi-1 ~=0)
                        for temp11=1:new-hi-1
                            temp=W{new-temp11}*temp;
                        end
                    end
                    X{new,hi,sub}=  temp*lambda{new,sub}*temp';
                end
            end
        end
    end
    
    % check for convergence
    error(l)=error(l)/(hierarchy*subjects);
    if (l >1500)
        temp_tol = (error(l-1000)-error(l))/error(l-1000);
        temp_tol1 = (error(l-100)-error(l))/error(l-100);
        temp_tol2 = (error(l-50)-error(l))/error(l-50);
        if((abs(temp_tol)<tole) && (abs(temp_tol1)<tole) && (abs(temp_tol2)<tole))
            return
        end
        
    end
    for hi=1:hierarchy
        
        % compute gradients for subject specific loadings
        for sub = 1:subjects
            grad_lambda{hi,sub} = diag(diag(-2*(Y{hi}'*(A{sub})*Y{hi}) + 2*(Y{hi}')*Y{hi}*lambda{hi,sub}*(Y{hi}')*Y{hi}));
            
            m_lambda{hi,sub} = beta1*past_m_lambda{hi,sub} + (1-beta1)*grad_lambda{hi,sub};
            
            v_lambda{hi,sub} = beta2*past_v_lambda{hi,sub} + (1-beta2)*grad_lambda{hi,sub}.^2;
            v_lambda{hi,sub} = max(v_lambda{hi,sub},past_v_lambda{hi,sub});
                        
            past_m_lambda{hi,sub} = m_lambda{hi,sub};
            past_v_lambda{hi,sub} = v_lambda{hi,sub};
            
        end
        
        % compute gradients for the components
        grad_W{hi} = zeros(k(hi),k(hi+1));
        for sub4 = 1:subjects
            for i=hi:hierarchy
                if (hi==1)
                    grad_W{hi} = grad_W{hi}+-4*A{sub4}*W{hi}*X{i,hi,sub4} ...
                    + 4*W{hi}*X{i,hi,sub4}*W{hi}'*W{hi}*X{i,hi,sub4};                    
                else
                    grad_W{hi} = grad_W{hi}+-4*Y{hi-1}'*A{sub4}*Y{hi-1}*W{hi}*X{i,hi,sub4} ...
                    + 4*Y{hi-1}'*Y{hi-1}*W{hi}*X{i,hi,sub4}*W{hi}'*Y{hi-1}'*Y{hi-1}*W{hi}*X{i,hi,sub4};
                end        
            end
        end
        
        m_W{hi} = beta1*past_m_W{hi} + (1-beta1)*grad_W{hi};
        
        v_W{hi} = beta2*past_v_W{hi} + (1-beta2)*grad_W{hi}.^2;
        v_W{hi} = max(v_W{hi},past_v_W{hi});
        
        past_m_W{hi} = m_W{hi};
        past_v_W{hi} = v_W{hi};
        
    end
    
    % update gradients

    for hi=1:hierarchy
        for sub = 1:subjects
            lambda{hi,sub} = lambda{hi,sub} - eta*(m_lambda{hi,sub})./(sqrt(v_lambda{hi,sub})+eps);
            temp_lamb = lambda{hi,sub};
            for lambe=1:k(hi+1)
                if (temp_lamb(lambe,lambe)<0)
                    temp_lamb(lambe,lambe) = 0; 
                end
            end
            lambda{hi,sub} = temp_lamb;
        end
        W{hi} = W{hi} - (eta)*(m_W{hi})./(sqrt(v_W{hi})+eps);
        temp = W{hi};
        for n_k1=1:k(hi+1)
            temp(:,n_k1) = proj_L1_Linf(squeeze(temp(:,n_k1)),alpha(hi));
            if hi>1
                temp(temp<0) = 0;
            end
            
        end
        W{hi} = temp;

    end
    
end


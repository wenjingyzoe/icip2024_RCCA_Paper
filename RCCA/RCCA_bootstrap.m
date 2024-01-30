%this funtion give us the optimal alpha base on bootstrap (non-parametric) method 

function [rho_a,A_a,B_a,u_a,S21_a,S11_a,S22_a,weight_a,optimal_alpha,MSE_vec]=RCCA_bootstrap(X,K,m1,m2,alpha_list,u0,Sigma0,maxiter,num_B)

    num_a=length(alpha_list);
    MSE_vec=NaN(num_a,1);
    n=length(X(:,1));

    for  a=1:num_a

       alpha=alpha_list(a);
       [rho_T,A_T,B_T,mu,~,~,~,~,S]=RCCA_alpha_Divergence(X,m1,m2,alpha,u0,Sigma0,K,maxiter);
  


        if anynan(A_T)

         MSE_vec(a)=NaN(1,1);  

        else

         b=1;
         error=0;
         
          while b < num_B+1

            X_bootstrap=datasample(X,n,'Replace',true); %non-parametric bootstrap step
            [~,A_b,B_b,~,~,~,~,~,~]=RCCA_alpha_Divergence(X_bootstrap,m1,m2,alpha,mu,S,K,maxiter);
                   
                   if anynan(A_b)

       
                     error=NaN;
                     b=num_B+1;
                   else

                    error=(compute_A_mse(A_T,A_b).*rho_T+compute_A_mse(B_T,B_b).*rho_T)+error;
                    disp(error);
                    b=b+1;
                   end

           end

          MSE_vec(a)=sum(error)/num_B;  
      
        end

    end

   index_optimal=find(MSE_vec==min(MSE_vec));
   optimal_alpha=alpha_list(index_optimal);


   disp("optimal alpha is");
   disp(optimal_alpha);
   save("bestalpha.mat","a");
   [rho_a,A_a,B_a,u_a,S21_a,S11_a,S22_a,weight_a,~]=RCCA_alpha_Divergence(X,m1,m2,optimal_alpha,u0,Sigma0,K,maxiter);
  
end
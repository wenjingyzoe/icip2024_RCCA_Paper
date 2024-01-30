%############# algorithm1 provides the robust estimated covariance matrix when the appropriate alpha value is given.########
%input: data--X, a--alpha; u0---initial value of u; Sigma0---initial value of Sigma; maxiter--the maximum iteration steps
%output: estimated mean and sample covariance (details in the paper, see Robust Canonical Correlation Estimation Algorithm)
%##########################################################################################################################

function [u_a,Sigma_a,weight]=algorithm1(X,a,u0,Sigma0,threshold,maxiter)


 u_a=u0; 
 Sigma_a=Sigma0;
 error=10;
 n=length(X(:,1));  
 m=length(X(1,:));
 weight=ones(n,1); 
 num_iter=0;
 
%#############if alpha=1, it gives MLE of mean and covariance matrix################
if a==1
     u_a=mean(X,1);
     ua_matrix=repmat(u_a,n,1);
     X_center=X-ua_matrix;
     Sigma_a=(transpose(X_center)*X_center)/n;
     weight=ones(n,1);

else


   while error>threshold

     num_iter=num_iter+1;
     
     %w_odd=weight;
     u_odd=u_a;
     Sigma_odd=Sigma_a;    
    
     for i=1:n

      xi=X(i,:);
      v=(xi-u_odd);
      C=(v/Sigma_odd);
      C=C*v';
      fac_a=-(1-a)/2;
      weight(i)=exp(fac_a*C);

      end
     

      u_a=zeros(1,m);
      weight_vec=zeros(1,n);

      for i=1:n
        weight_vec(i)=weight(i)/sum(weight);
        xi=X(i,:);
        xi=xi.*repmat(weight_vec(i),1,m);
        u_a=xi+u_a;
      end  


     Sigma_numerator=weight(1)*transpose(X(1,:)-u_a)*(X(1,:)-u_a);

     for i=2:n
        numerator=weight(i)*transpose(X(i,:)-u_a)*(X(i,:)-u_a);
        Sigma_numerator=Sigma_numerator+numerator;
     end
  



     Sigma_a=Sigma_numerator/sum(weight);

     if det(Sigma_a)<10e-18
        Sigma_a=Sigma_a+0.01*eye(m);
     end

    
     diff_u=u_a-u_odd;
     error1=norm(diff_u);
     diff_Sigma=abs(Sigma_a-Sigma_odd);
     error2=norm(diff_Sigma,"fro");
     error=(error1+error2)/2;
     
     %disp(error)

     if num_iter>maxiter

         disp("algorithm has not convergent, choose different alpha")
         u_a=NaN(1,1);
         Sigma_a=NaN(1,1);
         weight=NaN(1,1);
         break
         
     end

   end
 
 end

end







% This function consists of Algorithm 1, which, for a given alpha, gives the corresponding robust estimated covariance matrix.
% RCCA_alpha_Divergence uses the SVD to decompose the robust estimation covariance matrix and then find the CCA direction.


function [rho_vec,A,B,mu,S21,S11,S22,weight_vec,S]=RCCA_alpha_Divergence(X,m1,m2,alpha,u0,Sigma0,K,maxiter)
  


  [mu,S,weight_vec]=algorithm1(X,alpha,u0,Sigma0,10e-4,maxiter);

  if anynan(mu)
  
      rho_vec=NaN(K,1);
      A=NaN(m2,K);
      B=NaN(m1,K);
      S21=NaN(m2,m1);
      S11=NaN(m1,m1);
      S22=NaN(m2,m2);

  else

  S11=S(1:m1, 1:m1);
  S22=S(m1+1:m1+m2 , m1+1:m1+m2);
  S12=S(1:m1,  m1+1:m1+m2);

  C12=((S11^(0.5))\S12)/S22^(0.5);
  [V1,P,V2]=svd(real(C12));
  A=(S22^(0.5))\V2(:,1:K);
  B=(S11^(0.5))\V1(:,1:K);
  rho=diag(P);
  rho_vec=rho(1:K);
  S21=S22*A*diag(rho_vec)*transpose(B)*S11;
  end

end
 

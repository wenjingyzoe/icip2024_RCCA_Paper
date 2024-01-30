%Mean square error 

function A_MSE=compute_A_mse(A_T,A_estimated)

 K=length(A_T(1,:));
 A_MSE=zeros(K,1);


 for l=1:K
    al_T=A_T(:,l);
    al_hat=A_estimated(:,l); 
    cos_a=abs(al_T'*al_hat)/(vecnorm(al_hat)*vecnorm(al_T));
    A_MSE(l)=acosd(cos_a);    
 end

end

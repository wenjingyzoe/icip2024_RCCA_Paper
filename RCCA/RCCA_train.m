PCAK=30;


X1PCA=X1_Train*U2(:,1:PCAK);
X2PCA=X2_Train*U2(:,1:PCAK);

m1=size(X1PCA,2);

m2=size(X2PCA,2);

r1=PCAK;
r2=PCAK;

K=min(r1,r2);

XPCA=[X1PCA X2PCA];

%alpha_lowest=0.975;
%alpha_highest=0.985;

%step=(alpha_highest-alpha_lowest)/11;
%alpha_list=alpha_lowest:step:alpha_highest;
%num_B=10;

%[rho_vec,A,B,mu,S21,S11,S22,weight_a,optimal_alpha,MSE_vec]=RCCA_bootsrap(XPCA,K,m1,m2,alpha_list,mean(XPCA),cov(XPCA),100,num_B);

[rho_vec,A,B,mu,S21,S11,S22,weight_vec,~]=RCCA_alpha_Divergence(XPCA,m1,m2,optimal_alpha,mean(XPCA),cov(XPCA),K,200);
  
%optimal_alpha=alpha_list(10);
%optimal_alpha=0.9965;
%figure;
%scatter(find(label_train==-1),weight_vec(label_train==-1));
%hold on;
%scatter(find(label_train==1),weight_vec(label_train==1));

mu1=mu(1:m1);
mu2=mu(m1+1:m1+m2);

%step3: find canonical correlation variables
num_Train=size(XPCA,1);
X_TrainC=XPCA-repmat(mu,num_Train,1);
canonical_var_M=X_TrainC *[B  zeros(m1, K); zeros(m2, K) A ];




z_vecMatrix=canonical_var_augument_std_Nor(canonical_var_M,K,r1,r2);
z1_vecMatrix=z_vecMatrix(:,1:r1);
z2_vecMatrix=z_vecMatrix(:,r1+1:r1+r2);


%step4: find common var and distinctive var

[CommonVar1_MatrixData,DistinctiveVar1_MatrixData,CommonVar2_MatrixData,DistinctiveVar2_MatrixData]=ComputeCommonAndDistinct_VarMatrixData(z_vecMatrix,rho_vec,r1,r2,K);



%step5: Compute determinist matrix Beta1 and Beta2

[Beta1_Matrix, Beta2_Matrix]=ComputeDeterministicMatrix(X_TrainC,z_vecMatrix,m1,m2,r1,r2);



mu=mean(XPCA);
mu1=mu(1:m1);
mu2=mu(m1+1:m1+m2);

%step9:find CommonPattern and Distinctivepattern in each dataset(in PCA space)
CommonPart1=real(CommonVar1_MatrixData*Beta1_Matrix')+repmat(mu1*0.75,num_Train,1);
%muc1=mean(CommonPart1);
%CommonPart1=real(CommonVar1_MatrixData*Beta1_Matrix')+repmat(muc1,num_Train,1);

CommonPart2=real(CommonVar2_MatrixData*Beta2_Matrix')+repmat(mu2*0.75,num_Train,1);
%muc2=mean(CommonPart2);
%CommonPart2=real(CommonVar2_MatrixData*Beta2_Matrix')+repmat(muc2,num_Train,1);

%mud1=mu1-muc1;
%mud2=mu2-muc2;


DistinctivePart1=real(DistinctiveVar1_MatrixData*Beta1_Matrix')+repmat(mu1*0.25,num_Train,1);
DistinctivePart2=real(DistinctiveVar2_MatrixData*Beta2_Matrix')+repmat(mu1*0.25,num_Train,1);






%transfer to orginal space

CommonPart1=(CommonPart1)*U2(:,1:PCAK)';
CommonPart2=(CommonPart2)*U2(:,1:PCAK)';

DistinctivePart1=(DistinctivePart1)*U2(:,1:PCAK)';
DistinctivePart2=(DistinctivePart2)*U2(:,1:PCAK)';


X1rgb_commonpart=cell(num_Train,1);
X2rgb_commonpart=cell(num_Train,1);

X1rgb_distinctpart=cell(num_Train,1);
X2rgb_distinctpart=cell(num_Train,1);






for i=1:num_Train

    xic1=reshape(CommonPart1(i,:),blockSizeR,blockSizeC,[]);
    xic2=reshape(CommonPart2(i,:),blockSizeR,blockSizeC,[]);
    

    xid1=reshape(DistinctivePart1(i,:),blockSizeR,blockSizeC,[]);
    xid2=reshape(DistinctivePart2(i,:),blockSizeR,blockSizeC,[]);



    X1rgb_commonpart{i}=uint8(xic1);
    X2rgb_commonpart{i}=uint8(xic2);

    X1rgb_distinctpart{i}=uint8(xid1);
    X2rgb_distinctpart{i}=uint8(xid2);

end



X1rgb_Whole_CommonPartpatches=reshape(X1rgb_commonpart,TrainnumRows,TrainnumCols);
X2rgb_Whole_CommonPartpatches=reshape(X2rgb_commonpart,TrainnumRows,TrainnumCols);

X1rgb_Whole_DistinctPartpatches=reshape(X1rgb_distinctpart,TrainnumRows,TrainnumCols);
X2rgb_Whole_DistinctPartpatches=reshape(X2rgb_distinctpart,TrainnumRows,TrainnumCols);



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

X1rgb_Whole_CommonPartEstimated=cell2mat(X1rgb_Whole_CommonPartpatches);
X2rgb_Whole_CommonPartEstimated=cell2mat(X2rgb_Whole_CommonPartpatches);


X1rgb_Whole_DistinctPartEstimated=cell2mat(X1rgb_Whole_DistinctPartpatches);
X2rgb_Whole_DistinctPartEstimated=cell2mat(X2rgb_Whole_DistinctPartpatches);





figure;
imshow(X1rgb_Whole_CommonPartEstimated);
set(gcf, 'Position',  [10, 10, 500, 400])

figure;
imshow(X2rgb_Whole_CommonPartEstimated);
%title('RCCA ');
set(gcf, 'Position',  [10, 10, 500, 400])



figure;
imshow(X1rgb_Whole_DistinctPartEstimated);
%title('RCCA ');
set(gcf, 'Position',  [10, 10, 500, 400])



figure;
imshow(X2rgb_Whole_DistinctPartEstimated);
%title('RCCA ');
set(gcf, 'Position',  [10, 10, 500, 400])


meanOfX2Train=mean(X2_Train);

X2_TrainCenter=X2_Train-repmat(meanOfX2Train,num_Train,1);

X2rgb_TrainCenter=cell(num_Train,1);

for i=1:num_Train

    xi=reshape(X2_TrainCenter(i,:),blockSizeR,blockSizeC,[]);
  
    X2rgb_TrainCenter{i}=uint8(xi);

end



X2rgb_TrainCenter=reshape(X2rgb_TrainCenter,TrainnumRows,TrainnumCols);
X2rgb_TrainCenter=cell2mat(X2rgb_TrainCenter);

figure;
imshow(X2rgb_TrainCenter);
set(gcf, 'Position',  [10, 10, 500, 400])



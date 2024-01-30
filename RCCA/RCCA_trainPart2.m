%denosing Xd
%%%%%%%%%%%%%Estimated Xd=abs(DistinctiveParttern2-DistinctiveParttern1);

DistinctivePatternDiff=abs(DistinctivePattern2-DistinctivePattern1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X2rgbDouble_DistinctPattern=cell(num_Total,1);

for i=1:num_Total
    xi=reshape(DistinctivePatternDiff(i,:),blockSizeR,blockSizeC,[]);
    X2rgbDouble_DistinctPattern{i}=xi;

end

X2rgbDouble_DistinctPattern=reshape(X2rgbDouble_DistinctPattern,TotalnumRows,TotlanumCols);
X2rgbDouble_DistinctPattern=cell2mat(X2rgbDouble_DistinctPattern);

X2rgbDoubleVec_DistinctPattern=reshape(X2rgbDouble_DistinctPattern,[],3);
meanVec_DistinctPattern=mean(X2rgbDoubleVec_DistinctPattern);

%sigma matrix
X2CenteredrgbDoubleVec_DistinctPattern=X2rgbDoubleVec_DistinctPattern-repmat(meanVec_DistinctPattern,Rows*Columns,1);

C=(X2CenteredrgbDoubleVec_DistinctPattern'*X2CenteredrgbDoubleVec_DistinctPattern)/(Rows*Columns);

[~,Dc,Tc]=svd(C);


X2T=abs(X2CenteredrgbDoubleVec_DistinctPattern*Tc(:,1));
mean(X2T);
figure;
scatter(1:length(X2T),X2T);

index_reomve=find(X2T<100);



X2rgbDoubleVec_DistinctPattern(index_reomve,:)=zeros(length(index_reomve),3);
X2removedNoise=zeros(Rows,Columns,3,'uint8');

for i=1:3
    x2D=X2rgbDoubleVec_DistinctPattern(:,i);
    X2removedNoise(:,:,i)= reshape(x2D,Rows,Columns);
end
X2removedNoise=rgb2gray(X2removedNoise);
%X2removedNoise=imbinarize(X2removedNoise,"adaptive");
figure;
imshow(X2removedNoise)
title('Changed Map / RCCA estimated');

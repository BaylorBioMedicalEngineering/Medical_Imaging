I1=[ 0 0 0 0 0 0
    0 0 1 0 0 0
    0 0 1 1 0 0
    0 0 1 1 0 0
    0 0 0 0 0 0
    0 0 0 0 0 0];
I2=[0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0
    0 0 1 1 0 0 0 0 0 0 0 0
    0 0 1 1 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 1 1 0 0 0 0
    0 0 0 0 0 1 1 1 1 0 0 0
    0 0 0 0 0 1 1 1 1 1 0 0
    0 0 0 0 0 1 1 1 1 1 0 0
    0 0 0 0 0 1 1 4 4 1 0 0
    0 0 0 0 0 1 1 4 4 1 0 0
    0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0];
I=I2;

%set parameters, note these huristics were picked for small A
[rows,cols]=size(I);
n=180;
NperTheta=rows+cols;
thetaMin=0;
deltaTheta=pi/n;
thetaMax=thetaMin+(n-1)*deltaTheta;
dist=.5;

%find path
[A,b]=MakeCT(I,thetaMin,thetaMax,deltaTheta,NperTheta,dist);


%use ART



[rows_A,cols_A]=size(A);
%size of the figure window - basically how many iterations
figR=3;
figC=3;
iterations=figR*figC-1;
x0=rand(rows*cols,1);%rows*cols=cols_A
Ivec=reshape(I,rows*cols,1);

% stuff to change
%order of rows
ind1=1:rows_A;
ind2=ind1;
off=79;
for i=2:rows_A
    ind2(i)=mod(ind2(i-1)-1+off,rows_A)+1;
end

order=ind1;
relax=.25;
forget=1;

% ART
[x,error]=ART(A,b,x0,iterations,order,relax,forget);



%show results
True_error=zeros(1,iterations+1);
True_error(1,1)=norm(Ivec-x(:,1))/(rows*cols);
figure, subplot(figR, figC,1);
imshow(mat2gray(reshape(x(:,1),rows,cols)));
title('Initial Iterate');
for iter=1:iterations
    subplot(figR, figC,iter+1);
    imshow(mat2gray(reshape(x(:,iter+1),rows,cols)));
    title(sprintf('Iteration %d',iter));
    True_error(1,iter+1)=norm(Ivec-x(:,iter+1))/(rows*cols);
end
IhatART=reshape(x(:,iterations+1),rows,cols);
%figure, imshow(mat2gray(IhatART));

%BIP

%stuff to adjust
blocks=10;
ind_bip=zeros(blocks,ceil(rows_A/blocks));
for block=1:blocks
    ind_bip(block,:)=block:blocks:rows_A;
end
biprelax=1;
bipforget=.98;
bipw=ones(size(ind_bip))./(rows_A/10);

[xbip,errorbip]=BIP(A,b,x0,iterations,ind_bip,biprelax,bipforget,bipw);
%show results
True_errorbip=zeros(1,iterations+1);
True_errorbip(1,1)=norm(Ivec-xbip(:,1))/(rows*cols);
figure, subplot(figR, figC,1);
imshow(mat2gray(reshape(xbip(:,1),rows,cols)));
title('Initial Iterate');
for iter=1:iterations
    subplot(figR, figC,iter+1);
    imshow(mat2gray(reshape(xbip(:,iter+1),rows,cols)));
    title(sprintf('Iteration %d',iter));
    True_errorbip(1,iter+1)=norm(Ivec-xbip(:,iter+1))/(rows*cols);
end
IhatBIP=reshape(xbip(:,iterations+1),rows,cols);


%SAP

%stuff to adjust
strings=10;
ind_sap=zeros(strings,ceil(rows_A/strings));
for string=1:strings
    ind_sap(string,:)=string:strings:rows_A;
end
saprelax=.25;
sapforget=1;
sapw=ones(strings,1)/strings;

[xsap,errorsap]=SAP(A,b,x0,iterations,ind_sap,saprelax,sapforget,sapw);
%show results
True_errorsap=zeros(1,iterations+1);
True_errorsap(1,1)=norm(Ivec-xsap(:,1))/(rows*cols);
figure, subplot(figR, figC,1);
imshow(mat2gray(reshape(xsap(:,1),rows,cols)));
title('Initial Iterate');
for iter=1:iterations
    subplot(figR, figC,iter+1);
    imshow(mat2gray(reshape(xsap(:,iter+1),rows,cols)));
    title(sprintf('Iteration %d',iter));
    True_errorsap(1,iter+1)=norm(Ivec-xsap(:,iter+1))/(rows*cols);
end
IhatSAP=reshape(xbip(:,iterations+1),rows,cols);



%error plots
figure, 
plot(1:iterations+1,error,'b-',...
     1:iterations+1,True_error,'g-',...
     1:iterations+1,errorbip,'r-.',...
     1:iterations+1,True_errorbip,'m-.',...
     1:iterations+1,errorsap,'c:',...
     1:iterations+1,True_errorsap,'y:');
title('Error Comparison');
legend('ART Residuals','ART True Error',...
       'BIP Residuals','BIP True Error',...
       'SAP Residuals','SAP True Error');
xlabel('Iteration');
ylabel('Error');

fprintf('Error_{ART}-Error_{BIP}=%f\n',min(True_error)-min(True_errorbip));
fprintf('Error_{ART}-Error_{SAP}=%f\n',min(True_error)-min(True_errorsap));

%I=[1 3
%   2 4];
n=100;
I = phantom('Modified Shepp-Logan', n);
figure, subplot(2,2,1), imshow(I), title('original');

theta=0:179;

%use backprojection
sinogram=radon(I,theta);

x0=zeros(n*n,1);
iterations=10;
relax=.001;
forget=1;
tstart = tic;
Iart=ARTsino(sinogram,theta,x0,iterations,relax,forget);
tsimple = toc(tstart);
Iart=mat2gray(reshape(Iart,n,n));
subplot(2,2,2), imshow(Iart), title('ART');

tstart = tic;
Ihatbp2=FBP2(sinogram,theta,'ram-lak');
tlight = toc(tstart);
subplot(2,2,3), imshow(Ihatbp2), title('light fbp');

tstart = tic;
Ir=iradon(sinogram,theta);
tradon = toc(tstart);
subplot(2,2,4), imshow(Ir), title('MatLab');
disp(tsimple)
disp(tlight)
disp(tradon)
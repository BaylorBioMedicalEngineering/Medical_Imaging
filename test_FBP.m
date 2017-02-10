%I=[1 3
%   2 4];
I = phantom('Modified Shepp-Logan', 100);
figure, subplot(2,2,1), imshow(I), title('original');

theta=0:179;

%use backprojection
sinogram=radon(I,theta);

tstart = tic;
Ihatbp=FBP(sinogram,theta);
tsimple = toc(tstart);
subplot(2,2,2), imshow(Ihatbp), title('easy fbp');

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
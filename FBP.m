function I=FBP(sinogram,theta)
    [proj,angles]=size(sinogram);
    n=2*floor(proj/(2*sqrt(2)));
    pad=2^(ceil(log2(proj))+1)-proj;
    step=2/(proj+pad);
    if pad>0
        sino=[sinogram
              zeros(pad,angles)];
    else
        sino=sinogram;
    end
    if min(size(theta))==1
        if angles==max(size(theta))
            I=zeros(n,n);
            %make filter - ramp
            filter=[0:step:1 1-step:-step:step]';
            %do fbp for each angle
            for angle=1:angles
                % filter
                fsino=real(ifft(fft(sino(:,angle)).*filter));
                filtsino=fsino(1:proj);
                % backproject
                T=imrotate(filtsino*ones(1,2*n),theta(angle)+90,'bilinear');
                %sum up just part in image
                [tr,tc]=size(T);
                startr=ceil((tr-n)/2)+1;
                startc=ceil((tc-n)/2)+1;
                Temp=T(startr:startr+n-1,startc:startc+n-1);
                I=I+Temp;
            end
            I=I*pi/(2*angles);
        else
            I=[];
        end
    end
end
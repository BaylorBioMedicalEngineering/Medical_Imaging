function I=FBP2(sinogram,theta,filter)
    [proj,angles]=size(sinogram);
    convrad=pi/180;
    n=2*floor(proj/(2*sqrt(2)));
    ctr=ceil(n/2);
    midproj=ceil(proj/2);
    pad=2^(ceil(log2(proj))+1)-proj;
    step=2/(proj+pad);
    lin=0:step:1;
    w=pi*lin;
    Hcos=cos(w/2);
    Hsin=sin(w);
    switch filter
        case 'ram-lak'
            H=lin;
        case 'shep-logan'
            H=Hsin/pi;
        case 'cos'
            H=lin.*Hcos;
        case 'hann'
            H=lin.*(1+Hcos)/2;
        case 'hamming'
            H=lin.*(.54+.46.*Hcos);
        otherwise
            H=lin;
    end
    H=[H';H(end-1:-1:2)'];
    
    if pad>0
        sino=[sinogram
              zeros(pad,angles)];
    else
        sino=sinogram;
    end
    if min(size(theta))==1
        if angles==max(size(theta))
            I=zeros(n,n);
            %do fbp for each angle
            for angle=1:angles
                % filter
                fsino=real(ifft(fft(sino(:,angle)).*H));
                fsino(proj+1:end)=[];
                % backproject
                cth=cos(convrad*theta(angle));
                sth=sin(convrad*theta(angle));
                for row=1:n
                    sy=sth*(n+1-row-ctr);
                    for col=1:n
                        pt=cth*(col-ctr)+sy+midproj;
                        low=floor(pt);
                        high=low+1;
                        if((low>0)&(high<proj+1))
                            Wl=high-pt;
                            Wh=pt-low;
                            I(row,col)=I(row,col)+Wl*fsino(low)+Wh*fsino(high);
                        end
                    end
                end
            end
            I=I*(pi/(2*angles));
        else
            I=[];
        end
    end
end
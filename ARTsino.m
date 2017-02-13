function [x]=ARTsino(sino,theta,x0,iterations,relax,forget)
    cols_A=length(x0);
    rows=floor(sqrt(cols_A));
    cols=cols_A/rows;
    midrow=rows/2;
    midcol=cols/2;
    [projs,angles]=size(sino);
    midproj=projs/2;
    rows_A = projs*angles;
    theta=theta.*(pi/180);
    x=zeros(cols_A,1);
    xnext=x0;
    for iter=1:iterations
        relax=relax*forget;
        for angle=1:angles
            c=cos(theta(angle));
            s=sin(theta(angle));
            for proj=1:projs
                At=zeros(rows,cols);
                for row=1:rows
                    for col=1:cols
                        if( abs((col-midcol)*c+(rows+1-row-midrow)*s+midproj-proj) <0.5)
                            At(row,col)=1;
                        end
                    end
                end
                Ar=reshape(At,1,cols_A);
                xlast=xnext;
                normAr=norm(Ar)^2;
                if normAr>0
                    residual=sino(proj,angle)-Ar*xlast;
                    xnext=xlast+relax*(residual*Ar')/normAr;
                else
                    xnext=xlast;
                end
            end
        end
        x=xnext;
    end
end
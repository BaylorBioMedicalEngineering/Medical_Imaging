function [x,error]=BIP(A,b,x0,iterations,order,relax,forget,w)
    [rows_A,cols_A]=size(A);
    x=zeros(cols_A,iterations+1);
    error=zeros(1,iterations+1);
    x(:,1)=x0;
    [blocks,projections]=size(order);
    error(1,1)=norm(A*x(:,1)-b)/cols_A;
    for iter=1:iterations
        relax=relax*forget;
        xnext=x(:,iter);
        for i=1:blocks
            xlast=xnext;
            xnext=zeros(size(xnext));
            for bi=1:projections
                normAi=norm(A(order(i,bi),:))^2;
                if normAi>0
                    residual=b(order(i,bi))-A(order(i,bi),:)*xlast;
                    xnext=xnext+w(i,bi)*(residual*A(order(i,bi),:)')/normAi;
                else
                    xnext=xnext;
                end
            end
            xnext=xlast+relax*xnext;
        end
        x(:,iter+1)=xnext;
        error(1,iter+1)=norm(A*x(:,iter+1)-b)/cols_A;
    end
end
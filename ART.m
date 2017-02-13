function [x,error]=ART(A,b,x0,iterations,order,relax,forget)
    [rows_A,cols_A]=size(A);
    
    x=zeros(cols_A,iterations+1);
    error=zeros(1,iterations+1);
    x(:,1)=x0;
    error(1,1)=norm(A*x(:,1)-b)/cols_A;
    for iter=1:iterations
        xnext=x(:,iter);
        relax=relax*forget;
        for i=1:rows_A
            xlast=xnext;
            normAi=norm(A(order(i),:))^2;
            if normAi>0
                residual=b(order(i))-A(order(i),:)*xlast;
                xnext=xlast+relax*(residual*A(order(i),:)')/normAi;
            else
                xnext=xlast;
            end
        end
        x(:,iter+1)=xnext;
        error(1,iter+1)=norm(A*x(:,iter+1)-b)/cols_A;
    end
end
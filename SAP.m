function [x,error]=SAP(A,b,x0,iterations,order,relax,forget,w)
    [rows_A,cols_A]=size(A);
    x=zeros(cols_A,iterations+1);
    error=zeros(1,iterations+1);
    x(:,1)=x0;
    [strings,projections]=size(order);
    error(1,1)=norm(A*x(:,1)-b)/cols_A;
    for iter=1:iterations
        relax=relax*forget;
        xlast=x(:,iter);
        xnext=zeros(size(xlast));
        for string=1:strings
            xstring=xlast;
            for projection=1:projections
                normAi=norm(A(order(string,projection),:))^2;
                if normAi>0
                    residual=b(order(string,projection))-A(order(string,projection),:)*xstring;
                    xstring=xstring+relax*(residual*A(order(string,projection),:)')/normAi;
                else
                    xstring=xstring;
                end
            end
            xnext=xnext+w(string)*xstring;
        end
        x(:,iter+1)=xnext;
        error(1,iter+1)=norm(A*x(:,iter+1)-b)/cols_A;
    end
end
% 二分法

function [c, err, yc] = bisect(f, a, b, tol, f_val)

% 找到函数f=f_val在区间[a, b]上的根
% f: 函数句柄
% a, b: 包含根的区间
% tol: 精度要求

err = abs(b - a); % 初始误差

while err > tol

    c = (a + b) / 2;
    yc = feval(f, c);

    if yc == f_val
        break; % 找到精确根
    elseif yc > f_val
        a = c;
    else
        b = c;
    end
    err = abs(b - a);
end

c = b;
err = abs(b - a);
yc = feval(f,c);

end
%{
function x = bisect(f, a, b, tol, f_val)

    while 1
    
        t_0 = (a + b)/2;
        f_0 = feval(f, t_0);
        
        if f_0 < f_val
        
            b = t_0;
        else 
            a = t_0;
        end
    
        if abs(b-a) < tol
            x = (a+b)/2;
            break;
        end
    
    end
end
%}
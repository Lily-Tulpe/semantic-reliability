% 计算函数一阶导数
function fun_diff = calculate_fun_diff(fun, x)
    
    % fun是函数句柄，x是自变量
    h = 1e-6;
    len = length(x);
    fun_diff = zeros(len,1);
    
    for i=1:len
        temp=zeros(size(x));
        temp(i)=h;
        fun_diff(i)=(fun(x+temp)-fun(x))/h;
    end

end
% 通信错误率函数
function err_r = calculate_err_r(t_r, C, V, d, T_s)
    
    % tr为通信时间，C为信道容量，V为信道色散，d为数据量，Ts为一个符号持续时间
    
    L_r = t_r/T_s;
    err_r = qfunc(sqrt(L_r/V).*(C - d./L_r)*log(2));
   
end
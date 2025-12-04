% 感知错误率函数
function err_s = calculate_err_s(t_s, P_s, h_s, sigma2, delta, T_s)
    
    % ts为感知时间，Ps为感知功率，hs为感知信道系数，delta为虚警率，sigma2为噪声功率，Ts为单个符号持续时间
    
    L_s = t_s / T_s;
    err_s = qfunc((P_s * L_s * h_s^2 - (sigma2) * (-log(delta))) ./ (sqrt(2 * P_s * L_s * (sigma2) * h_s^2)));
 
end
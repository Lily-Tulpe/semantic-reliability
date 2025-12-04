% 感知错误率函数，h_s未知

function err_s = calculate_err_s_unh_s(t_s, P_s, gain, sigma2, delta, T_s)
    
    % ts为感知时间，Ps为感知功率，hs为感知信道系数，delta为虚警率，sigma2为噪声功率，Ts为单个符号持续时间
    
    L_s = t_s / T_s;
    gamma_s = P_s*gain/sigma2;

    a = sqrt(2*L_s*gamma_s);
    b = sqrt(2*(-log(delta)));
    m = 1;

    P_D = marcumq(a, b, m);
    err_s = 1 - P_D;

end
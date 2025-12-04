% 计算感知时间阈值的函数
function ts_th = calculate_ts_th(err_s_th, P_s, h_s, sigma2, delta, T_s)

    % err_s_th为感知错误率阈值, Ps为感知功率, hs为感知信道系数
    % sigma2为噪声功率, delta为虚警率, Ts为单个符号持续时间

    y = qfuncinv(err_s_th);
    Ls_sqrt = (y*sqrt(2*P_s*sigma2*(h_s.^2))+sqrt((y.^2)*2*P_s*sigma2*(h_s.^2)+4*P_s*(h_s.^2)*sigma2*(-log(delta))))./(2*P_s*(h_s.^2));
    L_s = ceil(Ls_sqrt.^2);
    ts_th = L_s*T_s;

end
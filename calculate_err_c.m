% 计算错误率函数
function err_c = calculate_err_c(a, b, t_c, c, f, t_th, F_W)

    % a为尺度参数，b为形状参数，tc为计算时间，c为计算任务量，f为CPU频率，t_th为阈值

    x = max((t_c - c/f - t_th),0);

    if b == 0
        G = exp(-x/a);
    else
        G = 1 - (1 + (b*x)/a).^(-1/b);
    end

    err_c = (1 - F_W)*(1 - G);

end
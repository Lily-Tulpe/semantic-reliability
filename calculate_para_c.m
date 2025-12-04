% 通信相关参数计算
function [C, V] = calculate_para_c(P_c, h_c, sigma2)

    % Pc为通信用户功率，hc为通信信道系数，sigma2为噪声功率

    gamma = (P_c*(h_c^2))/(sigma2); 
    V = 1-1/((1+gamma)^2); 
    C = log2(1+gamma); 

end
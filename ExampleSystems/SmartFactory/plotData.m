hold on;

plot(log2.I_w_HL_Weg(1:200));
plot(log2.I_w_HR_Weg(1:200));
plot(log2.I_w_BRU_Weg(1:200));
plot(log2.I_w_BHR_Weg(1:200));
plot(log2.I_w_BHL_Weg(1:200));
plot(log2.I_w_BLO_Weg(1:200));

% plot(log2.O_w_HAL_Ctrl(1:200));
% plot(log2.O_w_HAR_Ctrl(1:200));
% plot(100*log2.O_w_BRU_Axis_Ctrl(1:200));
% plot(100*log2.O_w_BHR_Axis_Ctrl(1:200));
% plot(100*log2.O_w_BHL_Axis_Ctrl(1:200));
% plot(100*log2.O_w_BLO_Axis_Ctrl(1:200));

legend('HL','HR','BRU','BHR','BHL','BLO');
% legend('Weg','Ctrl');

hold off;
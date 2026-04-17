%% 绘制 Ideal Power Grid 中所有 .mat 文件的导纳 Bode 图
% 达到科学研究（SCI论文）的绘图标准
clc;
clear;
close all;

%% 1. 文件夹与配置设置
Path_root_Results = 'Ideal Power Grid'; % 目标文件夹
if ~isfolder(Path_root_Results)
    error('未找到文件夹：%s', Path_root_Results);
end

matFiles = dir(fullfile(Path_root_Results, '*.mat'));
numFiles = numel(matFiles);
if numFiles == 0
    error('在 %s 中没有找到任何 .mat 文件！', Path_root_Results);
end
fprintf('共找到 %d 个 .mat 文件，开始绘制等效导纳 Y11 (从 v_alpha 到 i_galpha) 的Bode图...\n', numFiles);

%% 2. 频率范围设定
fmin = 1;       % 最小频率 Hz (根据并网变流器特性，可从 1Hz 开始)
fmax = 10000;    % 最大频率 Hz (如需观测更高频可改为 10000)
Npts = 1000;
f_vec = logspace(log10(fmin), log10(fmax), Npts);
w_vec = 2 * pi * f_vec;

%% 3. 创建高质量图窗
fig = figure('Name', 'Admittance Y11 Bode Plot', ...
             'Color', 'w', ...
             'Units', 'pixels', ...
             'Position', [100, 100, 800, 600]); 

% 使用好看的色图（如 parula 或 渐变色）避免线条混杂
colors = parula(max(numFiles, 2)); 

%% 4. 循环处理并绘制 Bode 图
ax_mag = subplot(2, 1, 1);
hold(ax_mag, 'on');
ax_phs = subplot(2, 1, 2);
hold(ax_phs, 'on');

poles_count = 0;
zeros_count = 0;

for k = 1:numFiles
    matPath = fullfile(matFiles(k).folder, matFiles(k).name);
    try
        data = load(matPath, 'sys');
        if ~isfield(data, 'sys')
            warning('文件 %s 中没有找到 sys 变量，已跳过。', matFiles(k).name);
            continue;
        end
        
        % 提取 Y11：输入1 (v_alpha) 到 输出1 (i_galpha)
        % 取决于 PMSGlinearization.m: u = [v_alpha, v_beta, omega_best, i_gqref]， y = [i_galpha, i_gbeta]
        sys_Y11 = data.sys(1, 1);
        
        % 只有第一个循环我们统计系统极点与零点个数，以便显示
        if poles_count == 0
            % 使用 minreal 去除不可观/不可控的状态，使得零极点分析更具物理意义
            sys_min = minreal(sys_Y11, 1e-4); 
            p_sys = pole(sys_min);
            z_sys = tzero(sys_min);
            poles_count = length(p_sys);
            zeros_count = length(z_sys);
        end
        
        % 计算频率响应
        [mag, phase, ~] = bode(sys_Y11, w_vec);
        mag = squeeze(mag);      % 转换为一维
        phase = squeeze(phase);  % 转换为一维
        
        % 将相位限制在 -180° 到 180° 之间
        phase = mod(phase + 180, 360) - 180;
        
        mag_dB = 20 * log10(mag);
        
        % 绘制幅频与相频
        plot(ax_mag, f_vec, mag_dB, 'Color', [colors(k, 1:3), 0.7], 'LineWidth', 1.5);
        plot(ax_phs, f_vec, phase, 'Color', [colors(k, 1:3), 0.7], 'LineWidth', 1.5);
        
    catch ME
        warning('读取或处理文件 %s 时出错：%s', matFiles(k).name, ME.message);
    end
end

%% 5. 设置坐标轴与美化 (科学文献标准)
% 幅频图设置
set(ax_mag, 'XScale', 'log', 'XMinorTick', 'on');
grid(ax_mag, 'on');
ylabel(ax_mag, 'Magnitude (dB)', 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
set(ax_mag, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.2);
ax_mag.XColor = 'k';
ax_mag.YColor = 'k';

% 相频图设置
set(ax_phs, 'XScale', 'log', 'XMinorTick', 'on');
grid(ax_phs, 'on');
ylim(ax_phs, [-180 180]);
yticks(ax_phs, -180:45:180);
xlabel(ax_phs, 'Frequency (Hz)', 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
ylabel(ax_phs, 'Phase (deg)', 'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
set(ax_phs, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.2);
ax_phs.XColor = 'k';
ax_phs.YColor = 'k';

%% 6. 显示零极点个数信息文本
% 将信息框放置在幅频图的左下角或右上角
info_str = sprintf('Number of Poles: %d\nNumber of Zeros: %d', poles_count, zeros_count);
annotation('textbox', [0.15 0.82 0.2 0.08], ...
    'String', info_str, ...
    'EdgeColor', 'k', ...
    'BackgroundColor', 'w', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'FitBoxToText', 'on');

%% 7. 调整布局和保存
linkaxes([ax_mag, ax_phs], 'x');
xlim(ax_mag, [fmin, fmax]);

% 避免刻度重叠
sgtitle('Admittance Y_{11} Bode Plot of Grid-Tied PMSG', 'FontName', 'Times New Roman', 'FontSize', 16, 'FontWeight', 'bold');

% 如果有必要，可解除下面两行的注释以自动保存高分辨率图片
% saveas(fig, fullfile(Path_root_Results, 'Admittance_Bode_Y11.png'));
% print(fig, fullfile(Path_root_Results, 'Admittance_Bode_Y11.eps'), '-dpng', '-r300');

disp('--> 绘图完成！');

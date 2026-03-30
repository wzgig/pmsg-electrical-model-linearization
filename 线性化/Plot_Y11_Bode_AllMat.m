function Plot_Y11_Bode_AllMat()
    clc; close all;

    % ====== 结果路径（与你主脚本一致）======
    Path_root_Results = "Single_variable_response_for_xi";

    if ~isfolder(Path_root_Results)
        error("Folder not found: %s", Path_root_Results);
    end

    matFiles = dir(fullfile(Path_root_Results, "*.mat"));
    if isempty(matFiles)
        error("No .mat files found in: %s", Path_root_Results);
    end

    % ====== 频率设置（你可以按需要修改）======
    % 建议用对数频率，单位 Hz
    fmin = 1e-1;     % Hz
    fmax = 1e5;     % Hz
    Npts = 2000;
    w = 2*pi*logspace(log10(fmin), log10(fmax), Npts);  % rad/s

    % bode 绘图选项
    opts = bodeoptions;
    opts.Grid = "on";
    opts.FreqUnits = "Hz";
    opts.PhaseWrapping = "on";

    % 输出子文件夹（可选）
    outDir = fullfile(Path_root_Results, "Y11_Bode");
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    fprintf("Found %d mat files in %s\n", numel(matFiles), Path_root_Results);
    fprintf("Saving figures to: %s\n\n", outDir);

    okCount = 0;
    failCount = 0;

    for k = 1:numel(matFiles)
        matPath = fullfile(matFiles(k).folder, matFiles(k).name);
        [~, baseName, ~] = fileparts(matFiles(k).name);

        try
            S = load(matPath);

            if ~isfield(S, "sys")
                error("Variable 'sys' not found in mat file.");
            end

            sys = S.sys;

            % ====== Y11 定义：默认取 sys(1,1) ======
            % 如果你的 Y11 通道不是(1,1)，改这里即可：
            Y11 = sys(1,1);

            % ====== 绘制并保存 ======
            h = figure('Visible','off','Color','w');
            bode(Y11, w, opts);
            titleStr = "Y_{11} Bode: " + baseName;
            title(titleStr, 'Interpreter','tex');

            % 保存格式（png + fig）
            pngPath = fullfile(outDir, baseName + "_Y11_Bode.png");
            figPath = fullfile(outDir, baseName + "_Y11_Bode.fig");

            exportgraphics(h, pngPath, "Resolution", 300);
            savefig(h, figPath);

            close(h);

            okCount = okCount + 1;
            fprintf("[OK] %s -> saved\n", baseName);

        catch ME
            failCount = failCount + 1;
            fprintf("[FAIL] %s : %s\n", baseName, ME.message);
        end
    end

    fprintf("\n====================================================\n");
    fprintf("Done. OK = %d, FAIL = %d\n", okCount, failCount);
    fprintf("Figures saved in: %s\n", outDir);
    fprintf("====================================================\n");
end

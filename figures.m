% Plot figures for Archibald et al. (2024), The American Naturalist
% 9.9.0.1538559 (R2020b) Update 3

%% Fig. 2

% Load simulation data
load fig2_data

% Load empirical data
exp_data.growth.high = readtable('Experimental Data\growth_highlight.csv', 'ReadRowNames', 1);
exp_data.photo.high = readtable('Experimental Data\photo_highlight.csv', 'ReadRowNames', 1);
exp_data.graze.high = readtable('Experimental Data\graze_highlight.csv', 'ReadRowNames', 1);
exp_data.resp.high = readtable('Experimental Data\resp_highlight.csv', 'ReadRowNames', 1);
exp_data.cue.high = readtable('Experimental Data\cue_highlight.csv', 'ReadRowNames', 1);

% Extract data indices
ki = repmat(Tevo,length(T),1) == repmat(T',1,length(Tevo));

%% Setup figure frame
figure('position', [1001 50 862 1306])
h = tiledlayout(5,2);
xlabel(h, 'Temperature (\circC)', 'fontname', 'arial', 'fontsize', 14)

% Plot 
nexttile(1)
hold on
errorbar(exp_data.resp.high.Acc(2:4), exp_data.resp.high.Resp(2:4), exp_data.resp.high.se(2:4), 'linestyle', 'none', 'color', [0.6 0.6 0.6], 'CapSize', 0)
errorbar(exp_data.resp.high.Acc(3), exp_data.resp.high.Resp(3), exp_data.resp.high.se(3), 'linestyle', 'none', 'color', 'k', 'CapSize', 0)
errorbar(exp_data.resp.high.Acc(1), exp_data.resp.high.Resp(1), exp_data.resp.high.se(1), 'linestyle', 'none', 'color', 'b', 'CapSize', 0)
errorbar(exp_data.resp.high.Acc(5), exp_data.resp.high.Resp(5), exp_data.resp.high.se(5), 'linestyle', 'none', 'color', 'r', 'CapSize', 0)
plot(exp_data.resp.high.Acc(2:4), exp_data.resp.high.Resp(2:4), '--o', 'color',  [0.6 0.6 0.6], 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.resp.high.Acc([1 3 5]), exp_data.resp.high.Resp([1 3 5]), '--ok', 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.resp.high.Acc(1), exp_data.resp.high.Resp(1), 'o', 'markerfacecolor', 'b', 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.resp.high.Acc(5), exp_data.resp.high.Resp(5), 'o', 'markerfacecolor', 'r', 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.resp.high.Acc(2), exp_data.resp.high.Resp(2), 'o', 'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.resp.high.Acc(4), exp_data.resp.high.Resp(4), 'o', 'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', 'none', 'linewidth', 1)
hold off
xlim([16 32])
ylim([0 3])
set(gca, 'xtick', [18 24 30])
set(gca, 'ytick', [0 1 2 3])
ylabel({'Respiration rate', '(gC gC^{-1} d^{-1})'})
title('Empirical Data', 'FontWeight', 'normal')
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(17, 3, 'A', 'fontname', 'arial', 'fontsize', 12, 'FontWeight', 'bold')

nexttile(2)
hold on
plot(T, cellfun(@(x) getfield(x, 'resp'), dn(:,Tevo==24)), '-', 'color', [0.6 0.6 0.6], 'linewidth', 1)
plot(Tevo,  cellfun(@(x) getfield(x, 'resp'), dn(ki)), '-k', 'linewidth', 1)
plot(18, dn{T==18,Tevo==18}.resp, 'o', 'markerfacecolor', 'b', 'markeredgecolor', 'none', 'linewidth', 1)
plot(24, dn{T==24,Tevo==24}.resp, 'o', 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'linewidth', 1)
plot(30, dn{T==30,Tevo==30}.resp, 'o', 'markerfacecolor', 'r', 'markeredgecolor', 'none', 'linewidth', 1)
plot(18, dn{T==18,Tevo==24}.resp, 'o', 'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', 'none', 'linewidth', 1)
plot(30, dn{T==30,Tevo==24}.resp, 'o', 'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', 'none', 'linewidth', 1)
hold off
xlim([16 32])
ylim([0 3])
set(gca, 'xtick', [18 24 30])
set(gca, 'ytick', [0 1 2 3])
title('Simulation Results', 'FontWeight', 'normal')
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(17, 3, 'B', 'fontname', 'arial', 'fontsize', 12, 'FontWeight', 'bold')

nexttile(3)
hold on
errorbar(exp_data.cue.high.Acc(4:6), exp_data.cue.high.CUE(4:6), exp_data.cue.high.se(4:6), 'linestyle', 'none', 'color', [0.6 0.6 0.6], 'CapSize', 0)
errorbar(exp_data.cue.high.Acc([1 5 9]), exp_data.cue.high.CUE([1 5 9]), exp_data.cue.high.se([1 5 9]), 'linestyle', 'none', 'color', 'k', 'CapSize', 0)
errorbar(exp_data.cue.high.Acc(1), exp_data.cue.high.CUE(1), exp_data.cue.high.se(1), 'linestyle', 'none', 'color', 'b', 'CapSize', 0)
errorbar(exp_data.cue.high.Acc(9), exp_data.cue.high.CUE(9), exp_data.cue.high.se(9), 'linestyle', 'none', 'color', 'r', 'CapSize', 0)
plot(exp_data.cue.high.Acc(4:6), exp_data.cue.high.CUE(4:6), '--o', 'color',  [0.6 0.6 0.6], 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.cue.high.Acc([1 5 9]), exp_data.cue.high.CUE([1 5 9]), '--ok', 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.cue.high.Acc(1), exp_data.cue.high.CUE(1), 'o', 'markerfacecolor', 'b', 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.cue.high.Acc(9), exp_data.cue.high.CUE(9), 'o', 'markerfacecolor', 'r', 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.cue.high.Acc(4), exp_data.cue.high.CUE(4), 'o', 'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.cue.high.Acc(6), exp_data.cue.high.CUE(6), 'o', 'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', 'none', 'linewidth', 1)
hold off
xlim([16 32])
ylim([0 0.6])
set(gca, 'xtick', [18 24 30])
set(gca, 'ytick', [0 0.2 0.4 0.6])
ylabel('Carbon use efficiency')
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(17, 0.6, 'C', 'fontname', 'arial', 'fontsize', 12, 'FontWeight', 'bold')

nexttile(4)
hold on
plot(T, cellfun(@(x) getfield(x, 'cue'), dn(:,Tevo==24)), '-', 'color', [0.6 0.6 0.6], 'linewidth', 1)
plot(Tevo,  cellfun(@(x) getfield(x, 'cue'), dn(ki)), '-k', 'linewidth', 1)
plot(18, dn{T==18,Tevo==18}.cue, 'o', 'markerfacecolor', 'b', 'markeredgecolor', 'none', 'linewidth', 1)
plot(24, dn{T==24,Tevo==24}.cue, 'o', 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'linewidth', 1)
plot(30, dn{T==30,Tevo==30}.cue, 'o', 'markerfacecolor', 'r', 'markeredgecolor', 'none', 'linewidth', 1)
plot(18, dn{T==18,Tevo==24}.cue, 'o', 'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', 'none', 'linewidth', 1)
plot(30, dn{T==30,Tevo==24}.cue, 'o', 'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', 'none', 'linewidth', 1)
hold off
xlim([16 32])
ylim([0 0.6])
set(gca, 'xtick', [18 24 30])
set(gca, 'ytick', [0 0.2 0.4 0.6])
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(17, 0.6, 'D', 'fontname', 'arial', 'fontsize', 12, 'FontWeight', 'bold')

nexttile(5)
hold on
errorbar(exp_data.photo.high.Acc(2:4), exp_data.photo.high.Photo(2:4), exp_data.photo.high.se(2:4), 'linestyle', 'none', 'color', [0.6 0.6 0.6], 'CapSize', 0)
errorbar(exp_data.photo.high.Acc([1 3 5]), exp_data.photo.high.Photo([1 3 5]), exp_data.photo.high.se([1 3 5]), 'linestyle', 'none', 'color', 'k', 'CapSize', 0)
errorbar(exp_data.photo.high.Acc(1), exp_data.photo.high.Photo(1), exp_data.photo.high.se(1), 'linestyle', 'none', 'color', 'b', 'CapSize', 0)
errorbar(exp_data.photo.high.Acc(5), exp_data.photo.high.Photo(5), exp_data.photo.high.se(5), 'linestyle', 'none', 'color', 'r', 'CapSize', 0)
plot(exp_data.photo.high.Acc(2:4), exp_data.photo.high.Photo(2:4), '--o', 'color',  [0.6 0.6 0.6], 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.photo.high.Acc([1 3 5]), exp_data.photo.high.Photo([1 3 5]), '--ok', 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.photo.high.Acc(1), exp_data.photo.high.Photo(1), 'o', 'markerfacecolor', 'b', 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.photo.high.Acc(5), exp_data.photo.high.Photo(5), 'o', 'markerfacecolor', 'r', 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.photo.high.Acc(2), exp_data.photo.high.Photo(2), 'o', 'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.photo.high.Acc(4), exp_data.photo.high.Photo(4), 'o', 'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', 'none', 'linewidth', 1)
hold off
xlim([16 32])
ylim([0 1.5])
set(gca, 'xtick', [18 24 30])
ylabel({'Photosynthesis', '(gC gC^{-1} d^{-1})'})
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(17, 1.5, 'E', 'fontname', 'arial', 'fontsize', 12, 'FontWeight', 'bold')

nexttile(6)
hold on
plot(T, cellfun(@(x) getfield(x, 'photo'), dn(:,Tevo==24)), '-', 'color', [0.6 0.6 0.6], 'linewidth', 1)
plot(Tevo,  cellfun(@(x) getfield(x, 'photo'), dn(ki)), '-k', 'linewidth', 1)
plot(18, dn{T==18,Tevo==18}.photo, 'o', 'markerfacecolor', 'b', 'markeredgecolor', 'none', 'linewidth', 1)
plot(24, dn{T==24,Tevo==24}.photo, 'o', 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'linewidth', 1)
plot(30, dn{T==30,Tevo==30}.photo, 'o', 'markerfacecolor', 'r', 'markeredgecolor', 'none', 'linewidth', 1)
plot(18, dn{T==18,Tevo==24}.photo, 'o', 'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', 'none', 'linewidth', 1)
plot(30, dn{T==30,Tevo==24}.photo, 'o', 'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', 'none', 'linewidth', 1)
hold off
xlim([16 32])
ylim([0 1.5])
set(gca, 'xtick', [18 24 30])
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(17, 1.5, 'F', 'fontname', 'arial', 'fontsize', 12, 'FontWeight', 'bold')

nexttile(7)
hold on
errorbar(exp_data.graze.high.Acc(2:4), exp_data.graze.high.graz(2:4), exp_data.graze.high.se(2:4), 'linestyle', 'none', 'color', [0.6 0.6 0.6], 'CapSize', 0)
errorbar(exp_data.graze.high.Acc([1 3 5]), exp_data.graze.high.graz([1 3 5]), exp_data.graze.high.se([1 3 5]), 'linestyle', 'none', 'color', 'k', 'CapSize', 0)
errorbar(exp_data.graze.high.Acc(1), exp_data.graze.high.graz(1), exp_data.graze.high.se(1), 'linestyle', 'none', 'color', 'b', 'CapSize', 0)
errorbar(exp_data.graze.high.Acc(5), exp_data.graze.high.graz(5), exp_data.graze.high.se(5), 'linestyle', 'none', 'color', 'r', 'CapSize', 0)
plot(exp_data.graze.high.Acc(2:4), exp_data.graze.high.graz(2:4), '--o', 'color',  [0.6 0.6 0.6], 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.graze.high.Acc([1 3 5]), exp_data.graze.high.graz([1 3 5]), '--ok', 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.graze.high.Acc(1), exp_data.graze.high.graz(1), 'o', 'markerfacecolor', 'b', 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.graze.high.Acc(5), exp_data.graze.high.graz(5), 'o', 'markerfacecolor', 'r', 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.graze.high.Acc(2), exp_data.graze.high.graz(2), 'o', 'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.graze.high.Acc(4), exp_data.graze.high.graz(4), 'o', 'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', 'none', 'linewidth', 1)
hold off
xlim([16 32])
ylim([0 0.4])
set(gca, 'xtick', [18 24 30])
ylabel('Grazing rate (d^{-1})')
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(17, 0.4, 'G', 'fontname', 'arial', 'fontsize', 12, 'FontWeight', 'bold')

nexttile(8)
hold on
plot(T, cellfun(@(x) getfield(x, 'graz'), dn(:,Tevo==24)), '-', 'color', [0.6 0.6 0.6], 'linewidth', 1)
plot(Tevo,  cellfun(@(x) getfield(x, 'graz'), dn(ki)), '-k', 'linewidth', 1)
plot(18, dn{T==18,Tevo==18}.graz, 'o', 'markerfacecolor', 'b', 'markeredgecolor', 'none', 'linewidth', 1)
plot(24, dn{T==24,Tevo==24}.graz, 'o', 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'linewidth', 1)
plot(30, dn{T==30,Tevo==30}.graz, 'o', 'markerfacecolor', 'r', 'markeredgecolor', 'none', 'linewidth', 1)
plot(18, dn{T==18,Tevo==24}.graz, 'o', 'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', 'none', 'linewidth', 1)
plot(30, dn{T==30,Tevo==24}.graz, 'o', 'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', 'none', 'linewidth', 1)
hold off
xlim([16 32])
ylim([0 0.4])
set(gca, 'xtick', [18 24 30])
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(17, 0.4, 'H', 'fontname', 'arial', 'fontsize', 12, 'FontWeight', 'bold')

nexttile(9)
hold on
errorbar(exp_data.growth.high.Acc(4:6), exp_data.growth.high.Growth_Rate(4:6), exp_data.growth.high.se(4:6), 'linestyle', 'none', 'color', [0.6 0.6 0.6], 'CapSize', 0)
errorbar(exp_data.growth.high.Acc(1), exp_data.growth.high.Growth_Rate(1), exp_data.growth.high.se(1), 'linestyle', 'none', 'color', 'k', 'CapSize', 0)
errorbar(exp_data.growth.high.Acc(1), exp_data.growth.high.Growth_Rate(1), exp_data.growth.high.se(1), 'linestyle', 'none', 'color', 'b', 'CapSize', 0)
errorbar(exp_data.growth.high.Acc(9), exp_data.growth.high.Growth_Rate(9), exp_data.growth.high.se(9), 'linestyle', 'none', 'color', 'r', 'CapSize', 0)
plot(exp_data.growth.high.Acc(4:6), exp_data.growth.high.Growth_Rate(4:6), '--o', 'color',  [0.6 0.6 0.6], 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.growth.high.Acc([1 5 9]), exp_data.growth.high.Growth_Rate([1 5 9]), '--ok', 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.growth.high.Acc(1), exp_data.growth.high.Growth_Rate(1), 'o', 'markerfacecolor', 'b', 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.growth.high.Acc(9), exp_data.growth.high.Growth_Rate(9), 'o', 'markerfacecolor', 'r', 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.growth.high.Acc(4), exp_data.growth.high.Growth_Rate(4), 'o', 'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', 'none', 'linewidth', 1)
plot(exp_data.growth.high.Acc(6), exp_data.growth.high.Growth_Rate(6), 'o', 'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', 'none', 'linewidth', 1)
hold off
xlim([16 32])
ylim([0 0.6])
set(gca, 'xtick', [18 24 30])
set(gca, 'ytick', [0 0.2 0.4 0.6])
ylabel('Net growth rate (d^{-1})')
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(17, 0.6, 'I', 'fontname', 'arial', 'fontsize', 12, 'FontWeight', 'bold')

nexttile(10)
hold on
plot(T, cellfun(@(x) getfield(x, 'ngr'), dn(:,Tevo==24)), '-', 'color', [0.6 0.6 0.6], 'linewidth', 1)
plot(Tevo,  cellfun(@(x) getfield(x, 'ngr'), dn(ki)), '-k', 'linewidth', 1)
plot(18, dn{T==18,Tevo==18}.ngr, 'o', 'markerfacecolor', 'b', 'markeredgecolor', 'none', 'linewidth', 1)
plot(24, dn{T==24,Tevo==24}.ngr, 'o', 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'linewidth', 1)
plot(30, dn{T==30,Tevo==30}.ngr, 'o', 'markerfacecolor', 'r', 'markeredgecolor', 'none', 'linewidth', 1)
plot(18, dn{T==18,Tevo==24}.ngr, 'o', 'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', 'none', 'linewidth', 1)
plot(30, dn{T==30,Tevo==24}.ngr, 'o', 'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', 'none', 'linewidth', 1)
hold off
xlim([16 32])
ylim([0 0.6])
set(gca, 'ytick', [0 0.2 0.4 0.6])
set(gca, 'xtick', [18 24 30])
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(17, 0.6, 'J', 'fontname', 'arial', 'fontsize', 12, 'FontWeight', 'bold')

%% Fig. 3

load fig3_data_tinc_con
load fig3_data_tinc_evo
load fig3_data_tdec_con
load fig3_data_tdec_evo

figure('position', [1001 159 1217 1180])
ht = tiledlayout(6,2);
xlabel(ht, 'Days')

nexttile(1)
hold on
plot(t, Tinc, 'color', 'red')
plot(dn_tinc_con.t, dn_tinc_con.that, 'color', [0.6 0.6 0.6])
plot(dn_tinc_evo.t, dn_tinc_evo.that, 'color', 'black')
hold off
ylim([17 31])
set(gca, 'ytick', [18 24 30])
ylabel('^\circC')
legend('Temperature', 'Control T_{evo}', 'Evolving T_{evo}', 'location', 'southeast', 'fontname', 'arial', 'fontsize', 8)
legend boxoff
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(5, max(ylim), 'A', 'fontweight', 'bold', 'fontsize', 12)

nexttile(3)
hold on
semilogy(dn_tinc_con.t, dn_tinc_con.x, 'color', [0.6 0.6 0.6])
semilogy(dn_tinc_evo.t, dn_tinc_evo.x, 'color', 'black')
hold off
set(gca, 'yticklabels', {'0', '5', '10', '15'});
ylabel({'Abundance', '(10^5 cells mL^{-1})'})
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(5, max(ylim), 'C', 'fontweight', 'bold', 'fontsize', 12)

nexttile(5)
hold on
plot(dn_tinc_con.t, dn_tinc_con.ngr, 'color', [0.6 0.6 0.6])
plot(dn_tinc_evo.t, dn_tinc_evo.ngr, 'color', 'black')
scatter(198, dn_tinc_evo.ngr(end), [], 'r', 'filled')
scatter(55, dn_tinc_con.ngr(end), [], [0.6 0.6 0.6], 'filled')
hold off
ylabel('Growth rate (d^{-1})')
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(5, max(ylim), 'E', 'fontweight', 'bold', 'fontsize', 12)

nexttile(7)
hold on
plot(dn_tinc_con.t, dn_tinc_con.resp, 'color', [0.6 0.6 0.6])
plot(dn_tinc_evo.t, dn_tinc_evo.resp, 'color', 'black')
scatter(198, dn_tinc_evo.resp(end), [], 'r', 'filled')
scatter(55, dn_tinc_con.resp(end), [], [0.6 0.6 0.6], 'filled')
hold off
ylim([0 1.5])
ylabel({'Respiration rate', '(gC gC^{-1} d^{-1})'})
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(5, max(ylim), 'G', 'fontweight', 'bold', 'fontsize', 12)

nexttile(9)
hold on
plot(dn_tinc_con.t, dn_tinc_con.strat(:,2), 'color', [0.6 0.6 0.6])
plot(dn_tinc_evo.t, dn_tinc_evo.strat(:,2), 'color', 'black')
hold off
ylim([0 0.5])
set(gca, 'ytick', [0 0.5])
ylabel({'Photosynthesis', 'investment (\rho_p)'})
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(5, max(ylim), 'I', 'fontweight', 'bold', 'fontsize', 12)

nexttile(11)
hold on
plot(dn_tinc_con.t, dn_tinc_con.strat(:,3), 'color', [0.6 0.6 0.6])
plot(dn_tinc_evo.t, dn_tinc_evo.strat(:,3), 'color', 'black')
hold off
ylim([0.5 1])
set(gca, 'ytick', [0.5 1])
ylabel({'Grazing', 'investment (\rho_g)'})
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(5, max(ylim), 'K', 'fontweight', 'bold', 'fontsize', 12)

nexttile(2)
hold on
plot(t, Tdec, 'color', 'blue')
plot(dn_tdec_con.t, dn_tdec_con.that, 'color', [0.6 0.6 0.6])
plot(dn_tdec_evo.t, dn_tdec_evo.that, 'color', 'black')
hold off
ylim([17 31])
set(gca, 'ytick', [18 24 30])
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(5, max(ylim), 'B', 'fontweight', 'bold', 'fontsize', 12)

nexttile(4)
hold on
semilogy(dn_tdec_con.t, dn_tdec_con.x, 'color', [0.6 0.6 0.6])
semilogy(dn_tdec_evo.t, dn_tdec_evo.x, 'color', 'black')
hold off
ylim([0 15e5])
set(gca, 'yticklabels', {'0', '5', '10', '15'})
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(5, max(ylim), 'D', 'fontweight', 'bold', 'fontsize', 12)

nexttile(6)
hold on
plot(dn_tdec_con.t, dn_tdec_con.ngr, 'color', [0.6 0.6 0.6])
plot(dn_tdec_evo.t, dn_tdec_evo.ngr, 'color', 'black')
scatter(198, dn_tdec_evo.ngr(end), [], 'b', 'filled')
scatter(55, dn_tdec_con.ngr(end), [], [0.6 0.6 0.6], 'filled')
hold off
ylim([-0.5 0.5])
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(5, max(ylim), 'F', 'fontweight', 'bold', 'fontsize', 12)

nexttile(8)
hold on
plot(dn_tdec_con.t, dn_tdec_con.resp, 'color', [0.6 0.6 0.6])
plot(dn_tdec_evo.t, dn_tdec_evo.resp, 'color', 'black')
scatter(198, dn_tdec_evo.resp(end), [], 'b', 'filled')
scatter(55, dn_tdec_con.resp(end), [], [0.6 0.6 0.6], 'filled')
hold off
ylim([0 1.5])
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(5, max(ylim), 'H', 'fontweight', 'bold', 'fontsize', 12)

nexttile(10)
hold on
plot(dn_tdec_con.t, dn_tdec_con.strat(:,2), 'color', [0.6 0.6 0.6])
plot(dn_tdec_evo.t, dn_tdec_evo.strat(:,2), 'color', 'black')
hold off
ylim([0 0.5])
set(gca, 'ytick', [0 0.5])
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(5, max(ylim), 'J', 'fontweight', 'bold', 'fontsize', 12)

nexttile(12)
hold on
plot(dn_tdec_con.t, dn_tdec_con.strat(:,3), 'color', [0.6 0.6 0.6])
plot(dn_tdec_evo.t, dn_tdec_evo.strat(:,3), 'color', 'black')
hold off
ylim([0.5 1])
set(gca, 'ytick', [0.5 1])
set(gca, 'fontname', 'arial', 'fontsize', 12)
text(5, max(ylim), 'L', 'fontweight', 'bold', 'fontsize', 12)

%% Fig. S1

load figs1

plabel = {'(a) u_B', '(b) Q_{10,B}', '(c) u_{max}', '(d) Q_{10,U}', '(e) v_P', '(f) Q_{10,P}', '(g) v_G', '(h) Q_{10,G}', '(i) r_a', '(j) r_0', '(k) Q_{10,R}', ...
    '(l) \lambda', '(m) \gamma', '(n) Q_{min,C}', '(o) Q_{min,N}', '(p) b_C', '(q) b_N', '(r) k_I', '(s) k_B', '(t) c'};

figure('position', [1000 310 631 1037])
ht = tiledlayout(7, 3, 'tilespacing', 'compact');
xlabel(ht, 'Temperature (\circC)', 'fontname', 'arial', 'fontsize', 12)
ylabel(ht, 'Net growth rate (d^{-1})', 'fontname', 'arial', 'fontsize', 12)
for i = 1:length(plabel)
    nexttile
    plot(T, diag(uplus(:,:,i+1)), '-g', T, diag(uminus(:,:,i+1)), '-m', T,  diag(u), '-k', 'linewidth', 1.5)
    text(18.5, 0.7, plabel{i}, 'fontname', 'arial', 'fontsize', 12)
    xlim([18 30])
    ylim([0 0.8])
    h = gca;
    h.YAxis.TickValues = [0 0.4 0.8];
    if ~ismember(i,1:3:19); h.YAxis.TickLabels = {}; end
    h.XAxis.TickValues = [ 18 24 30];
    if ~ismember(i, 19:21); h.XAxis.TickLabels = {}; end
    set(gca, 'fontname', 'arial', 'fontsize', 12)
end


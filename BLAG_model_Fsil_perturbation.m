%% =====================================================================
%  BLAG_model_Fsil_perturbation.m
%  =====================================================================
%  Two-phase BLAG model:
%    Phase 1 (570–424 Ma): Full BLAG feedback spinup — the system finds
%       its own self-consistent equilibrium with CO2-weathering feedback.
%    Phase 2 (424–300 Ma): Prescribed silicate weathering from Li isotope
%       data. The WI multiplier perturbs F_sil relative to the spinup
%       end-state. CO2 and other variables respond via the full BLAG
%       equations (Eqs. 59A-59H), except silicate weathering.
%
%  All equations from Berner, Lasaga & Garrels (1983).
%  Reference: Am. J. Sci. 283, 641-683.
%  Li forcing: Ghosh (2026) Li box model output.
% =====================================================================

clear; clc; close all;
fprintf('=== BLAG Model: Two-Phase Prescribed F_sil ===\n');
fprintf('    Phase 1: Spinup 570–424 Ma (BLAG feedback)\n');
fprintf('    Phase 2: Forcing 424–300 Ma (Li-prescribed F_sil)\n\n');

%% =====================================================================
%  SECTION 1: LOAD GEOCARB FORCING (fA, fSR)
%  =====================================================================
geocarb_file = fullfile('Data', 'GEOCARB_input_arrays_tMod.csv');
T_geocarb    = readtable(geocarb_file);

ages_raw = T_geocarb.age;
age_max  = max(ages_raw);  % 570 Ma

% Model time: t=0 at 570 Ma, t increases forward
mt = age_max - ages_raw;
[mt, si] = sort(mt);
flds = {'fA','fSR'};
for k = 1:numel(flds)
    raw.(flds{k}) = T_geocarb.(flds{k})(si);
end
[mt, ui] = unique(mt, 'stable');
for k = 1:numel(flds)
    raw.(flds{k}) = raw.(flds{k})(ui);
end

fprintf('Loaded %d GEOCARB forcing time points.\n', numel(mt));

%% =====================================================================
%  SECTION 2: LOAD Li BOX MODEL OUTPUT (WI DATA)
%  =====================================================================
Li_file = fullfile('Data', 'Li_model_output_Ghosh2026.xlsx');
T_Li    = readtable(Li_file);
fprintf('Loaded Li model output: %d samples.\n', height(T_Li));

Li_ages    = T_Li.Age_Ma;
Li_WI_high = T_Li.WI_high;
valid      = ~isnan(Li_WI_high);
Li_ages    = Li_ages(valid);
Li_WI_high = Li_WI_high(valid);

[Li_ages, sort_idx] = sort(Li_ages, 'descend');
Li_WI_high = Li_WI_high(sort_idx);

[Li_ages_unique, ~, ic] = unique(Li_ages, 'stable');
Li_WI_unique = accumarray(ic, Li_WI_high, [], @mean);
Li_ages    = Li_ages_unique;
Li_WI_high = Li_WI_unique;

WI_start = Li_WI_high(1);
WI_cumul = Li_WI_high / WI_start;

fprintf('Valid WI_high samples: %d\n', length(Li_ages));
fprintf('WI_high at %.0f Ma = %.4f (baseline)\n', Li_ages(1), WI_start);

%% =====================================================================
%  SECTION 3: BLAG PARAMETERS (Table 2, Eqs. 36-53B)
%  =====================================================================
% Present-day reservoir sizes
D0       = 1000;
C0       = 3000;
S_CaSi0  = 1e5;
S_MgSi0  = 1e5;
M_Mg0    = 75.0;
M_Ca0    = 14.0;
M_HCO3_0 = 2.8;
A_CO2_0  = 0.055;

% Fluxes (Table 2)
Fwsil_Ca0 = 2.72;
Fwsil_Mg0 = 3.02;
F_sil_total_0 = Fwsil_Ca0 + Fwsil_Mg0;

p.A_CO2_0      = A_CO2_0;
p.frac_CaSi    = Fwsil_Ca0 / F_sil_total_0;
p.frac_MgSi    = Fwsil_Mg0 / F_sil_total_0;
p.kw_D0        = 0.00198;
p.kw_C0        = 0.00261;
p.kw_CaSi0     = Fwsil_Ca0 / S_CaSi0;
p.kw_MgSi0     = Fwsil_Mg0 / S_MgSi0;
p.k_vsw0       = 5.00 / M_Mg0;    % 0.0667 my^-1
p.k_MD0        = 0.00144;
p.k_MC0        = 0.000953;
p.K_eq         = 1721;
p.k_prec0      = 1.1620;
p.runoff_coeff     = 0.038;
p.bicarb_coeff     = 0.049;
p.greenhouse_coeff = 0.347;
p.fB_coeff         = 0.25;
p.T0               = 288;

fprintf('\n=== BLAG Parameters ===\n');
fprintf('  K_eq=%.0f, k_prec=%.4f, k_prec*K_eq=%.0f (should be 2000)\n', ...
    p.K_eq, p.k_prec0, p.k_prec0*p.K_eq);
fprintf('  F_sil_total(0) = %.2f (Ca: %.2f + Mg: %.2f)\n', ...
    F_sil_total_0, Fwsil_Ca0, Fwsil_Mg0);

%% =====================================================================
%  SECTION 4: PHASE 1 — SPINUP (570–424 Ma, BLAG feedback)
%  =====================================================================
% Use full BLAG feedback so the system finds its own equilibrium.
% Model time: t=0 at 570 Ma, t=146 at 424 Ma.
%
% fA and fSR interpolants use the GEOCARB model time directly.

age_start_spinup = 570;
age_switch       = 424;
age_end          = 300;

t_switch = age_start_spinup - age_switch;  % 146 My
t_end    = age_start_spinup - age_end;     % 270 My

% Phase 1 parameters: use GEOCARB time directly (t=0 at 570 Ma)
p1 = p;
p1.interp_fA  = griddedInterpolant(mt, raw.fA,  'linear', 'nearest');
p1.interp_fSR = griddedInterpolant(mt, raw.fSR, 'linear', 'nearest');
p1.mode = 'feedback';  % full BLAG feedback

% Initial conditions at 570 Ma: present-day values (standard BLAG)
y0_spinup = [D0; C0; S_CaSi0; S_MgSi0; M_Mg0; M_Ca0; M_HCO3_0; A_CO2_0];

% Verify steady state at t=0
[dydt_check, ~] = BLAG_odes_Fsil(0, y0_spinup, p1);
fprintf('\n=== Steady-State Check at 570 Ma ===\n');
labels = {'dD/dt','dC/dt','dS_CaSi/dt','dS_MgSi/dt','dM_Mg/dt','dM_Ca/dt','dM_HCO3/dt','dA_CO2/dt'};
for i = 1:8
    fprintf('  %14s = %+.6e\n', labels{i}, dydt_check(i));
end

% Integrate Phase 1
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', 1.0, 'NonNegative', 1:8);
fprintf('\nPhase 1: Integrating 570–424 Ma (%d Myr, BLAG feedback)...\n', t_switch);
tic;
[t1, y1] = ode15s(@(t,y) BLAG_odes_Fsil(t, y, p1), [0, t_switch], y0_spinup, opts);
elapsed1 = toc;
fprintf('Done. %d steps in %.2f s.\n', length(t1), elapsed1);

% Report state at 424 Ma (end of spinup)
y_424 = y1(end,:)';
age1  = age_start_spinup - t1;
CO2_424_xPAL = y_424(8) / A_CO2_0;
CO2_424_ppm  = CO2_424_xPAL * 280;

fprintf('\n=== State at %.0f Ma (end of spinup) ===\n', age_switch);
fprintf('  D      = %.1f\n', y_424(1));
fprintf('  C      = %.1f\n', y_424(2));
fprintf('  S_CaSi = %.1e\n', y_424(3));
fprintf('  S_MgSi = %.1e\n', y_424(4));
fprintf('  M_Mg   = %.2f\n', y_424(5));
fprintf('  M_Ca   = %.4f\n', y_424(6));
fprintf('  M_HCO3 = %.4f\n', y_424(7));
fprintf('  A_CO2  = %.4f (%.0f ppm, %.1f xPAL)\n', y_424(8), CO2_424_ppm, CO2_424_xPAL);

% Post-process Phase 1
npts1   = length(t1);
n_flux  = 22;
fluxes1 = zeros(npts1, n_flux);
for i = 1:npts1
    [~, fi] = BLAG_odes_Fsil(t1(i), y1(i,:)', p1);
    fluxes1(i,:) = fi;
end

%% =====================================================================
%  SECTION 5: COMPUTE F_sil BASELINE AT 424 Ma
%  =====================================================================
% F_sil_baseline is the silicate weathering rate that the BLAG model
% naturally produced at 424 Ma during the spinup. We read it directly
% from the spinup end-state.

F_sil_baseline = fluxes1(end, 5);  % total silicate weathering at 424 Ma

fprintf('\n=== F_sil Baseline ===\n');
fprintf('  F_sil at 424 Ma (from spinup) = %.4f (10^18 mol/my)\n', F_sil_baseline);
fprintf('  Fw_CaSi = %.4f, Fw_MgSi = %.4f\n', fluxes1(end,3), fluxes1(end,4));
fprintf('  F_degas (met CO2) = %.4f\n', fluxes1(end,10));
fprintf('  F_prec  = %.4f\n', fluxes1(end,7));

%% =====================================================================
%  SECTION 6: BUILD WI INTERPOLANT FOR PHASE 2
%  =====================================================================
% Phase 2 model time: t2=0 at 424 Ma, t2 increases forward
Li_mt_phase2 = age_switch - Li_ages;

sim_mask     = (Li_mt_phase2 >= 0) & (Li_mt_phase2 <= (age_switch - age_end));
Li_mt_sim    = Li_mt_phase2(sim_mask);
WI_cumul_sim = WI_cumul(sim_mask);

t_end_phase2 = age_switch - age_end;  % 124 My

if Li_mt_sim(end) < t_end_phase2
    Li_mt_sim    = [Li_mt_sim;    t_end_phase2];
    WI_cumul_sim = [WI_cumul_sim; WI_cumul_sim(end)];
end
if Li_mt_sim(1) > 0
    Li_mt_sim    = [0;   Li_mt_sim];
    WI_cumul_sim = [1.0; WI_cumul_sim];
end
[Li_mt_sim, ui_mt] = unique(Li_mt_sim, 'stable');
WI_cumul_sim = WI_cumul_sim(ui_mt);

%% =====================================================================
%  SECTION 7: PHASE 2 — PRESCRIBED F_sil (424–300 Ma)
%  =====================================================================
% Use the spinup end-state as initial conditions.
% Switch silicate weathering to prescribed mode.

p2 = p;
p2.mode = 'prescribed';
p2.F_sil_baseline = F_sil_baseline;
p2.interp_WI_cumul = griddedInterpolant(Li_mt_sim, WI_cumul_sim, 'linear', 'nearest');

% Shift GEOCARB interpolants to Phase 2 local time (t2=0 at 424 Ma)
offset_phase2 = age_start_spinup - age_switch;  % 146
p2.interp_fA  = griddedInterpolant(mt - offset_phase2, raw.fA,  'linear', 'nearest');
p2.interp_fSR = griddedInterpolant(mt - offset_phase2, raw.fSR, 'linear', 'nearest');

% Verify time coordinates
fprintf('\n=== Phase 2 Time Check ===\n');
fprintf('  At t2=0 (%.0f Ma): fA=%.3f, fSR=%.3f\n', age_switch, p2.interp_fA(0), p2.interp_fSR(0));
fprintf('  At t2=%.0f (%.0f Ma): fA=%.3f, fSR=%.3f\n', ...
    t_end_phase2, age_end, p2.interp_fA(t_end_phase2), p2.interp_fSR(t_end_phase2));

% Initial conditions = spinup end-state
y0_phase2 = y_424;

% Verify near-steady-state at start of Phase 2
[dydt_check2, ~] = BLAG_odes_Fsil(0, y0_phase2, p2);
fprintf('\n=== Steady-State Check at 424 Ma (Phase 2 start) ===\n');
for i = 1:8
    fprintf('  %14s = %+.6e\n', labels{i}, dydt_check2(i));
end

% Integrate Phase 2
opts2 = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', 0.5, 'NonNegative', 1:8);
fprintf('\nPhase 2: Integrating 424–300 Ma (%d Myr, prescribed F_sil)...\n', t_end_phase2);
tic;
[t2, y2] = ode15s(@(t,y) BLAG_odes_Fsil(t, y, p2), [0, t_end_phase2], y0_phase2, opts2);
elapsed2 = toc;
fprintf('Done. %d steps in %.2f s.\n', length(t2), elapsed2);

% Post-process Phase 2
age2    = age_switch - t2;
npts2   = length(t2);
fluxes2 = zeros(npts2, n_flux);
for i = 1:npts2
    [~, fi] = BLAG_odes_Fsil(t2(i), y2(i,:)', p2);
    fluxes2(i,:) = fi;
end

CO2_abs2  = y2(:,8);
CO2_xPAL2 = CO2_abs2 / A_CO2_0;
CO2_ppm2  = CO2_xPAL2 * 280;
DeltaT2   = fluxes2(:,11);

%% =====================================================================
%  SECTION 8: COMBINE RESULTS
%  =====================================================================
% Merge Phase 1 and Phase 2 for full 570–300 Ma plots
age_all    = [age1; age2(2:end)];
y_all      = [y1; y2(2:end,:)];
fluxes_all = [fluxes1; fluxes2(2:end,:)];
CO2_abs_all  = y_all(:,8);
CO2_xPAL_all = CO2_abs_all / A_CO2_0;
CO2_ppm_all  = CO2_xPAL_all * 280;
DeltaT_all   = fluxes_all(:,11);

%% =====================================================================
%  SECTION 9: LOAD Li pCO2 FOR COMPARISON
%  =====================================================================
Li_pCO2     = T_Li.pCO2_median_ppm(valid);
Li_pCO2_p16 = T_Li.pCO2_p16_ppm(valid);
Li_pCO2_p84 = T_Li.pCO2_p84_ppm(valid);
Li_ages_all_v = T_Li.Age_Ma(valid);
Li_pCO2     = Li_pCO2(sort_idx);
Li_pCO2_p16 = Li_pCO2_p16(sort_idx);
Li_pCO2_p84 = Li_pCO2_p84(sort_idx);
Li_ages_all_v = Li_ages_all_v(sort_idx);

pCO2_cap = 20000;
Li_pCO2_p16 = min(Li_pCO2_p16, pCO2_cap);
Li_pCO2_p84 = min(Li_pCO2_p84, pCO2_cap);

Li_sim_mask     = (Li_ages_all_v >= age_end) & (Li_ages_all_v <= age_switch);
Li_ages_plot    = Li_ages_all_v(Li_sim_mask);
Li_pCO2_plot    = Li_pCO2(Li_sim_mask);
Li_pCO2_p16_plt = Li_pCO2_p16(Li_sim_mask);
Li_pCO2_p84_plt = Li_pCO2_p84(Li_sim_mask);

%% =====================================================================
%  SECTION 10: PLOTTING
%  =====================================================================

% --- Figure 1: Full 570–300 Ma Overview ---
figure('Name', 'BLAG Full Run (570–300 Ma)', 'Position', [40 40 1400 900]);

subplot(2,2,1);
plot(age_all, CO2_ppm_all, 'b-', 'LineWidth', 2); hold on;
xline(age_switch, 'r--', 'Li forcing on', 'LineWidth', 1.5, 'LabelOrientation', 'horizontal');
yline(280, 'k--', 'Present (280 ppm)', 'HandleVisibility', 'off');
scatter(Li_ages_plot, Li_pCO2_plot, 60, 'ro', 'filled', 'DisplayName', 'Li model');
for j = 1:length(Li_ages_plot)
    if ~isnan(Li_pCO2_p16_plt(j)) && ~isnan(Li_pCO2_p84_plt(j))
        plot([Li_ages_plot(j) Li_ages_plot(j)], ...
             [Li_pCO2_p16_plt(j) Li_pCO2_p84_plt(j)], 'r-', 'LineWidth', 0.8, ...
             'HandleVisibility', 'off');
    end
end
set(gca, 'XDir', 'reverse');
xlabel('Age (Ma)'); ylabel('pCO_2 (ppm)');
title('Atmospheric CO_2'); legend('Location', 'best'); grid on;

subplot(2,2,2);
plot(age_all, fluxes_all(:,5), 'b-', 'LineWidth', 2, 'DisplayName', 'F_{sil} (used)'); hold on;
plot(age_all, fluxes_all(:,19), 'b--', 'LineWidth', 1, 'DisplayName', 'F_{sil} (BLAG feedback)');
plot(age_all, fluxes_all(:,10), 'r-', 'LineWidth', 2, 'DisplayName', 'F_{degas}');
xline(age_switch, 'r--', 'HandleVisibility', 'off');
set(gca, 'XDir', 'reverse');
xlabel('Age (Ma)'); ylabel('Flux (10^{18} mol/Myr)');
title('Silicate Weathering & Degassing'); legend('Location', 'best'); grid on;

subplot(2,2,3);
plot(age_all, p.T0 + DeltaT_all - 273, 'r-', 'LineWidth', 2); hold on;
yline(15, 'k--', 'Present');
xline(age_switch, 'r--', 'HandleVisibility', 'off');
set(gca, 'XDir', 'reverse');
xlabel('Age (Ma)'); ylabel('Temperature (°C)');
title('Global Mean Temperature (Eq. 31)'); grid on;

subplot(2,2,4);
plot(age_all, fluxes_all(:,21), 'b-', 'LineWidth', 2); hold on;
xline(age_switch, 'r--', 'HandleVisibility', 'off');
set(gca, 'XDir', 'reverse');
xlabel('Age (Ma)'); ylabel('pH');
title('Ocean pH (Eq. 60A)'); grid on;

sgtitle('BLAG Model: Spinup (570–424 Ma) + Li Forcing (424–300 Ma)', 'FontSize', 14);

% --- Figure 2: Phase 2 CO2 Comparison (424–300 Ma) ---
figure('Name', 'Phase 2: pCO2 Comparison', 'Position', [100 100 800 450]);
plot(age2, CO2_ppm2, 'b-', 'LineWidth', 2.5, 'DisplayName', 'BLAG (prescribed F_{sil})');
hold on;
scatter(Li_ages_plot, Li_pCO2_plot, 80, 'ro', 'filled', ...
    'DisplayName', 'Li model pCO_2 (power-law)');
for j = 1:length(Li_ages_plot)
    if ~isnan(Li_pCO2_p16_plt(j)) && ~isnan(Li_pCO2_p84_plt(j))
        plot([Li_ages_plot(j) Li_ages_plot(j)], ...
             [Li_pCO2_p16_plt(j) Li_pCO2_p84_plt(j)], 'r-', 'LineWidth', 1, ...
             'HandleVisibility', 'off');
    end
end
yline(280, 'k--', 'Preindustrial', 'HandleVisibility', 'off');
yline(CO2_424_ppm, 'g:', sprintf('%.0f ppm (spinup)', CO2_424_ppm), 'HandleVisibility', 'off');
set(gca, 'XDir', 'reverse', 'FontSize', 12);
xlabel('Age (Ma)', 'FontSize', 13); ylabel('pCO_2 (ppm)', 'FontSize', 13);
title(sprintf('Atmospheric CO_2 (%.0f–%.0f Ma) — BLAG Equations', age_switch, age_end), ...
    'FontSize', 14);
legend('Location', 'best', 'FontSize', 11);
grid on; xlim([age_end age_switch]);
ylim([0 max(8000, max(CO2_ppm2)*1.1)]);

% --- Figure 3: Phase 2 CO2 Budget ---
figure('Name', 'Phase 2: CO2 Budget', 'Position', [140 140 1400 500]);

subplot(1,3,1);
CO2_src  = 2*fluxes2(:,8) + fluxes2(:,9) + fluxes2(:,7);
CO2_sink_carb = 2*fluxes2(:,1) + fluxes2(:,2);
CO2_sink_sil  = 2*fluxes2(:,5);
plot(age2, CO2_src, 'r-', 'LineWidth', 2, 'DisplayName', 'CO_2 sources'); hold on;
plot(age2, CO2_sink_carb, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Carb weath sink');
plot(age2, CO2_sink_sil, 'b-', 'LineWidth', 2, 'DisplayName', 'Sil weath sink');
plot(age2, CO2_src - CO2_sink_carb - CO2_sink_sil, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'Net dCO_2/dt');
yline(0, 'k:'); set(gca, 'XDir', 'reverse');
xlabel('Age (Ma)'); ylabel('Flux (10^{18} mol/Myr)');
title('CO_2 Budget (Eq. 59H)'); legend('Location', 'best'); grid on;

subplot(1,3,2);
plot(age2, fluxes2(:,18), 'k-', 'LineWidth', 2, 'DisplayName', 'WI cumul.'); hold on;
Li_mt_data = age_switch - Li_ages;
WI_mask = (Li_mt_data >= 0) & (Li_mt_data <= t_end_phase2);
scatter(Li_ages(WI_mask), WI_cumul(WI_mask), 80, 'rs', 'filled', 'DisplayName', 'WI data');
yline(1, 'k:', 'Baseline');
set(gca, 'XDir', 'reverse');
xlabel('Age (Ma)'); ylabel('WI(t) / WI(424 Ma)');
title('Weathering Intensity Multiplier'); legend('Location', 'best'); grid on;

subplot(1,3,3);
plot(age2, fluxes2(:,5), 'b-', 'LineWidth', 2, 'DisplayName', 'F_{sil} (prescribed)'); hold on;
plot(age2, fluxes2(:,19), 'b--', 'LineWidth', 1, 'DisplayName', 'F_{sil} (feedback)');
plot(age2, fluxes2(:,10), 'r-', 'LineWidth', 2, 'DisplayName', 'F_{degas}');
set(gca, 'XDir', 'reverse');
xlabel('Age (Ma)'); ylabel('Flux (10^{18} mol/Myr)');
title('Weathering & Degassing'); legend('Location', 'best'); grid on;

sgtitle(sprintf('Phase 2 Details (%.0f–%.0f Ma)', age_switch, age_end), 'FontSize', 13);

% --- Figure 4: Ocean Chemistry ---
figure('Name', 'Ocean Chemistry', 'Position', [180 180 1400 500]);

subplot(1,4,1);
plot(age_all, y_all(:,6), 'b-', 'LineWidth', 1.5); hold on;
yline(M_Ca0, 'b--', 'Present'); xline(age_switch, 'r--');
set(gca, 'XDir', 'reverse');
xlabel('Age (Ma)'); ylabel('M_{Ca} (10^{18} mol)'); title('Ocean Ca'); grid on;

subplot(1,4,2);
plot(age_all, y_all(:,7), 'g-', 'LineWidth', 1.5); hold on;
yline(M_HCO3_0, 'g--', 'Present'); xline(age_switch, 'r--');
set(gca, 'XDir', 'reverse');
xlabel('Age (Ma)'); ylabel('M_{HCO3} (10^{18} mol)'); title('Ocean HCO_3^-'); grid on;

subplot(1,4,3);
plot(age_all, y_all(:,5), 'r-', 'LineWidth', 1.5); hold on;
yline(M_Mg0, 'r--', 'Present'); xline(age_switch, 'r--');
set(gca, 'XDir', 'reverse');
xlabel('Age (Ma)'); ylabel('M_{Mg} (10^{18} mol)'); title('Ocean Mg'); grid on;

subplot(1,4,4);
plot(age_all, y_all(:,1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Dolomite'); hold on;
plot(age_all, y_all(:,2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Calcite');
xline(age_switch, 'r--'); set(gca, 'XDir', 'reverse');
xlabel('Age (Ma)'); ylabel('10^{18} mol');
title('Carbonate Reservoirs'); legend('Location', 'best'); grid on;

sgtitle('Ocean Chemistry & Carbonate Reservoirs (570–300 Ma)', 'FontSize', 13);

% --- Figure 5: Correction Factors (Phase 2 only) ---
figure('Name', 'Correction Factors', 'Position', [220 220 1200 400]);
subplot(1,4,1);
plot(age2, fluxes2(:,12), 'b-', 'LineWidth', 2);
set(gca, 'XDir', 'reverse'); xlabel('Age (Ma)'); ylabel('f_A');
title('Land Area (Eq. 26)'); grid on;
subplot(1,4,2);
plot(age2, fluxes2(:,14), 'g-', 'LineWidth', 2);
set(gca, 'XDir', 'reverse'); xlabel('Age (Ma)'); ylabel('f_T');
title('Temp Factor (Eq. 28)'); grid on;
subplot(1,4,3);
plot(age2, fluxes2(:,15), 'r-', 'LineWidth', 2);
set(gca, 'XDir', 'reverse'); xlabel('Age (Ma)'); ylabel('f_B');
title('CO_2 Factor (Eq. 32)'); grid on;
subplot(1,4,4);
plot(age2, fluxes2(:,22), 'm-', 'LineWidth', 2);
set(gca, 'XDir', 'reverse'); xlabel('Age (Ma)'); ylabel('f_A \times f_T \times f_B');
title('Combined Modifier'); grid on;
sgtitle('BLAG Correction Factors — Phase 2', 'FontSize', 13);

%% =====================================================================
%  SECTION 11: SUMMARY TABLE
%  =====================================================================
fprintf('\n=== Phase 2 Results Summary ===\n');
fprintf('%6s  %8s  %8s  %6s  %8s  %8s  %8s  %8s  %8s  %6s\n', ...
    'Age', 'CO2xPAL', 'CO2ppm', 'dT(C)', 'F_sil', 'F_degas', 'F_prec', ...
    'Fw_D', 'Fw_C', 'WI');

ages_check = [424 420 415 410 405 400 395 390 385 380 375 370 365 360 ...
              355 350 345 340 335 330 325 320 315 310 305 300];
for j = 1:length(ages_check)
    [~, idx] = min(abs(age2 - ages_check(j)));
    fprintf('%6.0f  %8.2f  %8.0f  %+6.1f  %8.3f  %8.3f  %8.3f  %8.3f  %8.3f  %6.3f\n', ...
        age2(idx), CO2_xPAL2(idx), CO2_ppm2(idx), DeltaT2(idx), ...
        fluxes2(idx,5), fluxes2(idx,10), fluxes2(idx,7), ...
        fluxes2(idx,1), fluxes2(idx,2), fluxes2(idx,18));
end

%% =====================================================================
%  SECTION 12: SAVE RESULTS
%  =====================================================================
ages_out = (age_switch:-1:age_end)';
nout = length(ages_out);
results = table();
results.Age_Ma = ages_out;

t_out = age_switch - ages_out;
results.CO2_xPAL     = interp1(t2, CO2_xPAL2, t_out, 'linear');
results.CO2_ppm      = results.CO2_xPAL * 280;
results.DeltaT_K     = interp1(t2, DeltaT2, t_out, 'linear');
results.T_C          = p.T0 + results.DeltaT_K - 273;
results.F_sil_total  = interp1(t2, fluxes2(:,5),  t_out, 'linear');
results.F_degas      = interp1(t2, fluxes2(:,10), t_out, 'linear');
results.F_prec       = interp1(t2, fluxes2(:,7),  t_out, 'linear');
results.Fw_D         = interp1(t2, fluxes2(:,1),  t_out, 'linear');
results.Fw_C         = interp1(t2, fluxes2(:,2),  t_out, 'linear');
results.F_vsw_Mg     = interp1(t2, fluxes2(:,6),  t_out, 'linear');
results.WI_mult      = interp1(t2, fluxes2(:,18), t_out, 'linear');
results.F_sil_blag   = interp1(t2, fluxes2(:,19), t_out, 'linear');
results.pH           = interp1(t2, fluxes2(:,21), t_out, 'linear');
results.Dolomite     = interp1(t2, y2(:,1), t_out, 'linear');
results.Calcite      = interp1(t2, y2(:,2), t_out, 'linear');
results.M_Mg         = interp1(t2, y2(:,5), t_out, 'linear');
results.M_Ca         = interp1(t2, y2(:,6), t_out, 'linear');
results.M_HCO3       = interp1(t2, y2(:,7), t_out, 'linear');

outfile = 'BLAG_Fsil_perturbation_results.csv';
writetable(results, outfile);
fprintf('\nResults saved to %s (%d rows)\n', outfile, nout);
fprintf('\n=== Done ===\n');

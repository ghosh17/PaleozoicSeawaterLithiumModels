function [dydt, flux_out] = BLAG_odes_Fsil(t, y, p)
% BLAG_ODES_FSIL  ODE system for the BLAG model (Berner, Lasaga & Garrels 1983).
%
%   Supports two modes controlled by p.mode:
%     'feedback'   — Full BLAG with CO2-weathering feedback (Eqs. 33-35)
%     'prescribed' — Silicate weathering prescribed from Li isotope data
%
%   State vector y (all in 10^18 mol):
%     y(1) = D         Dolomite (Ca-equivalent moles)
%     y(2) = C         Calcite (CaCO3)
%     y(3) = S_CaSi    Ca-silicate minerals
%     y(4) = S_MgSi    Mg-silicate minerals
%     y(5) = M_Mg      Ocean dissolved Mg
%     y(6) = M_Ca      Ocean dissolved Ca
%     y(7) = M_HCO3    Ocean dissolved HCO3
%     y(8) = A_CO2     Atmospheric CO2 (10^18 mol)
%
%   Implements Eqs. 23-60A from BLAG 1983.
%   Reference: Berner, Lasaga & Garrels (1983), Am. J. Sci. 283, 641-683.

    %% Unpack state (enforce positivity)
    D      = max(y(1), 1e-10);
    C      = max(y(2), 1e-10);
    S_CaSi = max(y(3), 1e-10);
    S_MgSi = max(y(4), 1e-10);
    M_Mg   = max(y(5), 1e-10);
    M_Ca   = max(y(6), 1e-10);
    M_HCO3 = max(y(7), 1e-10);
    A_CO2  = max(y(8), 1e-10);

    if ~isfield(p, 'mode') || ~ismember(p.mode, {'feedback', 'prescribed'})
        error('BLAG:InvalidMode', 'p.mode must be ''feedback'' or ''prescribed''.');
    end

    %% ================================================================
    %  SPREADING RATE AND LAND AREA
    %  ================================================================
    fA  = max(p.interp_fA(t),  0.01);
    fSR = max(p.interp_fSR(t), 0.01);

    %% ================================================================
    %  TEMPERATURE FROM CO2 (Eq. 31)
    %  ================================================================
    rCO2   = A_CO2 / p.A_CO2_0;
    DeltaT = log(max(rCO2, 1e-6)) / p.greenhouse_coeff;
    DeltaT = max(min(DeltaT, 30), -30);

    %% ================================================================
    %  CORRECTION FACTORS (Eqs. 27, 28, 29, 32)
    %  ================================================================
    R_ratio      = max(1 + p.runoff_coeff * DeltaT, 0.01);   % Eq. 27
    C_HCO3_ratio = max(1 + p.bicarb_coeff * DeltaT, 0.01);  % Eq. 29
    fT           = C_HCO3_ratio * R_ratio;                     % Eq. 28
    fB           = max(1.0 + p.fB_coeff * log(max(rCO2, 1e-6)), 0.01);  % Eq. 32
    fw_combined  = fA * fT * fB;

    %% ================================================================
    %  CARBONATE WEATHERING (Eqs. 33-34)
    %  ================================================================
    Fw_D = p.kw_D0 * fw_combined * D;
    Fw_C = p.kw_C0 * fw_combined * C;

    %% ================================================================
    %  SILICATE WEATHERING
    %  ================================================================
    % BLAG feedback values (Eq. 35)
    Fw_CaSi_blag = p.kw_CaSi0 * fw_combined * S_CaSi;
    Fw_MgSi_blag = p.kw_MgSi0 * fw_combined * S_MgSi;
    F_sil_blag   = Fw_CaSi_blag + Fw_MgSi_blag;

    if strcmp(p.mode, 'prescribed')
        % Prescribed from Li isotope data
        if ~isfield(p, 'F_sil_baseline') || ~isfield(p, 'interp_WI_cumul')
            error('BLAG:MissingParam', ...
                'Prescribed mode requires p.F_sil_baseline and p.interp_WI_cumul.');
        end
        WI_mult     = max(p.interp_WI_cumul(t), 1e-6);
        F_sil_total = p.F_sil_baseline * WI_mult;
        Fw_CaSi     = F_sil_total * p.frac_CaSi;
        Fw_MgSi     = F_sil_total * p.frac_MgSi;
    else
        % Standard BLAG feedback (Eq. 35)
        WI_mult     = 1.0;
        Fw_CaSi     = Fw_CaSi_blag;
        Fw_MgSi     = Fw_MgSi_blag;
        F_sil_total = F_sil_blag;
    end

    %% ================================================================
    %  VOLCANIC-SEAWATER Mg EXCHANGE (Eqs. 39-40)
    %  ================================================================
    k_vsw    = p.k_vsw0 * fSR;
    F_vsw_Mg = k_vsw * M_Mg;

    %% ================================================================
    %  METAMORPHISM (Eqs. 41-44)
    %  ================================================================
    k_MD   = p.k_MD0 * fSR;
    k_MC   = p.k_MC0 * fSR;
    Fmet_D = k_MD * D;
    Fmet_C = k_MC * C;

    F_degas_equiv = 2*Fmet_D + Fmet_C;

    %% ================================================================
    %  CALCITE PRECIPITATION (Eqs. 45-47, 53A-53B)
    %  ================================================================
    prec_arg = M_Ca * M_HCO3^2 - p.K_eq * A_CO2;
    prec_arg = max(prec_arg, 0);
    F_prec   = p.k_prec0 * prec_arg;

    %% ================================================================
    %  DIAGNOSTICS
    %  ================================================================
    F_river_Ca = Fw_D + Fw_C + Fw_CaSi + F_vsw_Mg;
    pH_ocean   = 6.29 - log10(max(A_CO2, 1e-10) / max(M_HCO3, 1e-10));

    %% ================================================================
    %  ODE SYSTEM (Eqs. 59A-59H)
    %  ================================================================
    dydt = zeros(8, 1);

    dydt(1) = -(Fw_D + Fmet_D);                                          % 59A
    dydt(2) = -(Fw_C + Fmet_C) + F_prec;                                 % 59B
    dydt(3) = Fmet_C + Fmet_D - Fw_CaSi - F_vsw_Mg;                     % 59C
    dydt(4) = Fmet_D - Fw_MgSi + F_vsw_Mg;                              % 59D
    dydt(5) = Fw_D + Fw_MgSi - F_vsw_Mg;                                % 59E
    dydt(6) = Fw_D + Fw_C + Fw_CaSi + F_vsw_Mg - F_prec;               % 59F
    dydt(7) = 4*Fw_D + 2*Fw_C + 2*Fw_CaSi + 2*Fw_MgSi - 2*F_prec;     % 59G
    dydt(8) = 2*Fmet_D - 2*Fw_D ...                                      % 59H
            - 2*Fw_CaSi - 2*Fw_MgSi ...
            + Fmet_C - Fw_C ...
            + F_prec;

    %% ================================================================
    %  DIAGNOSTIC FLUX OUTPUT
    %  ================================================================
    if nargout > 1
        flux_out = [Fw_D, ...              %  1: dolomite weathering
                    Fw_C, ...              %  2: calcite weathering
                    Fw_CaSi, ...           %  3: Ca-silicate weathering
                    Fw_MgSi, ...           %  4: Mg-silicate weathering
                    F_sil_total, ...       %  5: total silicate weathering
                    F_vsw_Mg, ...          %  6: volcanic-seawater Mg exchange
                    F_prec, ...            %  7: calcite precipitation
                    Fmet_D, ...            %  8: dolomite metamorphism
                    Fmet_C, ...            %  9: calcite metamorphism
                    F_degas_equiv, ...     % 10: total metamorphic CO2
                    DeltaT, ...            % 11: temperature anomaly (K)
                    fA, fSR, ...           % 12-13: correction factors
                    fT, fB, ...            % 14-15: temperature & CO2 factors
                    R_ratio, ...           % 16: runoff ratio (Eq. 27)
                    C_HCO3_ratio, ...      % 17: HCO3 ratio (Eq. 29)
                    WI_mult, ...           % 18: WI cumulative multiplier
                    F_sil_blag, ...        % 19: BLAG feedback F_sil
                    F_river_Ca, ...        % 20: total river Ca
                    pH_ocean, ...          % 21: ocean pH (Eq. 60A)
                    fw_combined];          % 22: combined weathering modifier
    end
end

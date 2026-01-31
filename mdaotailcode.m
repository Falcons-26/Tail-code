clc; clear;

% ==================================================
% USER INPUTS
% ==================================================

disp('--- INPUT PARAMETERS (CONSISTENT UNITS) ---');
T = input('Enter applied torque T (N·mm): ');
M = input('Enter bending moment M (N·mm): ');
IY_MAX = input('Enter maximum allowable Iy (mm^4): ');

file_name = input('Enter output Excel file name (without .xlsx): ', 's');
if ~endsWith(lower(file_name), '.xlsx')
    file_name = strcat(file_name, '.xlsx');
end

% Beam length (for deflection)
L = 348.31;   % mm

% ==================================================
% MATERIAL PROPERTIES
% ==================================================

E = 3120.0;   % N/mm^2 (Young''s modulus for balsa)

density_map = containers.Map( ...
    [4.0, 6.0], ...
    [97.5, 93.75] ...
); % kg/m^3

% ==================================================
% GEOMETRIC OPTIONS
% ==================================================

b_values = [4.0, 6.0];   % flange thickness (mm)
c_values = [4.0, 6.0];   % web thickness (mm)

results = [];

row = 1;

% ==================================================
% PARAMETRIC STUDY
% ==================================================

for b = b_values
    for c = c_values

        rho_b = density_map(b);
        rho_c = density_map(c);

        for a = 10:19   % web height (mm)

            if (a + b) >= 20
                continue
            end

            for d = 10:40  % flange width (mm)

                % ----------------------------------
                % SECTION PROPERTIES
                % ----------------------------------

                Af = d * b;
                Aw = c * a;

                yf = b / 2.0;
                yw = b + a / 2.0;

                y_bar = (Af * yf + Aw * yw) / (Af + Aw);

                If = (d * b^3) / 12.0;
                Iw = (c * a^3) / 12.0;

                Iy = If + Af * (yf - y_bar)^2 + ...
                     Iw + Aw * (yw - y_bar)^2;

                % ----------------------------------
                % FILTER: MAX INERTIA
                % ----------------------------------

                if Iy > IY_MAX
                    continue
                end

                % ----------------------------------
                % BENDING STRESS
                % ----------------------------------

                y_max = max(abs(y_bar), abs((b + a) - y_bar));
                sigma_b = M * y_max / Iy;

                % ----------------------------------
                % TORSIONAL SHEAR STRESS
                % ----------------------------------

                J = (1/3) * (b * d^3 + c * a^3);
                tau = T / J;

                % ----------------------------------
                % PRINCIPAL STRESS
                % ----------------------------------

                sigma_max = ...
                    sigma_b / 2 + sqrt((sigma_b / 2)^2 + tau^2);

                % ----------------------------------
                % DEFLECTION
                % ----------------------------------

                delta = (M * L^2) / (12 * E * Iy);

                % ----------------------------------
                % MASS (per 1 m length)
                % ----------------------------------

                a_m = a / 1000;
                b_m = b / 1000;
                c_m = c / 1000;
                d_m = d / 1000;

                volume_flange = d_m * b_m;
                volume_web = a_m * c_m;

                mass_kg = rho_b * volume_flange + rho_c * volume_web;
                mass_g = mass_kg * 1000;

                % ----------------------------------
                % SAVE RESULTS
                % ----------------------------------

                results(row,:) = [ ...
                    a, b, c, d, Iy, ...
                    sigma_b, tau, sigma_max, ...
                    delta, mass_g ...
                ];
                row = row + 1;

            end
        end
    end
end

% ==================================================
% OUTPUT
% ==================================================

varNames = { ...
    'a_mm', 'b_mm', 'c_mm', 'd_mm', ...
    'Iy_mm4', ...
    'BendingStress_Nmm2', ...
    'ShearStress_Nmm2', ...
    'PrincipalStress_Nmm2', ...
    'Deflection_mm', ...
    'Mass_g' ...
};

T_results = array2table(results, 'VariableNames', varNames);

T_sorted = sortrows(T_results, 'Iy_mm4', 'descend');
T_final = T_sorted(1:min(150, height(T_sorted)), :);

save_path = 'C:\Users\atuly\Documents';
if ~exist(save_path, 'dir')
    mkdir(save_path);
end

file_path = fullfile(save_path, file_name);
writetable(T_final, file_path);

disp('Simulation complete.');
disp(['Top 150 sections (highest Iy first) saved to:']);
disp(file_path);

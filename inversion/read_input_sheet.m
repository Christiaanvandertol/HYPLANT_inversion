function tab = read_input_sheet(pathInput)

    if nargin == 0
        pathInput = fullfile('../../input');
    end
    tab = readtable(fullfile(pathInput, 'inversion_input_data.csv'));

    expected_cols = {'tune', 'variable', 'value', 'upper', 'lower', 'uncertainty', 'description', 'units'};
    colnames = tab.Properties.VariableNames;
    assert(all(ismember(expected_cols, colnames)), ...
           ['wrong column names in input data table; expected: ' sprintf('%s, ', expected_cols{:})])

    tab.include = (tab.tune == 1);

    %% to avoid fitting LIDFa without LIDFb
    i_lidfa = strcmp(tab.variable, 'LIDFa');
    i_lidfb = strcmp(tab.variable, 'LIDFb');
    if tab.include(i_lidfa) || tab.include(i_lidfb)
        tab.include(i_lidfa) = 1;
        tab.include(i_lidfb) = 1;
    end

    %% smc to 1e-2
    i_smc = strcmp(tab.variable, 'SMC');
    tab.value(i_smc) = tab.value(i_smc) * 1e-2;
    tab.lower(i_smc) = tab.lower(i_smc) * 1e-2;
    tab.upper(i_smc) = tab.upper(i_smc) * 1e-2;
    tab.uncertainty(i_smc) = tab.uncertainty(i_smc) * 1e-2;

    %% sif include in fit all 4 components
%     i_sif = contains(tab.variable, 'SIF');  % >= 2016a
    i_sif = ~cellfun(@isempty, strfind(tab.variable, 'SIF')); % <= 2015b
    if any(tab.include(i_sif))
        tab.include(i_sif) = 1;
    end

end
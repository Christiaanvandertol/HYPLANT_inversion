function tab = readInputSheetOctave(pathInput)
    tab = struct();
    % pathInput = fullfile('input', 'inversion_input_data.csv');

    fid = fopen([pathInput 'inversion_input_data.csv']);
    fgetl(fid);
    a = textscan(fid,'%d%s%f%f%f%f%s%s','Delimiter', ',');
    fclose(fid);

    tab.tune = a{1};
    tab.variable = a{2};
    tab.value = a{3};
    tab.lower = a{4};
    tab.upper = a{5};
    tab.uncertainty = a{6};
    tab.description = a{7};
    tab.unit = a{8};
    tab.include = (tab.tune == 1);

%keyboard

  %  [tab.tune,tab.variable,tab.value,tab.lower, ...
  %  talb.upper,tab.uncertainty,tab.description.tab.unit] = ...
  %  fscanf(fid,'%d,%s,%f,%f,%f,%f%s,%s\r')
  %  fclose(fid)

   % data = importdata([pathInput 'inversion_input_data.csv']);
  %  data = csvread([pathInput 'inversion_input_data.csv']);
  %  keyboardr
 %  num = data.data;
 %   txt = data.textdata;

  %  [r, ~] = size(num);
  %  tab.tune = num(:, 1);
  %  tab.include = (tab.tune == 1);
  %  tab.variable = txt(2: (r + 1));
  %  tab.value = num(:, 3);
  %  tab.lower = num(:, 4);
  %  tab.upper = num(:, 5);
  %  tab.uncertainty = num(:, 6);
end




    %fid = fopen([pathInput 'inversion_input_data.csv']);
    %fgetl(fid)
    %[A,B,C,D,E,F,G,H] = ...
    %fscanf(fid,'%d,%s,%f,%f,%f,%f,%s,%s\r')
    %fclose(fid)

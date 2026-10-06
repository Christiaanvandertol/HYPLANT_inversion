function [wl, spec, tts] = parseFlox(wl, numUp, numDown, txt, date, lon, lat)


    %% compute reflectance

    num = numUp ./ numDown;
    num = num(:, 2:end);
    
    i = (wl > 400) & (wl < 900);
    wl = wl(i);
    spec = num(i, :);


    %% compute SZA
    % expected date format yymmdd
    % excpected header format XHH_MM_SS

    y = str2num(['20' date(1:2)]);
    soy = datenum(y, 1, 1);

    dates = cellfun(@(x) datenum([date x(1:9)], 'yymmddXHH_MM_SS'), txt(2:end));

    julDoy = dates - soy;
    doy = floor(julDoy);
    t = (julDoy - doy) .* 24;
    % lon = 5;
    % lat = 52;
    
    sza = calczenithangle(doy, t, 0, 0, lon, lat);
    tts = sza / pi * 180;

end

function [v_catalog, star_names, ra_deg, dec_deg, mag_catalog] = catalog_matching(filename)

catalog = readtable(filename, 'VariableNamingRule', 'preserve');

% --- RA/Dec conversion ---
ra_deg = raStrToDeg(catalog.RightAscension);
dec_deg = decStrToDeg(catalog.Declination);

ra = deg2rad(ra_deg);
dec = deg2rad(dec_deg);

x = cos(dec).*cos(ra);
y = cos(dec).*sin(ra);
z = sin(dec);

v_catalog = [x y z];

% --- names ---
star_names = catalog.Name;

% --- magnitude (NEW) ---
if ismember('Apparent Magnitude', catalog.Properties.VariableNames)
    mag_catalog = catalog.("Apparent Magnitude");
else
    error('Magnitude column not found in catalog');
end

end


function deg = raStrToDeg(raStrings)

n = numel(raStrings);
deg = zeros(n,1);

for i = 1:n
    s = raStrings{i};
    s = regexprep(s, char(160), ' ');
    s = strrep(s,'h',' ');
    s = strrep(s,'m',' ');
    s = strrep(s,'s','');

    parts = sscanf(s,'%f %f %f');

    if numel(parts) == 3
        deg(i) = 15*(parts(1) + parts(2)/60 + parts(3)/3600);
    else
        deg(i) = NaN;
    end
end

end


function deg = decStrToDeg(decStrings)

n = numel(decStrings);
deg = zeros(n,1);

for i = 1:n
    s = decStrings{i};
    s = regexprep(s, char(160), ' ');
    s = strrep(s,'°',' ');
    s = strrep(s,'′',' ');
    s = strrep(s,'″','');

    sign = 1;
    if contains(s,'-')
        sign = -1;
    end

    s = strrep(s,'+','');
    parts = sscanf(s,'%f %f %f');

    if numel(parts) == 3
        deg(i) = sign*(abs(parts(1)) + parts(2)/60 + parts(3)/3600);
    else
        deg(i) = NaN;
    end
end

end
function outstrand = comp_str(in_strd)

lth = length(in_strd);
for i = 1:lth
    outstrand(i) = 'N';
end

for i = 1:lth
    if in_strd(i) == 'A'
    outstrand(lth-i+1) = 'T';
    end
    if in_strd(i) == 'a'
    outstrand(lth-i+1) = 't';
    end
    if in_strd(i) == 'T' || in_strd(i) == 'U'
    outstrand(lth-i+1) = 'A';
    end
    if in_strd(i) == 't' || in_strd(i) == 'u'
    outstrand(lth-i+1) = 'a';
    end
    if in_strd(i) == 'C'
    outstrand(lth-i+1) = 'G';
    end
    if in_strd(i) == 'c'
    outstrand(lth-i+1) = 'g';
    end
    if in_strd(i) == 'G'
    outstrand(lth-i+1) = 'C';
    end
    if in_strd(i) == 'g'
    outstrand(lth-i+1) = 'c';
    end
    if in_strd(i) == ' '
    outstrand(lth-i+1) = ' ';
    end
end

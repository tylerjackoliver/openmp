strarr = zeros(6, 2);
cnt = 1;
chunks = [1 2 3 6 12 24];
    
for j = 1:6

        a = table2array((import_file(sprintf("BEST_%d/output.out", chunks(j)))));
        a = cellstr(a);
        str_l1 = a{1}(40:end);
        str_l2 = a{2}(40:end);
        strarr(cnt, 1) = str2num(str_l1);
        strarr(cnt, 2) = str2num(str_l2);
        cnt = cnt + 1;
        
end
        
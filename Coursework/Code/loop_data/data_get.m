strarr = zeros(23, 2);
cnt = 1;
for i = 1:3
    
    for j = 1:7

        a = table2array((import_file(sprintf("%s%d/output.out", folder_names(i), chunk_size(j)))));
        a = cellstr(a);
        str_l1 = a{1}(40:end);
        str_l2 = a{2}(40:end);
        strarr(cnt, 1) = str2num(str_l1);
        strarr(cnt, 2) = str2num(str_l2);
        cnt = cnt + 1;
        
    end
    
end
        
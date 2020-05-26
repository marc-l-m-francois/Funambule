function disptitle(string)
    le = 72-2-size(string,2);
    disp(' ');
    disp([char(61*ones(1,floor(le/2))),' ',string,' ',char(61*ones(1,ceil(le/2)))]);
    disp(' ');
end

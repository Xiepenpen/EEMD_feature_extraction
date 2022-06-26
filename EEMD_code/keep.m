%This function is to keep a few decimal places
function aaa = keep(a,n)

format long

b = a*10^n;
aa = floor(b);
need_num = aa/10^n;
need_str = num2str(need_num);


aaa = str2num(need_str);
    
end

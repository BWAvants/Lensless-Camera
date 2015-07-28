function y = fhandle_mode(x,mode, f_h, ft_h) 

if sum(mode == 1) || strcmp(mode,'notransp');
    y = f_h(x);
elseif sum(mode == 0) 
    y = ft_h(f_h(x));
else
    y = ft_h(x);
end
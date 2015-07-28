function flatCamera1
resp = inputdlg('Enter "n" dimension for Hadamard patterns','flatCamera1');
if isempty(resp)
    return;
end
Hn = str2double(resp{1});
if isnan(Hn)
    return;
end
Hrow = fwht(eye(n)*n,n,'hadamard');
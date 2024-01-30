
function [logical]=check_psd(A)
d = eig(A);
logical= all(d>0);
end
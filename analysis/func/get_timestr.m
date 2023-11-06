function [yrstr, dystr, hrstr] = get_timestr(nday, yr_s)
%Get strings of YEAR, DAY and HR for a given 'nday' (double) and 'yr_s' 
% (double; year at the beginning).
% 
% E.g., nday=270.25 and yr_s = 21 gives
%    yrstr = '0021', dystr = '270', hrstr = '06'
% 

yr2dy = 365;
dy2hr = 24;

nyr = floor((nday - 1)/yr2dy) + yr_s;
yrstr = num2str(nyr, '%4.4i');
dystr = num2str(floor(nday - (nyr-yr_s)*yr2dy), '%3.3i');
hrstr = num2str(round(mod(nday - yr2dy, 1)*dy2hr), '%2.2i');
    


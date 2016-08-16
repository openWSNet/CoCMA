function output=wi(pav,a,b,rs)
pref=a*(2*rs)^(-b);
y=(pav/pref)^(1/-b);
output=2/pi*(acos(y)-y*sqrt(1-y*y));
end
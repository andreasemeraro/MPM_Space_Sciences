% ----------------------------------------------------------------------
function density = atmosphere2(z)
%
% Mars Atmosphere heihgt scales throgh https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016JA023641
% Mars Atmosphere densities trhough https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/1998JE900041
% ---------------------------------------------------------------------
%...Geometric altitudes (km):
h = ...
[ 10 20 30 40 50 60 70 ...
80 90 100 110 120 130 140 ];

%...Corresponding densities (kg/m^3) from USSA76:
r = ...
[7.3e-3 2.6e-3 9.7e-4 3.4e-4 1.15e-4 4.06e-5 8.93e-6 ...
2.17e-6 3.68e-7 7.03e-8 1.66e-8 3.86e-9 7.0e-10 1.8e-10 ];
%...Scale heights (km):
H1=10.2 + 0.1*z; %we are not considering the solar influence (see web link below)
%H2 = 100.0;
%...Handle altitudes outside of the range:
if z > 140
z = 140;
elseif z < 0
z = 0;
end
%...Determine the interpolation interval:
for j = 1:13
if z >= h(j) && z < h(j+1)
i = j;
end
end
if z == 140
i = 13;
end
%...Exponential interpolation:
%if i>=1 && i<13
density = r(i)*exp(-(z - h(i))/H1);
% end 
% 
% if i==13
%     density = r(i)*exp(-(z-h(i))/H2);
% end
end%atmopshere
%-----------------------------------------------------------------------
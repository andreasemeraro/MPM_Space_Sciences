%Calcolo deltaV_drag
function DeltaV_Drag

%Questo script calcola il dV dovuto al drag atmosferico, come indicato in 
%eq 25 dell'articolo https://doi.org/10.1016/j.actaastro.2021.07.002

%ovviamente può essere riadattato per la Terra et al.
mu=42828; %mu Mars (kn3/s2)
a0=300;
p0=300;
rMars=3388;
an=(a0+p0+rMars)/2
p1=250;
a1=250;
af=(a1+p1+rMars)/2
delta_v = sqrt(2*mu/af - 2*mu/(an + af)) - sqrt(mu/af) + sqrt(mu/an) - sqrt(2*mu/an - 2*mu/(an + af)) %suppongo in km/s
end
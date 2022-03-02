function [a, ecc, inc, w, nu, RAAN] = Carts2Kepler(x, y, z, Vx, Vy, Vz)
%
%		    DESCIPTION:      This function is created to convert the classical 
%               		     orbital elements to cartesian position and velocity 
%		       		         parameters of any satellite orbit in the geocentric-equatorial
%                            reference system.
% 			INPUT:
%			Position Components: 		
% 			[X Y Z]...(Km)
%
%			Velocity Components:
% 			[Vx Vy Vz]...(Km/s) 
%
% 			OUTPUT:
% 			alt:    Altitude.....................(Km)							
% 			ecc:    Eccentricity											    
% 			inc:	Inclination..................(rad)							
% 			w:	    Argument of perigee..........(rad)	
% 			nu:	    Satellite position...........(rad)							
% 			RAAN:	Right Asc. of Ascending Node.(rad)							
%
%---------------------------- Constants ----------------------------------------------%
mu_earth = 3.986 * 10^5; % Earth Gravitational Constant (km3/s2)
re = 6378.1; % earth radius (km)

%------h vector, eccenticity vector and ascending node vector(n)--------------------

h = cross([x , y, z],[Vx, Vy, Vz]);

ecc = cross([Vx, Vy, Vz], h)/mu_earth - [x , y, z]/norm([x , y, z]);

n = [ -h(2), h(1), 0]';

%-----nu, inclination, RAAN, w, ecc(scalar)---------------------------------------------

if dot([x , y, z],[Vx, Vy, Vz])>=0
     
   nu = acos(dot(ecc,[x , y, z])/(norm(ecc)*norm([x , y, z])));
   
else
   
    nu = 2*pi - acos(dot(ecc,[x , y, z])/(norm(ecc)*norm([x , y, z])));
end

inc = acos(h(3)/norm(h));

if n(2)>=0
    RAAN = acos(n(1)/norm(n));
else
     RAAN =2*pi - acos(n(1)/norm(n));
end

if ecc(3)>=0
    w = acos(dot(n,ecc)/(norm(n)*norm(ecc)));
else
    w =2*pi - acos(dot(n,ecc)/(norm(n)*norm(ecc)));
end

ecc = norm(ecc);

a = 1/(2/norm([x, y, z])- (norm([Vx, Vy, Vz]))^2/mu_earth);

end


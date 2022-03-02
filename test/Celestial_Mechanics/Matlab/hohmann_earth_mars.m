%This script would like to perform an estimate of
%velocities occuring in the Hohmann transfer orbit equation
%from Earth to Mars
m=10^8; %kg initial mass of spacecraft
I_sp=400; %s specific impulse
g_0=9.81*10^-3; %m/s^2
muSun=132712440018; %gravitational constant of Sun (km^3/s^2)
R1=1.5*10^8; %Earth average distance from Sun (km)
R2=2.28*10^8; %Mars average distance from Sun (km)

%The circular orbital speed of planet 1 relative
%to the sun

V1=sqrt(muSun/R1);

%The circular orbital speed of planet 2 relative
%to the sun
V2=sqrt(muSun/R2);
formatSpec='La velocità eliocentrica media della Terra è %4.2f km/s. \n';
fprintf(formatSpec,V1);
formatSpec1='La velocità eliocentrica media di Marte è %4.2f km/s. \n';
fprintf(formatSpec1,V2);

%Hohmann transfer orbit velocities at departure and arrival
Vdeparture=sqrt((2*muSun*R2)/(R1*(R1+R2)));
formatSpec2=['La velocità eliocentrica dello spacecraft alla partenza' ...
     ' è %4.2f km/s. \n'];
fprintf(formatSpec2,Vdeparture);
dV1=Vdeparture-V1;
formatSpec3='L eccesso di velocità alla partenza è %4.2f km/s. \n';
fprintf(formatSpec3,dV1);
dV2=V2*(1-sqrt(2*R1/(R1+R2)));
formatSpec4='L eccesso di velocità all arrivo è %4.2f km/s. \n';
fprintf(formatSpec4,dV2);

%Total dV for a Hohmann transfer orbit
dVtot=dV1+dV2;
formatSpec5='Il dV totale richiesto è %4.2f km/s. \n';
fprintf(formatSpec5,dVtot);

%An estimate of the mass of the propellant used
dm=m*(1-exp(-dVtot/(g_0*I_sp)));
formatSpec6='La massa totale di propellente usato è %2.2f kg. \n';
fprintf(formatSpec6,dm);
formatSpec7=['La percentuale di massa di propellente vs massa iniziale' ...
    ' è del %4.2f %%. \n'];
fprintf(formatSpec7,dm/m*100);

%Give an estimate of time of flight

T=((2*((R1+R2)/2)^(3/2)/(sqrt(muSun)))/86400)/30; %s
formatSpec8=['Il tempo di volo dalla Terra a Marte è %4.2f mesi. \n'];
fprintf(formatSpec8,T);

import math


def app_pos(dt, phi, ksi, tz):
    
    b=2.-int(dt[0]/100.)+int(int(dt[0]/100.)/4.)        
    if dt[1]==1 or dt[1]==2:
        JD= int(365.25*(dt[0]+4715.)) + int(30.6001*(dt[1]+13.)) + dt[2] + b + dt[3]/24.+dt[4]/1440.+dt[5]/86400. -tz/24. - 1524.5

    else:
        JD= int(365.25*(dt[0]+4716.)) + int(30.6001*(dt[1]+1.)) + dt[2] + b + dt[3]/24.+dt[4]/1440.+dt[5]/86400. -tz/24. - 1524.5

    T=(JD-2451545.0)/36525.
        
    

#Mean longitude of the Sun
    L0=(280.46645+36000.76983*T+0.0003032*T**2)%360

#Eccentricity of the Earth's orbit
    e=0.016708617-0.000042037*T-0.0000001236*T**2

#Mean anomaly of the Sun
    M=357.52910+35999.05030*T-0.0001559*T**2-0.00000048*T**3

#Equation of center of the Sun
    C=(1.914600-0.004817*T-0.000014*T**2)*math.sin(math.radians(M)) +(0.019993 - 0.000101*T)*math.sin(2*math.radians(M)) + 0.000290*math.sin(3*math.radians(M))

#True longitude of the Sun
    theta=L0 + C

#Nutation and aberration correction of theta
    omega= 125.04452 - 1934.136261*T + 0.0020708*T**2 + T**3/450000.

#Apparent longitude of the Sun
    lambda1 = theta - 0.00569-0.00478*math.sin(math.radians(omega))

#Mean Obliquity of the ecliptic
    epsilon = 23+(26./60.)+(21.448/3600) - (46.8150/3600)*T - (0.00059/3600)*T**2+(0.001813/3600)*T**3

#Mean Latitude of the Moon      
    Lp=(218.3165 + 481267.8813*T)%360

#Nutiation in the Obliquity of the ecliptic
    Deltaepsilon = 9.20*math.cos(math.radians(omega)) + 0.57*math.cos(2*math.radians(L0)) + 0.1*math.cos(2*math.radians(Lp)) - 0.09*math.cos(2*math.radians(omega))

#Obliquity of the ecliptic corrected
    epsilon1=epsilon+Deltaepsilon/3600.

#Sun's Right ascension
    alpha=math.degrees(math.atan2(math.cos(math.radians(epsilon1))*math.sin(math.radians(lambda1)), math.cos(math.radians(lambda1))))

    vary= (math.tan(math.radians(epsilon1)/2))**2

#Equation of time
    rc=4*math.degrees((vary*math.sin(2*math.radians(L0)) - 2*e*math.sin(math.radians(M)) + 4*e*vary*math.sin(math.radians(M))*math.cos(2*math.radians(L0)) - 0.5*(vary)**2*math.sin(4*math.radians(L0)) - 1.25*e**2*math.sin(2*math.radians
        (M))))

#True solar time
    tst=(dt[3]*60+dt[4]+dt[5]/60.+rc +4*ksi -60*tz)%1440

#Hour angle
    H=tst/4.-(abs(tst)/tst)*180

#Declination of the Sun's
    delta = math.degrees(math.asin(math.sin(math.radians(epsilon1))*math.sin(math.radians(lambda1))))

#Solar zenith angle
    S=math.degrees(math.acos(math.sin(math.radians(phi))*math.sin(math.radians(delta))+math.cos(math.radians(phi))*math.cos(math.radians(delta))*math.cos(math.radians(H))))

#Nutiation in longitude
    Deltapsi=-17.2*math.sin(math.radians(omega))-1.32*math.sin(2*math.radians(L0))-0.23*math.sin(2*math.radians(Lp))+0.21*math.sin(2*math.radians(omega))

        

    
    Elevation=90-S
        

#Azymut 
    A=(math.degrees(math.acos(((math.sin(math.radians(phi))*math.cos(math.radians(S)))-math.sin(math.radians(delta)))/(math.cos(math.radians(phi))*math.sin(math.radians(S))))))%360*(abs(H)/H)

#Sun rise, set
    Srss=math.acos(-math.tan(math.radians(phi))*math.tan(math.radians(delta)))
#Sun rise, set corrected
    Srss_c=math.acos(math.cos(math.radians(90.833))/(math.cos(math.radians(phi))*math.cos(math.radians(delta)))-math.tan(math.radians(phi))*math.tan(math.radians(delta)))
#Solar noon
    sn=(720.-4*ksi-rc+60*tz)

    sr=math.degrees(Srss)
    ss=abs(math.degrees(Srss))
        

    return A, Elevation, sn, sr, ss

        


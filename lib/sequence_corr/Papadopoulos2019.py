import math

def Papadopoulos2019(T_ms,T_as):
    # applicable period range
    Trange = [0.01,5]    
    Tmax = max(T_ms,T_as)
    Tmin = min(T_ms,T_as)
    if Tmax > Trange[1] or Tmin < Trange[0]:
        return 0.0
    
    c1 = 1-math.cos((math.pi/3.55)*min(2.84,math.log(max(Tmax,0.075)/0.075)))
    
    if Tmax == Tmin:
        rho = 0.6+0.06*math.cos(math.pi+math.pi*(min(max(math.log(Tmax),-3.68),-1.38)+4.6)/1.76)
        
    elif Tmin>0.075:
        rho_not_Tmax = 0.6+0.06*math.cos(math.pi+math.pi*(min(max(math.log(Tmax),-3.68),-1.38)+4.6)/1.76)
        rho = rho_not_Tmax * (1-0.25*c1+0.25*c1*math.cos(math.pi*math.log(Tmin/Tmax)/math.log(Tmax/0.075))) 
    else:
        rho_not_Tmax = 0.6+0.06*math.cos(math.pi+math.pi*(min(max(math.log(Tmax),-3.68),-1.38)+4.6)/1.76)
        rho = rho_not_Tmax * (1-0.38*c1+0.12*c1*math.cos(math.pi + math.pi * math.log(max(Tmin,0.02)/0.075)/math.log(max(Tmax/0.075,3.75))))
    return rho
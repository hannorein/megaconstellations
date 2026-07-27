import sys
sys.path.insert(0,"./model")
import mega as m
import numpy as np
import struct

def TLE(n, inc, Omega, M):
    TLE = """STARLINK-FUTURE
1 49460U 21104BE  21337.58334491  .00068382  00000+0  48921-6 0  9997
2 49460 %8.4f %8.4f 0000000   0.0000 %8.4f %11.8f  4453
""" %(inc,Omega,M, n)
    return TLE
di = {"Starlink":"starlinkfuture","OneWeb":"onewebfuture", "StarNet/GW":"starnetfuture","Kuiper":"kuiperfuture", "SXODC":"sxodc", "Sunrise":"sunrise"}
di = {"SXODC":"sxodc", "Sunrise":"sunrise"}

def TLEbinary(n, inc, Omega, M):
    numbers = [
            n,
            inc,
            Omega,
            0.0,
            M,
            0.0,
            0.0
            ]
    format_string = f'{len(numbers)}f' 
    binary_data = struct.pack(format_string, *numbers)
    return binary_data

for k in di:
    count = 0
    count2 = 0
    name = di[k]
    ICs = m.constellations_all[k]
    with open(name+".txt", "wb") as f:
        for IC in ICs:
            nplanes = IC["NPLANES"]
            nsat = IC["SATPP"]
            count2 += nsat * nplanes
            a = IC["ALT"]*1000.0 + m.REarth ## in m
            GM = 2.9755363e+24 ## m^3/day^2
            n = np.sqrt(GM/(a*a*a))/np.pi/2.
            Omegas = np.linspace(0.,2.*np.pi,nplanes,endpoint=False)
            inc = IC["INC"]
            incsynch = np.arccos(-np.power(a/12352000,7./2.))/np.pi*180
            issynch = False
            if (np.abs(inc - incsynch)) <10:
                inc = incsynch
                issynch = True
            #print("ConstellationPlane(alt: %.1f, inc: %.1f, nplanes: %d, satpp: %d)," %(IC["ALT"], IC["INC"], IC["NPLANES"], IC["SATPP"]))

            for i, Omega in enumerate(Omegas):
                Ms = np.linspace(0.,2.*np.pi,nsat) + 2.*np.pi/nsat*0.25*np.random.normal(size=nsat)
                Omega = np.fmod(Omega/np.pi*180.+360.,360.)
                for M in Ms:
                    M = np.fmod(M/np.pi*180.+360.,360.)
                    _Omega = Omega
                    if issynch:
                        _Omega += 1.0*np.random.normal()
                    t = TLEbinary(n=n, inc=inc, Omega=_Omega, M=M)
                    count +=1
                    f.write(t)
    print(k, count, count2)

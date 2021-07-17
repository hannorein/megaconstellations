import KeplerTools as KT
import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D

####################
# Written by Aaron Boley, June 2021
# Edited by Sam Lawler, June 2021
# ####################

twopi=np.pi*2
MEarth = 5.97e24
REarth = 6378.135e3
REkm = 6378.135

aukm=1.496e8
au=aukm*1e3


#FIGNAME="JS_hor30"
#FIGNAME="EQ_hor30"
#FIGNAME="DS_hor30"
#FIGNAME="JS_hor0"
#FIGNAME="EQ_hor0"
FIGNAME="DS_hor0"

#TITLE1="June Solstice, >30 Degrees Above Horizon"
#TITLE1="Equinox, >30 Degrees Above Horizon"
#TITLE1="December Solstice, >30 Degrees Above Horizon"
#TITLE1="June Solstice, Above Horizon"
#TITLE1="Equinox, Above Horizon"
TITLE1="December Solstice, Above Horizon"
#tilt=23.4 * np.pi/180
#tilt=0 * np.pi/180
tilt=-23.4 * np.pi/180

NAVG=1

#ELMIN=30
ELMIN=0

#place sun on negative x axis, meaning we need to reverse the sign of the tilt
tilt*=-1


# first three VLEO
# next 

ICs=[]

ICs.append({'NPLANES':7178,'SATPP':1,'INC':30,'ALT':328})
ICs.append({'NPLANES':7178,'SATPP':1,'INC':40,'ALT':334})
ICs.append({'NPLANES':7178,'SATPP':1,'INC':53,'ALT':345})
ICs.append({'NPLANES':40,'SATPP':50,'INC':96.9,'ALT':360})
ICs.append({'NPLANES':1998,'SATPP':1,'INC':75,'ALT':373})
ICs.append({'NPLANES':4000,'SATPP':1,'INC':53,'ALT':499})
ICs.append({'NPLANES':12,'SATPP':12,'INC':148,'ALT':604})
ICs.append({'NPLANES':18,'SATPP':18,'INC':115.7,'ALT':614})

ICs.append({'NPLANES':2547,'SATPP':1,'INC':53,'ALT':345.6})
ICs.append({'NPLANES':2478,'SATPP':1,'INC':48,'ALT':340.8})
ICs.append({'NPLANES':2493,'SATPP':1,'INC':42,'ALT':335.9})
ICs.append({'NPLANES':32,'SATPP':50,'INC':53,'ALT':550})
ICs.append({'NPLANES':72,'SATPP':22,'INC':53.2,'ALT':540})
ICs.append({'NPLANES':36,'SATPP':20,'INC':70,'ALT':570})
ICs.append({'NPLANES':6,'SATPP':58,'INC':97.6,'ALT':560})
ICs.append({'NPLANES':4,'SATPP':43,'INC':97.6,'ALT':560.1})
ICs.append({'NPLANES':18,'SATPP':40,'INC':87.9,'ALT':1200})
ICs.append({'NPLANES':36,'SATPP':49,'INC':87.9,'ALT':1201})
ICs.append({'NPLANES':32,'SATPP':72,'INC':40,'ALT':1202})
ICs.append({'NPLANES':32,'SATPP':72,'INC':55,'ALT':1203})

ICs.append({'NPLANES':16,'SATPP':30,'INC':85,'ALT':590})
ICs.append({'NPLANES':40,'SATPP':50,'INC':50,'ALT':600})
ICs.append({'NPLANES':60,'SATPP':60,'INC':55,'ALT':508})
ICs.append({'NPLANES':48,'SATPP':36,'INC':30,'ALT':1145})
ICs.append({'NPLANES':48,'SATPP':36,'INC':40,'ALT':1145})
ICs.append({'NPLANES':48,'SATPP':36,'INC':50,'ALT':1145})
ICs.append({'NPLANES':48,'SATPP':36,'INC':60,'ALT':1145})

ICs.append({'NPLANES':34,'SATPP':34,'INC':51.9,'ALT':630})
ICs.append({'NPLANES':36,'SATPP':36,'INC':42,'ALT':610})
ICs.append({'NPLANES':28,'SATPP':28,'INC':33,'ALT':509})

#lats = [-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60]
lats = np.linspace(-70,70,141)
print(lats)

lats = np.array(lats)*twopi/360.
hours = np.linspace(-12,12,48)
print(hours)
hang = hours*twopi/24.
print(hang)

horizonlimit=(90-ELMIN)*twopi/360.



xobs0 = REarth*1
yobs0 = 0.
zobs0 = 0.

xl=[]
yl=[]
zl=[]
xa=[]
ya=[]
za=[]

count = np.zeros( (len(lats),len(hang)))

def RotateZ(x,y,z,a):
    xp = x*np.cos(a)-y*np.sin(a)
    yp = x*np.sin(a)+y*np.cos(a)
    return xp,yp,z

def RotateY(x,y,z,a):
    xp=x*np.cos(a)+z*np.sin(a)
    zp=-x*np.sin(a)+z*np.cos(a)
    return xp,y,zp


for ic in ICs:
   nplanes=ic['NPLANES']
   nsat=ic['SATPP']
   ntot=nplanes*nsat
   print(nsat)
   sma = (ic['ALT']+REkm)/aukm

   if nsat==1:
     dplane = twopi/ntot
     ma = np.linspace(0,twopi,ntot,endpoint=False)
     colat = np.pi*0.5-ic['INC']*np.pi/180
     angOpt = np.sqrt( (np.cos(colat)-np.cos(colat+np.pi*0.5))*twopi/ntot)

     mashuf=np.zeros(ntot)
  
     nmix=int(np.sqrt(ntot))
  
     iref=0
     iskip=0
     for i in range(ntot):
       imix = iref+iskip*nmix
       if imix>ntot-1:
          iref+=1
          iskip=0
          imix=iref*1
       mashuf[i] = ma[imix]*1
       iskip+=1

     isum=0
     for ilocal in range(nsat):
        for iplane in range(nplanes):
            Omega = iplane*dplane
            x0,y0,z0,vx,vy,vz = KT.getXYZVVV(mashuf[isum],sma,0,0,Omega,ic['INC']*np.pi/180.,m0=MEarth,m1=0.)
            x0*=au
            y0*=au
            z0*=au

            x,y,z= RotateY(x0,y0,z0,tilt)
            xa.append(x)
            ya.append(y)
            za.append(z)
            OK=1
            if np.sqrt(y*y+z*z)<REarth:
              if x > 0: OK=0
            if OK:
              xl.append(x)
              yl.append(y)
              zl.append(z)
            isum+=1
   else:
     dplane=twopi/nplanes
     mashuf = np.linspace(0,twopi,nsat,endpoint=False)
     print("Simple approach")
     for ilocal in range(nsat):
          for iplane in range(nplanes):
              Omega = iplane*dplane
              print(ilocal,iplane,Omega,nsat,nplanes)
              x0,y0,z0,vx,vy,vz = KT.getXYZVVV(mashuf[ilocal]+iplane*dplane/nplanes,sma,0,0,Omega,ic['INC']*np.pi/180.,m0=MEarth,m1=0.)
              x0*=au
              y0*=au
              z0*=au

              x,y,z= RotateY(x0,y0,z0,tilt)
              xa.append(x)
              ya.append(y)
              za.append(z)
              OK=1
              if np.sqrt(y*y+z*z)<REarth:
                if x > 0: OK=0
              if OK:
                xl.append(x)
                yl.append(y)
                zl.append(z)


#print("POS",x,y,z)

print(len(xa)/NAVG)

xl0=np.array(xl)
yl0=np.array(yl)
zl0=np.array(zl)



for ilat,lat in enumerate(lats):
    #print("Working on lat {}".format(lat))
    # again, need to reverse sign of lat tilt so positive latitude lifts up
    xobs1,yobs1,zobs1 = RotateY(xobs0,yobs0,zobs0,-lat)
    for ihr,hr in enumerate(hang):
        #print("Working on hour angle {}".format(hr))
        xobs2,yobs2,zobs2= RotateZ(xobs1,yobs1,zobs1,hr)
        xobs,yobs,zobs= RotateY(xobs2,yobs2,zobs2,tilt)

        satx =xl-xobs
        saty =yl-yobs
        satz =zl-zobs
            
        obslen = np.sqrt(xobs**2+yobs**2+zobs**2)
        satlen = np.sqrt(satx**2+saty**2+satz**2)
        dotprod = xobs*satx+yobs*saty+zobs*satz
        angle = np.abs(np.arccos(dotprod/(obslen*satlen)))
        flag = angle < horizonlimit
        for fg in flag: 
          if fg: count[ilat][ihr]+=1.

count/=NAVG

#print("Total count for {} {} is {}".format(ilat,ihr,count[ilat][ihr]))

#plt.figure()
#plt.scatter(yl,zl,s=1)

xl=np.array(xl)
yl=np.array(yl)
zl=np.array(zl)

#plt.figure()
#plt.scatter(xl/6378e3,zl/6378e3,s=1)
#plt.xlabel("x [RE]")
#plt.ylabel("y [RE]")
#plt.axes().set_aspect('equal')
#plt.savefig("satsXYprojection.png")

#fig=plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.set_aspect('equal')
#ax.scatter(xl,yl,zl,s=1)


print(len(hours),len(lats),len(count),len(count[0]))
#add in sunset times
lat1 = np.arange(-65.6,65.6,0.1)
lat = lat1*np.pi/180.
hrrise = -np.arccos( np.cos(90.833*np.pi/180.)/(np.cos(lat) * np.cos(tilt))-np.tan(lat) * np.tan(tilt) ) * 12./np.pi
if tilt<0:
  hrrise = np.append(-12,np.append(hrrise,0.))
  lat2 = np.append(lat1[0],np.append(lat1,lat1[-1]))
elif tilt>0: 
  hrrise = np.append(0.,np.append(hrrise,-12.))
  lat2 = np.append(lat1[0],np.append(lat1,lat1[-1]))
else:
  hrrise = np.append([-12,hrrise[0]],np.append(hrrise,[hrrise[-1],0.]))
  lat2 = np.append([-70.,-70.],np.append(lat1,[70.,70.]))
print(lat2)
print(hrrise)
print(tilt)
plt.figure()
lats*=180/np.pi
cb=plt.contourf(hours,lats,count,levels=255)
plt.contour(hours,lats,count,levels=20,colors='white',linewidths=0.5)
if tilt<0:
  print('tilt<0')
  plt.plot(hrrise, lat2, '-', c='#6E6E6E', lw=2.)
  plt.plot(-hrrise, lat2, '-', c='#6E6E6E', lw=2.)
  plt.fill_between(-hrrise,lat2,np.zeros(len(lat2))+70.,facecolor='#BDBDBD',edgecolor='None',alpha=0.6)
  plt.fill_between(hrrise,np.zeros(len(lat2))+70.,lat2,facecolor='#BDBDBD',edgecolor='None',alpha=0.6)
if tilt>0: 
  print('tilt>0')
  plt.plot(hrrise, lat2, '-', c='#6E6E6E', lw=2.)
  plt.plot(-hrrise, lat2, '-', c='#6E6E6E', lw=2.)
  plt.fill_between(-hrrise,np.zeros(len(lat2))-70.,lat2,facecolor='#BDBDBD',edgecolor='None',alpha=0.5)
  plt.fill_between(hrrise,lat2,np.zeros(len(lat2))-70.,facecolor='#BDBDBD',edgecolor='None',alpha=0.5)
if tilt==0: 
  print('tilt=0')
  plt.plot(hrrise, lat2, '-', c='#6E6E6E', lw=2.)
  plt.plot(-hrrise, lat2, '-', c='#6E6E6E', lw=2.)
  plt.fill_between(-hrrise,np.zeros(len(lat2))+70.,lat2,facecolor='#BDBDBD',edgecolor='None',alpha=0.5)
  plt.fill_between(hrrise,lat2,np.zeros(len(lat2))+70.,facecolor='#BDBDBD',edgecolor='None',alpha=0.5)

plt.title(TITLE1)
plt.xlabel("Hours Since Midnight")
plt.ylabel("Latitude (N+/S-)")
plt.hlines(50,-12,12,linestyle="-",color="black")
plt.hlines(29,-12,12,linestyle="--",color="blue")
plt.hlines(20,-12,12,linestyle="-.",color="green")
plt.hlines(-30,-12,12,linestyle=":",color="purple")
plt.colorbar(cb,label="Number of Illuminated Sats")
plt.savefig(FIGNAME+"heat.png")

plt.figure()
plt.title(TITLE1)
#print(lats[40])
legm30,=plt.plot(hours,count[40,:],label="Chilean Obs.",linestyle=":",color="purple")
#numattime=count[40,:]
#hrrisetime = -np.arccos( np.cos(90.833*np.pi/180.)/(np.cos(-30*np.pi/180.) * np.cos(tilt))-np.tan(-30*np.pi/180.) * np.tan(tilt) ) * 12./np.pi
#numattime1 = numattime[hours<hrrisetime]
#numattime2 = numattime[hours>hrrisetime]
#numattimeplot = (np.float(numattime1[-1])+np.float(numattime2[0]))/2.
#plt.plot([hrrisetime,-hrrisetime],[numattimeplot,numattimeplot],'o',mec="purple",mfc="#BDBDBD",ms=5.)
#print(lats[90])
leg20,=plt.plot(hours,count[90,:],label="Hawaiian Obs.",linestyle="-.",color="green")
#numattime=count[90,:]
#hrrisetime = -np.arccos( np.cos(90.833*np.pi/180.)/(np.cos(20.*np.pi/180.) * np.cos(tilt))-np.tan(20.*np.pi/180.) * np.tan(tilt) ) * 12./np.pi
#numattime1 = numattime[hours<hrrisetime]
#numattime2 = numattime[hours>hrrisetime]
#numattimeplot = (np.float(numattime1[-1])+np.float(numattime2[0]))/2.
#plt.plot([hrrisetime,-hrrisetime],[numattimeplot,numattimeplot],'o',mec="green",mfc="#BDBDBD",ms=5.)
#print(lats[99])
leg29,=plt.plot(hours,count[99,:],label="Canary Isl.",linestyle="--",color="blue")
#numattime=count[99,:]
#hrrisetime = -np.arccos( np.cos(90.833*np.pi/180.)/(np.cos(29.*np.pi/180.) * np.cos(tilt))-np.tan(29.*np.pi/180.) * np.tan(tilt) ) * 12./np.pi
#numattime1 = numattime[hours<hrrisetime]
#numattime2 = numattime[hours>hrrisetime]
#numattimeplot = (np.float(numattime1[-1])+np.float(numattime2[0]))/2.
#plt.plot([hrrisetime,-hrrisetime],[numattimeplot,numattimeplot],'o',mec="blue",mfc="#BDBDBD",ms=5.)
#print(lats[120])
leg50,=plt.plot(hours,count[120,:],label="Canadian Obs.",linestyle="-",color="black")
#numattime=count[120,:]
#hrrisetime = -np.arccos( np.cos(90.833*np.pi/180.)/(np.cos(50.*np.pi/180.) * np.cos(tilt))-np.tan(50.*np.pi/180.) * np.tan(tilt) ) * 12./np.pi
#numattime1 = numattime[hours<hrrisetime]
#numattime2 = numattime[hours>hrrisetime]
#numattimeplot = (np.float(numattime1[-1])+np.float(numattime2[0]))/2.
#plt.plot([hrrisetime,-hrrisetime],[numattimeplot,numattimeplot],'o',mec="black",mfc="#BDBDBD",ms=5.)
plt.xlabel("Hours Since Midnight")
plt.ylabel("Number of Illuminated Sats")
plt.xticks(np.linspace(-12,12,9))
plt.legend(handles=[leg50,leg29,leg20,legm30])
plt.savefig(FIGNAME+"curves.png")

#plt.figure()
#plt.plot(count[99,:]-count[90,:])

plt.show()





""
# ?keplertools

""


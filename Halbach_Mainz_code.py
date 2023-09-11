import matplotlib.pyplot as plt
import numpy as np

## INPUT SEGMENT
k=2;            #Halbach coefficient (1=dipole, 2= quadrupole,...)
ri=19.4/100;   #inner radius [m]
ro=20/100;     #outer radius [m](magnet diameter is 6mm)
N=54;           #number of segments(magnets)
typ=1;          #typ of Mandhala (1=cylindrical segment, 2 = circular, 3= trigonal,..., m= number of vertices...)
a=6/1000;       #side length of magnet [m], if = 0, a is calculated from d (next) 
d=15/1000;       #distance between magnets [m], if a~=0 d is ignored
L=10/1000;       #axial length [m] 

z0=(ro+ri)/2/np.sqrt(2*k+4);
Z0=[-z0,z0];    #position of multiple identical rings in z-dimension [m]
BR=0.950;         #remanence [T]
mr=1.05;        #relative permeability 

nx=100;         #number of points for field plots
ny=100;
scale=5;        #scaling for vector plot, pts = nx/scale 

##CALCULATION 
phi=np.pi/typ; #rotation of each Mandhala element
R=(ro+ri)/2;      #central radius

#calculation of ideal Halbach
if (k==1): 
    B0=BR*np.log(ro/ri)*1000;
    unit=' mT';
else:
    B0=BR*k/(k-1)*(ri**(1-k)-ro**(1-k));
    if (k==2):
        unit=' T/m';
    else:
       unit=' T/m^'+str(k-1);
    
#calculation of permittivity effect
theta=(mr+1)/(mr-1);
fm=(1-theta)*theta*ro**2/(ri**2-theta**2*ro**2);

#calculation of segmentation
fs=np.sin((k+1)*np.pi/N)/((k+1)*np.pi/N);

#calculation of polygonal effect
if (a!=0):
  if (typ==1):
    a=0;
  elif (typ==2):
     d=ro-ri-2*a;
  else:
    d=2*(ro-ri-a/np.sin(np.pi/typ));

if (typ ==1):
    A=np.pi*(ro**2-ri**2)-d*N*(ro-ri);
    fM=A/np.pi/(ro**2-ri**2);
elif (typ==2):
    a=(ro-ri)/2-d/2;
    A=a**2*np.pi;
    fM=N*A/np.pi/(ro**2-ri**2);
else:
    a=(ro-ri-d/2)*np.sin(np.pi/typ);
    A=typ*a**2/4/np.tan(np.pi/typ);
    fM=N*A/np.pi/(ro**2-ri**2);

##calc truncation in 3rd dimension
denomi=(4*R**2+L**2);
if (k==1):
    fL=L*(6*R**2+L**2)*denomi**(-3/2);
elif (k==2):
    fL=L*(L**4+10*L**2*R**2+30*R**4)*denomi**(-5/2);
elif (k==3):
    fL=L*(L**6+14*L**4*R**2+70*L**2*R**4+140*R**6)*denomi**(-7/2);
elif (k==4):
    fL=L*(L**8+18*L**6*R**2+126*L**4*R**4+420*L**2*R**6+630*R**8)*denomi**(-9/2);
elif (k==5):
    fL=L*(L**10+22*L**8*R**2+198*L**6*R**4+924*L**4*R**6+2310*L**2*R**8+2772*R**10)*denomi**(-11/2);
elif (k==6):
    fL=L*(L**12+26*L**10*R**2+286*L**8*R**4+1716*L**6*R**6+6006*L**4*R**8+12012*L**2*R**10+12012*R**12)*denomi**(-13/2)
else:
    fL=0;

## OUTPUT VALUES
alignvalue=90
print('__________________________________________________________________________________________')
print(' ')
pole=['DIPOLE','QUADRUPOLE','HEXUPOLE','OCTUPOLE','DECAPOLE','DODECAPOLE'];
if (k<=6):
   s='cylindrical inner Halbach '+pole[k-1]+':  with k = {0:.0f} ({1:0.0f} poles)'.format(k,2*k);
   print(s.center(alignvalue));
else:
   s='cylindrical inner Halbach {0:.0f}-UPOLE:  with k = {1:.0f} ({0:.0f} poles)'.format(2*k,k);
   print(s.center(alignvalue));
print(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -');
s='inner radius = {0:.3f} mm, outer radius = {1:.3f} mm, central radius = {2:.3f} mm'.format(ri*1000,ro*1000,R*1000);
print(s.rjust(alignvalue));
s='distance between magnets = {0:.3f} mm, magnet side length = {1:.3f} mm;  remanence = {2:.3f} T'.format(d*1000,a*1000,BR);
print(s.rjust(alignvalue));
print(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -');
flux='B';
for jj in range(k-1):
    flux=flux+"'";
s='ideal  >> '+flux+' = {0:.3f}'.format(B0);
s=s+unit;
print(s.rjust(alignvalue));
s='discretized in {0:.0f} segments  >> '.format(N)
s=s+flux;
s=s+' = {0:.3f}'.format(B0*fs);
s=s+unit;
print(s.rjust(alignvalue));

s='using a rel. permeability of {0:.0f} >> '.format(mr);
s=s+flux;
s=s+' = {0:.3f}'.format(B0*fs*fm);
s=s+unit;
print(s.rjust(alignvalue));

shape=['cylindrical segments','circles','triangles','squares','pentagons','hexagons','heptagons','octagons','nonagons','decagons'];
if (typ<=10):
    s='using {0:.0f} '.format(N);
    s=s+shape[typ-1];
    s=s+' with side length {0:.3f} mm >> '.format(a*1000);
    s=s+flux;
    s=s+' = {0:.3f}'.format(B0*fs*fm*fM);
    s=s+unit;
    print(s.rjust(alignvalue));
else:
    s='using {0:.0f} '.format(N);
    s=s+shape[typ-1];
    s=s+'-ogons with side length {0:.3f} mm >> '.format(a*1000);
    s=s+flux;
    s=s+' = {0:.3f}'.format(B0*fs*fm*fM);
    s=s+unit;
    print(s.rjust(alignvalue));

if (fL==0):
     s='truncated to a length of {0:.3f} mm >> '.format(L*1000);
     s=s+flux;
     s=s+'  **** for this k no factor could be calculated ';
     print(s.rjust(alignvalue));
else:
     s='truncated to a length of {0:.3f} mm >> '.format(L*1000);
     s=s+flux;
     s=s+' = {0:.3f}'.format(B0*fs*fm*fM*fL);
     s=s+unit;
     print(s.rjust(alignvalue));
                                                                  
##CALCULATE 2D B-FIELD
Xx=np.linspace(-ri,ri,nx); #coordinates
Yy=np.linspace(-ri,ri,ny);
Bx=np.zeros((nx,ny));  #Bx-component
By=np.zeros((nx,ny));  #By-component
Ba=np.zeros((nx,ny));  #magnitude 
for ix in range(nx):
    for iy in range(ny):
        if ((Xx[ix]**2+Yy[iy]**2)<=ri**2):  #only inside a circle with radius ri
            rk=complex(Xx[ix],Yy[iy]);
            if (k==1):
                Bx[ix,iy]=B0*fs*fm*fM*fL;
                Ba[ix,iy]=abs(np.sqrt(Bx[ix,iy]**2+By[ix,iy]**2));
            else:
                Rk=rk**(k-1);
                Bx[ix,iy]=B0*fs*fm*fM*fL*Rk.real;
                By[ix,iy]=-B0*fs*fm*fM*fL*Rk.imag;
                Ba[ix,iy]=abs(np.sqrt(Bx[ix,iy]**2+By[ix,iy]**2));


## PLOTS
fig, axs = plt.subplots(2,3);

# plot geometry in xy plane
pa=np.linspace(0,2*np.pi,1000);
axs[0, 0].plot(ro*np.cos(pa),ro*np.sin(pa),'c');  #plot outer diameter
axs[0, 0].plot(ri*np.cos(pa),ri*np.sin(pa),'c');  #plot inner diameter
axs[0, 0].plot(R*np.cos(pa),R*np.sin(pa),'g--');  #plot central diameter R
Px=np.zeros((typ,N));   #this is where the vertices are stored
Py=np.zeros((typ,N));
for jj in range(N):
    ux=R*np.cos((jj+1)*2*np.pi/N);   #center of magnets on R
    uy=R*np.sin((jj+1)*2*np.pi/N);
    La=R/10;                         #length of magnetization arrow
    axs[0, 0].plot(ux,uy,'ro');            #plot centers of magnets
    alph=(k+1)*(jj+1)*2*np.pi/N;     #orientation angle of magnetization
    axs[0, 0].plot([ux-La*np.cos(alph),ux+La*np.cos(alph)],[uy-La*np.sin(alph),uy+La*np.sin(alph)],'r'); #plot magnetization arrow
    axs[0, 0].plot(ux+La*np.cos(alph),uy+La*np.sin(alph),'r.');    #plot arrow tip
    if (typ==1):    #plot segments
         pia=2*np.pi*(jj+1)/N-np.pi/N+np.arctan(d/2/ri);
         pie=2*np.pi*(jj+1)/N+np.pi/N-np.arctan(d/2/ri);
         pci=np.linspace(pia,pie,50);
         axs[0, 0].plot(ri*np.cos(pci),ri*np.sin(pci),'k')
         poa=2*np.pi*(jj+1)/N-np.pi/N+np.arctan(d/2/ro);
         poe=2*np.pi*(jj+1)/N+np.pi/N-np.arctan(d/2/ro);
         pco=np.linspace(poa,poe,50);
         axs[0, 0].plot(ro*np.cos(pco),ro*np.sin(pco),'k')
         axs[0, 0].plot([ri*np.cos(pci[0]),ro*np.cos(pco[0])],[ri*np.sin(pci[0]),ro*np.sin(pco[0])],'k')
         axs[0, 0].plot([ri*np.cos(pci[-1]),ro*np.cos(pco[-1])],[ri*np.sin(pci[-1]),ro*np.sin(pco[-1])],'k')
    elif (typ==2): #plot circles
        pc=np.linspace(0,2*np.pi,100);
        axs[0, 0].plot(ux+a*np.cos(pc),uy+a*np.sin(pc),'k')
    else: #plot polygons
      ru=a/2/np.sin(np.pi/typ);   #side length 
      px=[ux-ru*np.cos(phi-alph)]; #start with last vertex
      py=[uy+ru*np.sin(phi-alph)];
      for kk in range(typ): #calc all vertices and store them in Px, Py
         Px[kk,jj]=ux-ru*np.cos(2*(kk+1)*np.pi/typ)*np.cos(phi-alph)+ru*np.sin(2*(kk+1)*np.pi/typ)*np.sin(phi-alph);
         Py[kk,jj]=uy+ru*np.cos(2*(kk+1)*np.pi/typ)*np.sin(phi-alph)+ru*np.sin(2*(kk+1)*np.pi/typ)*np.cos(phi-alph);
         px.append(Px[kk,jj]);
         py.append(Py[kk,jj]);
      axs[0, 0].plot(px,py,'k')
axs[0, 0].set_title('Arrangement of magnets in xy plane')
axs[0, 0].set_xlabel('x [m]')
axs[0, 0].set_ylabel('y [m]')
axs[0, 0].axis('equal')

#plot of geometry and field along z-direction
z=np.linspace(np.min(Z0)-L*2, np.max(Z0)+L*2,2000);
Bz=np.zeros(len(z));
for kk in range(len(Z0)):
    B=R**(2*k+3)*(R**2+(z-Z0[kk])**2)**(-k-1.5)*B0*fs*fm*fM*fL;
    axs[0, 1].plot(z,B,'c');
    Bz=Bz+B;
#output of central field in stacked rings
s='stacking {0:.0f} rings with Z0 distances:  in center >> '.format(len(Z0)); 
s=s+flux;
s=s+' = {0:.3f}'.format(Bz[int(len(Bz)/2)]);
s=s+unit;
print(s.rjust(alignvalue));
#continue plotting 
axs[0, 1].plot(z,Bz,'b',lw=2);
A=B0*fs*fm*fM*fL/2;
for kk in range(len(Z0)): #plot rectangles at position of magnet rings
  rec_x =[Z0[kk]-L/2,Z0[kk]-L/2,Z0[kk]+L/2,Z0[kk]+L/2,Z0[kk]-L/2];
  rec_y =[0,A,A,0,0];
  axs[0,1].plot(rec_x,rec_y,'k:');
axs[0,1].plot([np.min(Z0)-L*2,np.max(Z0)+L*2],[0,0],'k--'); #plot baseline 
s=flux+' along cylinder axis, z';
axs[0, 1].set_title(s);
s=flux+' ['+unit+']';
axs[0, 1].set_ylabel(s);
axs[0, 1].set_xlabel('z [m]')
axs[0, 1].axis('tight');


#plot magnitude of flux with overlayed vector plot
extent=-ri, ri, -ri, ri;
axs[0, 2].imshow(np.rot90(Ba,1), extent=extent);
axs[0, 2].set_title('magnitude with vectors');
axs[0, 2].set_xlabel('x [m]');
axs[0, 2].set_ylabel('y [m]');
#reduce resolution for vector plot by scale
dx=int(np.floor(nx/scale));
dy=int(np.floor(ny/scale));
Qx=np.zeros((dx,dy));
Qy=np.zeros((dx,dy));
Xq=np.zeros(dx);
Yq=np.zeros(dy);
jx=-1;
for ix in range(0,nx,scale): #now sort B values in a reduced matrix
    jx=jx+1;
    Xq[jx]=Xx[ix];
    jy=-1;
    for iy in range(0,ny,scale):
        jy=jy+1;
        Qx[jx,jy]=Bx[iy,ix];  #calc the transpose that's why x/y are swapped
        Qy[jx,jy]=By[iy,ix];
        Yq[jy]=Yy[iy];
axs[0, 2].quiver(Xq,Yq,Qx,Qy,color='w');

#plot magnitude of flux
extent=-ri, ri, -ri, ri
fig1=axs[1, 0].imshow(np.rot90(Ba,1), extent=extent)
axs[1, 0].set_title('magnitude of flux')
axs[1, 0].set_xlabel('x [m]')
axs[1, 0].set_ylabel('y [m]')
fig.colorbar(fig1, ax=axs[1,0]);

#plot Bx
fig2=axs[1, 1].imshow(np.rot90(Bx,1), extent=extent)
axs[1, 1].set_title('Bx component')
axs[1, 1].set_xlabel('x [m]')
axs[1, 1].set_ylabel('y [m]')
fig.colorbar(fig2, ax=axs[1,1]);

#plot By
fig3=axs[1, 2].imshow(np.rot90(By,1), extent=extent)
axs[1, 2].set_title('By component')
axs[1, 2].set_xlabel('x [m]')
axs[1, 2].set_ylabel('y [m]')
fig.colorbar(fig3, ax=axs[1,2]);

print(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -');
s='the vertices of the polygons are saved in the variables Px and Py';
print(s.center(alignvalue));
s='2nd plot shows {0:.0f} rings at positions in variable Z0'.format(len(Z0));
print(s.center(alignvalue));
print(' ');
print('__________________________________________________________________________________________');
plt.show()


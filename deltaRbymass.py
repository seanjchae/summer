import uproot as ROOT
import vector
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
from itertools import combinations

tree=ROOT.open('unweighted_events.root')['LHEF']

pt=tree['Particle.PT'].array()
eta=tree['Particle.Eta'].array()
phi=tree['Particle.Phi'].array()
mass=tree['Particle.M'].array()
pid=tree['Particle.PID'].array()
status=tree['Particle.Status'].array()

#mkw=(pid==24) & (status ==2)
mkj=(status ==1) & (abs(pid) <=4)

#wpt=pt[mkw]
jpt=pt[mkj]
jeta=eta[mkj]
jphi=phi[mkj]
jmass=mass[mkj]

vec_w1=[]
vec_w2=[]
j1=[]
j2=[]

for i in range(len(jpt)):
	
	#define quark contents
	q0=vector.obj(pt=jpt[i,0],phi=jphi[i,0],eta=jeta[i,0],mass=jmass[i,0])
	q1=vector.obj(pt=jpt[i,1],phi=jphi[i,1],eta=jeta[i,1],mass=jmass[i,1])
	q2=vector.obj(pt=jpt[i,2],phi=jphi[i,2],eta=jeta[i,2],mass=jmass[i,2])
	q3=vector.obj(pt=jpt[i,3],phi=jphi[i,3],eta=jeta[i,3],mass=jmass[i,3])
	#make list for 4jets
	cand=[q0,q1,q2,q3]
	
	#finding all possible combinations
	comb=list(combinations(cand,2))
	
	temp=0
	khi=[]
	
	for j in range(len(comb)):
		v1=comb[j][0]
		v2=comb[j][1]
		v3=comb[-1-j][0]
		v4=comb[-1-j][1]
		M1=(v1+v2).mass
		M2=(v3+v4).mass
		mw=80.379
		G=2.04759951
		khi.append(((M1-mw)**2/G + (M2-mw)**2/G)/2)
		temp=temp+1

		if j > 2:
			break

	for k in range(temp):
		if khi[k] == min(khi):
			vec_w1.append(comb[k][0]+comb[k][1])
			#vec_w2.append(comb[-1-k][0]+comb[-1-k][1])
			j1.append(v1)
			j2.append(v2)						
	if i%100000==99999:
		break


wmass=[]
for i in range(len(vec_w1)):
	wmass.append(vec_w1[i].mass)


eta1=[]
eta2=[]
phi1=[]
phi2=[]

for i in range(len(j1)):
	eta1.append(j1[i].eta)
	phi1.append(j1[i].phi)

for i in range(len(j2)):
	eta2.append(j2[i].eta)
	phi2.append(j2[i].phi)

dRresult=[]
for i in range(len(phi1)):
	if abs(phi1[i]-phi2[i]) >3.14 :
		deltaphi=6.28-abs(phi1[i]-phi2[i])
	else :
		deltaphi=abs(phi1[i]-phi2[i])

for i in range(len(eta1)):
	dR=((eta1[i]-eta2[i])**2+(deltaphi)**2)**(1/2)
	dRresult.append(dR)




'''
m1=[]
m2=[]
for i in range(len(vec_w1)):
	m1.append(vec_w1[i].mass)
	m2.append(vec_w2[i].mass)		
'''



import mplhep as hep
import matplotlib as mpl
plt.style.use(hep.style.CMS)
plt.rcParams['figure.figsize']=(12,10)
counts,xbins,ybins,image=plt.hist2d(wmass,dRresult,range=[[40,120],[0,10]],density=True,bins=40,norm=mpl.colors.LogNorm(),cmap=mpl.cm.jet)
plt.colorbar()

plt.xlim(40,120)
plt.ylim(0,10)
plt.minorticks_on()
plt.xlabel('mass(W)[GeV]')
plt.ylabel('dR')
plt.legend()
plt.show()

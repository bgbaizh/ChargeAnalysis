import sys
#sys1.path.append("../src")
import pyscal as pc
import pyscal.traj_process as ptp
import matplotlib.pyplot as plt
import numpy as np
import cChargeAnalysis 
import seaborn as sns


#file=files[40]
format="poscar"
sys1 = pc.System()
sys1.read_inputfile('POSCAR', format=format)#format{'lammps-dump', 'poscar', 'ase', 'mdtraj'}

bohr=0.529177208607388
threadnum=14
cChargcal=cChargeAnalysis.ChargeAnalysis()
cChargcal.threadnum=threadnum
cChargcal.debug=True
cChargcal.ChargeAnalysisInit('atlas.den',np.sum(sys1.box,axis=0)/bohr,sys1.natoms)
a=np.array(cChargcal.gridtest)
print(np.sum(a))
print(np.histogram(a.flatten()))
plt.hist(a.flatten())
plt.show()
'''
modes=[2]
#modes=[3]# 0 for voronoisum 1 for cutoffsum 2 for pdf 3 for boo
modevoronoisum=0
modecutoffsum=1
modepdf=2
modeboo=3
allmodenum=4

CDresults=[[] for i in range(allmodenum)]

for mode in modes:
    if mode==3:
        boocut=3.645/3/bohr
        histlow=0
        histbin=1000
        boolist=[6]
        atomscut = sys1.atoms[0:1]
        den=[0 for i in range(len(atomscut))]
        for atom in atomscut:
            den[atom.loc]=cChargcal.calculate_q_sumYdotCharge(boolist,np.array(atom.pos)/bohr,boocut,histlow,histbin)
        CDresults[mode]=den
                
    if mode==2:
        pdfcut=3.645/3/bohr
        histlow=0
        histbin=1000
        atomscut = sys1.atoms[0:1]
        den=[0 for i in range(len(atomscut))]
        for atom in atomscut:
            den[atom.loc]=cChargcal.CD_pdf(np.array(atom.pos)/bohr,pdfcut,histlow,histbin)
            #print(den[atom.loc])
        CDresults[mode]=den
    if mode==1:
        CDSumCut=3.645/2/bohr
        atomscut = sys1.atoms
        den=[0 for i in range(len(atomscut))]
        for atom in atomscut:
            den[atom.loc]=cChargcal.SumMethod_cut(np.array(atom.pos)/bohr,CDSumCut)
            #print(den[atom.loc])
        CDresults[mode]=den
    if mode==0:
        sys1.find_neighbors(method="voronoi")#method can be selected among "cutoff voronoi number", details can 
        atomsvoro = sys1.atoms
        den=[0 for i in range(len(atomsvoro))]

        for atom in atomsvoro:
            v3s=[]
            st = 1
            for vno in atom.face_vertices:
                vphase = atom.vertex_numbers[st:st+vno]
                ipos = atom.vertex_positions[vphase[0]]
                jpos = atom.vertex_positions[vphase[int((len(vphase)-1)/2)]]
                kpos = atom.vertex_positions[vphase[len(vphase)-1]]
                v3s.append([ipos,jpos,kpos])
                st += (vno+1)
            
            den[atom.loc]=cChargcal.SumMethod_voro(np.array(v3s)/bohr,np.array(atom.vertex_positions)/bohr)
            #print(den[atom.loc])
        CDresults[mode]=den
    '''
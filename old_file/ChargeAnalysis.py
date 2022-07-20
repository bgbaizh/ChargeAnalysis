from scipy import integrate
import numpy as np
#from sympy import *
import random
from ChargeAnalysisJacbian import jac
class ChargeAnalysis:
    def __init__(self) -> None:
        pass
    def ChargeAnalysisInit(self,rhofilename,boxs,natoms):
        
        self.boxs=boxs #self.boxs=[12.15,12.15,12.15]
        self.natoms=natoms
        self.readrho(rhofilename) 


    def set_vertex(self,v3s,vertex_positions):
        '''
        v3s is n*3*3* matrix where n is the number of faces
        Must run this before any Method_Function.
        vertex_positions is a parameter containing all vertexes, which is compatible with Glassviewer  
        '''
        
        self.v3s=v3s
        self.facesnum=len(self.v3s)
        self.X=[[i[0] for i in n]for n in self.v3s]
        self.Y=[[i[1] for i in n]for n in self.v3s]
        self.Z=[[i[2] for i in n]for n in self.v3s]
        #self.vol4=np.abs(np.cross((self.v3s[2]-self.v3s[0]),(self.v3s[1]-self.v3s[0])).dot(self.v3s[3]-self.v3s[0]))/6
        self.minxyz=[np.min([pos[0] for pos in vertex_positions]),np.min([pos[1] for pos in vertex_positions]),np.min([pos[2] for pos in vertex_positions])]
        self.maxxyz=[np.max([pos[0] for pos in vertex_positions]),np.max([pos[1] for pos in vertex_positions]),np.max([pos[2] for pos in vertex_positions])]
    def readrho(self,filename):
        '''
        this function should be used after self.boxs and natoms are correctly set
        '''
        rho0=[]
        with open(filename,'r') as f:
            nx,ny,nz=f.readline().split()
            nx=int(nx)
            ny=int(ny)
            nz=int(nz)
            for i in f:
                for j in i.split():
                    rho0.append(float(j))
        self.rho=np.zeros((nx,ny,nz))   
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    self.rho[i][j][k]=rho0[k*nx*ny+j*nx+i]
        self.nxyz= [nx,ny,nz]
        self.dx3s=[self.boxs[i]/self.nxyz[i] for i in range(3)]
        self.dxvol=self.dx3s[0]*self.dx3s[1]*self.dx3s[2]
        
    def IntegrationMethod(self,bins_scale=1):

        guessbin=(self.nxyz[0]*self.nxyz[1]*self.nxyz[2]/self.natoms)**(1/3)
        integratebins=int(guessbin*bins_scale)
        
        result=0
        for r in np.linspace(-1,1,integratebins,endpoint=False):
            for s in np.linspace(-1,1,integratebins,endpoint=False):
                for q in np.linspace(-1,1,integratebins,endpoint=False):
                    jacval=jac(r,s,q,self.a3s,self.b3s,self.c3s,self.d3s)
                    x=[0,0,0]
                    n=[0,0,0]
                    for i in range(3):
                        x[i]=self.a3s[i]*(1-r)*(1-s)*(1-q)/8+self.b3s[i]*(1+r)*(7-2*s-2*q+s*q)/24+self.c3s[i]*(1+s)*(7-2*r-2*q+r*q)/24+self.d3s[i]*(1+q)*(7-2*s-2*r+s*r)/24
                        n[i]=int(np.floor(x[i]/self.dx3s[i]))
                    result+=self.getden(n)*jacval*(2/integratebins)**3
        return result
        
    def SumMethod(self):

        minxyzn=[int(np.floor(self.minxyz[i]/self.dx3s[i])) for i in range(3)]
        maxxyzn=[int(np.floor(self.maxxyz[i]/self.dx3s[i])) for i in range(3)]
        centerofvetex=np.sum(np.sum(self.v3s,axis=0),axis=0)/self.facesnum/3
        nvector=[]
        surfaceD=[]
        #test=np.zeros(self.nxyz)
        for n in range(self.facesnum):
            i=0
            j=1
            k=2
            nvector.append([(self.Y[n][j]*self.Z[n][k] - self.Y[n][j]*self.Z[n][i] - self.Y[n][i]*self.Z[n][k] - self.Y[n][k]*self.Z[n][j] + self.Y[n][i]*self.Z[n][j] + self.Y[n][k]*self.Z[n][i]),
                            (self.X[n][k]*self.Z[n][j] - self.X[n][i]*self.Z[n][j] - self.X[n][k]*self.Z[n][i] - self.X[n][j]*self.Z[n][k] + self.X[n][j]*self.Z[n][i] + self.X[n][i]*self.Z[n][k]),
                            (self.X[n][j]*self.Y[n][k] - self.X[n][j]*self.Y[n][i] - self.X[n][i]*self.Y[n][k] - self.X[n][k]*self.Y[n][j] + self.X[n][k]*self.Y[n][i] + self.X[n][i]*self.Y[n][j])
                            ])
            surfaceD.append(-np.dot(nvector[n],[self.X[n][i],self.Y[n][i],self.Z[n][i]]))
        nvector=np.array(nvector)
        surfaceD=np.array(surfaceD)
        nvectornorm=np.array([np.linalg.norm(n) for n in nvector])

        centerDistance=[]
        for i in range(self.facesnum):

            distemp=(np.dot(nvector[i],centerofvetex)+surfaceD[i])/nvectornorm[i]
            #法向量指向外，这个距离为负，法向量指向内部，这个距离为正
            if distemp>0:
                nvector[i]=-nvector[i]#让所有法向量指向外
                surfaceD[i]=-surfaceD[i]
            if distemp<0:
                distemp=-distemp
            centerDistance.append(distemp)
        centerDistance=np.array(centerDistance)
        #test1=0
        #test2=0
        result=0
        for i in range(minxyzn[0],maxxyzn[0]+1):
            for j in range(minxyzn[1],maxxyzn[1]+1):
                kstart='nan'
                kend='nan'
                for k in range(minxyzn[2],maxxyzn[2]+1):
                    
                    gridmidpoint=[i*self.dx3s[0]+self.dx3s[0]/2,j*self.dx3s[1]+self.dx3s[1]/2,k*self.dx3s[2]+self.dx3s[2]/2]
                    centertopoint_project_nvector=[]
                    for n in range(self.facesnum):
                        #if j%96==1:
                        #    print('here',i,j,k,n)
                        centertopoint_project_nvector.append(np.dot((gridmidpoint-centerofvetex),nvector[n])/nvectornorm[n])
                    boollist=np.array([centertopoint_project_nvector[n]<=centerDistance[n] for n in range(self.facesnum)])
                    #test2+=1
                    if boollist.all():
                        result+=self.getden([i,j,k])*self.dxvol
                        #test1+=1
                        #test[i%self.nxyz[0]][j%self.nxyz[1]][k%self.nxyz[2]]+=1
                        if k+1<=maxxyzn[2]:
                            kstart=int(k+1)
                            break
                if kstart!='nan':
                    for k in range(kstart,maxxyzn[2]+1)[::-1]:
                        gridmidpoint=[i*self.dx3s[0]+self.dx3s[0]/2,j*self.dx3s[1]+self.dx3s[1]/2,k*self.dx3s[2]+self.dx3s[2]/2]
                        centertopoint_project_nvector=[]
                        for n in range(self.facesnum):
                            centertopoint_project_nvector.append(np.dot((gridmidpoint-centerofvetex),nvector[n])/nvectornorm[n])
                        boollist=np.array([centertopoint_project_nvector[n]<=centerDistance[n] for n in range(self.facesnum)])
                        #test2+=1
                        if boollist.all():
                            result+=self.getden([i,j,k])*self.dxvol
                            #test1+=1
                            #test[i%self.nxyz[0]][j%self.nxyz[1]][k%self.nxyz[2]]+=1
                            if k-1>=kstart:
                                kend=int(k-1)
                                break
                    if kend!='nan':
                        for k in range(kstart,kend+1):
                            result+=self.getden([i,j,k])*self.dxvol
                            #test1+=1
                            #test2+=1
                            #test[i%self.nxyz[0]][j%self.nxyz[1]][k%self.nxyz[2]]+=1
                            
                            
        #print(test1,test2)
        return result#,test
    
    def MCMethod(self,nMC=100):
        minxyz=[np.min([i[j] for i in self.v3s]) for j in range(3)]
        maxxyz=[np.max([i[j] for i in self.v3s]) for j in range(3)]
        minxyzn=[int(np.floor(minxyz[i]/self.dx3s[i])) for i in range(3)]
        maxxyzn=[int(np.floor(maxxyz[i]/self.dx3s[i])) for i in range(3)]
        centerofvetex=np.sum(self.v3s,axis=0)/4
        nvector=[]
        surfaceD=[]
        for i in range(4):
            j=(i+1)%4
            k=(i+2)%4
            nvector.append([(self.Y[j]*self.Z[k] - self.Y[j]*self.Z[i] - self.Y[i]*self.Z[k] - self.Y[k]*self.Z[j] + self.Y[i]*self.Z[j] + self.Y[k]*self.Z[i]),
                            (self.X[k]*self.Z[j] - self.X[i]*self.Z[j] - self.X[k]*self.Z[i] - self.X[j]*self.Z[k] + self.X[j]*self.Z[i] + self.X[i]*self.Z[k]),
                            (self.X[j]*self.Y[k] - self.X[j]*self.Y[i] - self.X[i]*self.Y[k] - self.X[k]*self.Y[j] + self.X[k]*self.Y[i] + self.X[i]*self.Y[j])
                            ])
            surfaceD.append(-np.dot(nvector[i],[self.X[i],self.Y[i],self.Z[i]]))
        nvector=np.array(nvector)
        surfaceD=np.array(surfaceD)
        nvectornorm=np.array([np.linalg.norm(i) for i in nvector])

        centerDistance=[]
        for i in range(4):

            distemp=(np.dot(nvector[i],centerofvetex)+surfaceD[i])/nvectornorm[i]
            #法向量指向外，这个距离为负，法向量指向内部，这个距离为正
            if distemp>0:
                nvector[i]=-nvector[i]#让所有法向量指向外
                surfaceD[i]=-surfaceD[i]
            if distemp<0:
                distemp=-distemp
            centerDistance.append(distemp)
        centerDistance=np.array(centerDistance)

        result=0
        for i in range(minxyzn[0],maxxyzn[0]+1):
            for j in range(minxyzn[1],maxxyzn[1]+1):
                for k in range(minxyzn[2],maxxyzn[2]+1):
                    p=np.array([[i*self.dx3s[0]+self.dx3s[0]*l1,j*self.dx3s[1]+self.dx3s[1]*l2,k*self.dx3s[2]+self.dx3s[2]*l3] for l1,l2,l3 in [
                        [0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],[1,0,1],[1,1,0],[1,1,1]]])
                    
                    p_in_num=0
                    for i2 in range(8):
                        centertopoint_project_nvector=[]
                        for m in range(4):
                            centertopoint_project_nvector.append(np.dot((p[i2]-centerofvetex),nvector[m])/nvectornorm[m])
                        boollist=np.array([centertopoint_project_nvector[z]<=centerDistance[z] for z in range(4)])
                        if boollist.all():
                            p_in_num+=1
                    if p_in_num==8:
                        result+=self.getden([i,j,k])*self.dxvol
                    elif p_in_num==0:
                        continue
                    else:
                        
                        count=0
                        for i3 in range(0,nMC):
                            boollist=[]
                            centertopoint_project_nvector_random=[]
                            x=random.uniform(i*self.dx3s[0],i*self.dx3s[0]+self.dx3s[0])
                            y=random.uniform(j*self.dx3s[1],j*self.dx3s[1]+self.dx3s[1])
                            z=random.uniform(k*self.dx3s[2],k*self.dx3s[2]+self.dx3s[2])
                            for m in range(4):
                                centertopoint_project_nvector_random.append(np.dot((np.array([x,y,z])-centerofvetex),nvector[m])/nvectornorm[m])
                            boollist=np.array([centertopoint_project_nvector_random[z]<=centerDistance[z] for z in range(4)])
                            if boollist.all():
                                count+=1
                        result+=self.getden([i,j,k])*self.dxvol*count/nMC
                        
        return result
    
    def getden(self,i3s):
        for i in range(3):
            if i3s[i]>self.nxyz[i]-1 or i3s[i]+1<0:
                i3s[i]=i3s[i]%self.nxyz[i] 
        return self.rho[i3s[0]][i3s[1]][i3s[2]]

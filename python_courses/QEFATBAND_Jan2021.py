#!/usr/bin/env python37
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 21 12:13:53 2018

@author: yavar
"""


#######INSTRUCTION
####>>  a=QElog("test.out")
####>>  a.plotband()

import atexit,os
from os import path
import sys, string
import numpy as np
import scipy.stats as stats
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import re





new_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']

_PW_ELEC = '     number of electrons       ='
_PW_NBND = 'number of Kohn-Sham states='
_PW_NK = '     number of k points= '
_PW_KVECT = '  k ='


_PDOS_STATE = '     state #'
_PDOS_K = ' k = '
_PDOS_PSI = '     psi ='
_PDOS_abs = '    |psi|^2 ='
_PDOS_EIG = '===='
_PDOS_LOWDIN = '===='

class QElog(object):
    
    def __init__(self,log):
        self.log=log
        
##SPIN FORMAT CHANGED IN NEW VERSION
        
    
 #   @staticmethod
 #   def validate_log(self):
        if os.path.exists(log):
            print("your file is in path")
            if os.stat(self.log).st_size<5:
                raise ValueError("File Size NOT Normal")
            else:
                print("it seems normal file")
        else:
            print("!!!!!!!!!!file not found, Enter correct path..........")
       
    
#    def e_numb(self):
#        with open(self.log) as f:
#            for line in f:
#                if _PW_ELEC in line:
#                    enumb=line.split()[-1]
#            enum=float(enumb)
#        return enum
    @staticmethod
    def read_output(logname):
         indexes = {
                 _PW_ELEC: [],
                 _PW_NBND: [],
                 _PW_NK: [],
                 _PW_KVECT: []
                 }
         # OPEN File
         fileobj=open(logname, 'r')
         pwo_lines = fileobj.readlines()
         
         
         
         # If the file in spin polarized or not
         
         spinval=1
         for lin in pwo_lines:
             if re.search('SPIN UP', lin):
                 spinval=2
         
         #Extract index of lines contating important data
         for idx, line in enumerate(pwo_lines):
             for identifier in indexes:
                 if identifier in line:
                     indexes[identifier].append(idx)
 #                else:
 #                    raise ValueError("! check your file again, something went wrong...")
         
         
         #mahe dictionary to write values in it
         for elec_index in indexes[_PW_ELEC]:
             electnum = float(pwo_lines[elec_index].split()[-1])
             
         data={'nelec':electnum}
         
         
         for nbnd_index in indexes[_PW_NBND]:
             nbnd = int(pwo_lines[nbnd_index].split()[-1])
             
             
             
         data.update({'nband':nbnd})
         data.update({'SPIN':spinval})
         
#         print(data)
         for nk_index in indexes[_PW_NK]:
             knum = int(pwo_lines[nk_index].split()[4])
             
         data.update({'nk':knum})
 #        print(data)
         
        # Extract k vectors and update dictionary         
         kvector=[]
         for kvect_index in indexes[_PW_KVECT]:
             kvect = pwo_lines[kvect_index].replace('  k =','')
             kvect = kvect.replace('-', ' -')
             kvect = kvect.split()[0:3]
             kvect = list(map(float, kvect))
             kvector.append(kvect)
         karray=np.array(kvector)
         
         data.update({'kvectors':karray})
         
         
         
         j=0
         banddata=np.zeros((nbnd,knum*spinval))
         lnum=int(nbnd/8)
         for kvect_index in indexes[_PW_KVECT]:
             j=j+1
             b=np.empty((0))
             for i in range(lnum+2):
                 a=pwo_lines[kvect_index+i+1].strip().split()
                 a=np.array(a,dtype=float)
                 b=np.append(b,a)
             banddata[:,j-1]=b
                 
         data.update({'bandstruct':banddata})
         fileobj.close()
         return data                
   
         
    @property     
    def nelect(self):
        """contains nelec, nbnd, nkpoint, kvectors
        """
        return self.read_output(self.log)['nelec']
    
    @property
    def nband(self):
        return self.read_output(self.log)['nband']
    
    @property
    def nkpoints(self):
        return self.read_output(self.log)['nk']
    
    
    @property
    def spin(self):
        return self.read_output(self.log)['SPIN']
    
    @property
    def kvectors(self):
        return self.read_output(self.log)['kvectors']
    
    @property
    def eigenvalues(self):
        return self.read_output(self.log)['bandstruct']


    @property
    def vband(self):
        nelec=self.nelect
        nh=int((nelec/2))
        homos=self.eigenvalues[nh-1,:]
#       print('spin polarization is', self.spin)
        return homos
    @property     
    def cband(self):
        nelec=self.nelect
        nh=int((nelec/2))
        lumos=self.eigenvalues[nh,:]
#       print('spin polarization is', self.spin)
        return lumos#

    @property    
    def HOMO(self):
        homo=max(self.vband)
        return homo
    
    @property  
    def LUMO(self):
        lumo=min(self.cband)
        return lumo
    
    
  
    @property   
    def GAP(self):
        gap=self.LUMO-self.HOMO
        return gap
    
    @staticmethod
    
    
    
    
    
    def read_splabel(inpfile):
        labels=[]
        intervals=[0]
        with open(inpfile, "r") as ifile:
            for line in ifile:
                if "crystal_b" in line:
                    num_sp=int(next(ifile).strip())
                    for x in range(num_sp):
                        text=next(ifile).split()
                        a=(int(text[3])+intervals[x])
                        labels.append(text[5])
                        intervals.append(a)
                    intervals=intervals[:num_sp]
        intemp=intervals[1:]
        
        #  check if there is a discontin  in k vector
        discont=[]
        for i in range(len(intemp)):
            discont.append(intemp[i]-intervals[i])
        indices = [i for i, x in enumerate(discont) if x == 1]
        # return indices for discontiniuties
        discontind=[intervals[i] for i in indices]
        return intervals, discontind , labels
                        
    
    
    
    def plotband(self, fermi, name="test.png", low=2, up=2):
        
        
        outname=self.log
        inpname=re.sub(".out", ".in", outname)
        
        
        
        if path.exists(inpname):
            inpfile=inpname
        else:
            inpfile=input("""  The proper input file can not be found  
                          Enter name of band input Or file with QE crystal_b format  """)

        
        
        fig, ax = plt.subplots()
        
        
        nk=self.nkpoints
        x=np.zeros((nk,1))
        kvector=self.kvectors
        banddata=self.eigenvalues
        spin=self.spin
        
        
        splist=self.read_splabel(inpfile)
        ind0=splist[0]
   

               
        disind=splist[1]
        xindex=splist[2]
        
        kdim=np.shape(kvector)[0]       
        
        checker=np.ones(kdim)
        
        for iind in disind:
            checker[iind]=0
        
        
        
        for i in range(nk-1):
            x[i+1]=checker[i]*np.abs(np.linalg.norm(kvector[i+1]-kvector[i]))+x[i]
            
        ind=x[ind0]
        ind1d=ind.flatten()

        xindex=splist[2]
        print(ind1d)
        print(xindex)
        

        
               
        for k in range(0, int(self.nband)):
            yup=banddata[k,0:nk]
            ydn=banddata[k,nk*(spin-1):(spin)*nk]
            ax.plot(x,yup-fermi,  linewidth=1.5 , c= new_colors[9],zorder=int(1+(spin-1)*2))
            ax.plot(x,ydn-fermi,  linewidth=3 , color="#1a214c", zorder=2)            
        ax.set_xticks(ind1d)
        ax.set_xticklabels(xindex, fontsize=18)
        ax.fill_between(x.flatten(), self.LUMO-fermi, y2=self.HOMO-fermi, alpha=0.3, color="#616a6b" )
        for xc in ind:
            ax.axvline(x=xc, color="gray", linewidth=2.1, zorder=3)
            
        y1=-low
        y2=self.LUMO-fermi+up    
        
        y1int=int(y1)
        yvals=[]
        for i in np.arange(y1int, y2, 1.0):
            yvals.append(i)
        
        ax.set_yticks(yvals)
        ax.set_yticklabels(yvals, fontsize=20)
        ax.yaxis.set_minor_locator(MultipleLocator(5))
        ax.axhline(y=0, linewidth=1.5, color="#4C3F64")
        ax.axhline(y=0, linewidth=0.75, color="#D0BDF2")
        y1=-low
        y2=self.LUMO-fermi+up
        ax.set_ylim(y1,y2)          
        ax.set_xlim(x[0],x[nk-1]) 
        ax.set_ylabel("Energy (eV)", fontsize=20)
        
        fig.set_size_inches(6,6)
        fig.savefig(name,dpi=300, bbox_inches='tight')
    

   
    @staticmethod
    def read_projout(projout):
        """
        Open the output file of projwfc and returns a DIC (projdata)
        which contains states, kvectors, and wfc coefficients for each psi_k
        
        """
               
        indexes = {
                _PDOS_STATE : [],
                _PDOS_K : [],
                _PDOS_PSI : [],
                _PDOS_EIG :[],
                _PDOS_abs : [],
                _PDOS_LOWDIN : []
                 }        
        
        
               
        fileobj=open(projout, 'r')
        proj_lines = fileobj.readlines()
             
             
        #Extract index of lines contating important data
        for idx, line in enumerate(proj_lines):
            for identifier in indexes:
                if identifier in line:
                    indexes[identifier].append(idx)
        i=0
        state=[]
        char_list = ['state #', 'atom', '\(', '\)', '\:', ', wfc']
        for stateidx in indexes[_PDOS_STATE]:
            modfline=re.sub("|".join(char_list), " ", proj_lines[stateidx])
            mlinelist=modfline.split()
            state.append(mlinelist)
            i=i+1
        pdosdata=state
                 
        filter1=[x[0:3] for x in state]
        
        
        filterdata=np.array(filter1)
        
        natoms=int(pdosdata[-1][1])
        symbol=[x[2] for x in pdosdata]
        symbols=list(set(symbol))
        projdata={'natoms':natoms} 
    
        
        atomnum=[x[1] for x in pdosdata]
        atomnumbsl=[int(x) for x in atomnum]
        
        
        stat_separator=[]
        for i in range(natoms):
            sep=atomnumbsl.index(i+1)
            stat_separator.append(sep)
     
    
        nstates=len(atomnum)
      
        stat_separator.append(nstates)
        
          
        projdata.update({'nstates':nstates})
        projdata.update({'statesofatoms':filterdata})
# Find start and end index for each atom 
#  for example for atom #1 stat_separator[0] is start and stat_separator[1] is end    
        projdata.update({'statesind':stat_separator})
    
    
        projdata.update({'symbols':symbols})


        eig=[]
        for eig_index in indexes[_PDOS_EIG]:
            eigen=proj_lines[eig_index].replace('====', '').split()[3]
            eigen=float(eigen)
            eig.append(eigen)


        eig=np.array(eig)
        projdata.update({'eigvalues':eig})
               
        
        
        kvector=[]
        for k_index in indexes[_PDOS_K]:
            kvect=proj_lines[k_index].replace(' k =  ', '').split()[0:3]
            kvect=list(map(float,kvect))
            kvector.append(kvect)


        kvector=np.array(kvector)
        projdata.update({'kvectors':kvector})
        
        
#extract lines contatinig coeffiecients        
        incoeffline=zip(indexes[_PDOS_PSI],indexes[_PDOS_abs])
        coefficients=[]
        for indxll in  incoeffline:
            a=indxll[0]
            b=indxll[1]
            test=proj_lines[a:b]
            g=[]
            for lines in test:
                coeffs=re.findall(r"[-+]?\d*\.\d+|\d+", lines)
                g=g+coeffs
            coefficients.append(g)
        
        atomgroup=[]
        b=filter1
        for atom in symbols:
            res=[atom]
            indexes=[s[0] for s in b if s[2]==atom]
            res.extend(indexes)
            atomgroup.append(res)
            
            
        with open('groups.txt', 'w') as f:
            for item in atomgroup:
                for mem in item:
                    f.write("%s " % mem)
                f.write("\n")
        
        projdata.update({'proj_wfc':coefficients})
        projdata.update({'indexofatoms':atomgroup})
        return projdata
    
    def findcoef(self, projdosout):
        
        
        #extract contribution of each atom in this 
        data00=self.read_projout(projdosout)
        #read coefficients
        states=data00['statesofatoms']
        
        statnumbs=len(states)
        
        atomlist=[states[0][2]]
        atomindlist=[states[0][1]]
        
        for i in range(statnumbs-1):
            if ( states[i+1][1] != atomindlist[-1] ):
                atomlist.append(states[i+1][2])
                atomindlist.append(states[i+1][1])

        

        
        homoind=int(self.nelect-1+(np.where(self.vband==self.HOMO)[0][0])*self.nband)
        lumoind=int(self.nelect+(np.where(self.cband==self.LUMO)[0][0])*self.nband)
        
     
        
        homoval=data00['proj_wfc'][homoind][0::2]
        lumoval=data00['proj_wfc'][lumoind][0::2]
        
        
        homocoef=data00['proj_wfc'][homoind][1::2]
        lumocoef=data00['proj_wfc'][lumoind][1::2]
        
        homocoef=np.array(homocoef, dtype=int)
        lumocoef=np.array(lumocoef, dtype=int)
           
               
        homoval=np.array(homoval, dtype=float)
        lumoval=np.array(lumoval, dtype=float)
        
        natoms=int(states[-1][1])

        homodata=np.vstack((homocoef, homoval)).T
        lumodata=np.vstack((lumocoef, lumoval)).T
        
        
        strings=[]
        for i in range(natoms):
            strings.append(atomlist[i]+"#"+atomindlist[i])
        
        strings01=np.array(strings)
        strings01=np.reshape(strings01,(natoms,1))
        
        homoatom=np.zeros((natoms,1))
           
           
        for i in range(np.shape(homodata)[0]):
            ind=int(homodata[i,0])
            value=homodata[i,1]
            atomind=int(states[ind-1][1])
            homoatom[atomind-1]+=value
        
        homodata00=np.hstack((homoatom, strings01))
        
        homodata01=homodata00[homodata00[:,0].argsort()[::-1]]
        
        
        
        lumoatom=np.zeros((natoms,1))
        
        
           
        for i in range(np.shape(lumodata)[0]):
            ind=int(lumodata[i,0])
            value=lumodata[i,1]
            atomind=int(states[ind-1][1])
            lumoatom[atomind-1]=value+lumoatom[atomind-1]
            
            
        lumodata00=np.hstack((lumoatom, strings01))  
        
        lumodata01=lumodata00[lumodata00[:,0].argsort()[::-1]]

        
        homomask=[homoval > homoval[0]/10.0]
        lumomask=[lumoval > lumoval[0]/10.0]

        
        homodelocorb=stats.variation(homoval)
        lumodelocorb=stats.variation(lumoval)
               
 
        homodelocatoms=stats.variation(homoatom)
        lumodelocatoms=stats.variation(lumoatom)

        
        
        
        print("delocalization of the homo over orbitals is   :  ", homodelocorb)
        print("delocalization of the homo over atoms is   :  ", homodelocatoms)

        
        print("delocalization of the lumo is   :  ", lumodelocorb)

        print("number of atoms    :"  , natoms)
        print("sum all contrib:   ", np.sum(homoatom))
        
        return homocoef, lumocoef , homodata01,  lumodata01
        
       
        
    def checkproj(self, projdosout):
        """
        this function extract projection of band J on different states
        
        """
        nkout=self.nkpoints
        kvectorout=self.kvectors
        nbandout=self.nband
        coffnumb=nbandout*nkout
        
        data0=self.read_projout(projdosout)
        
        coefnumproj=len(data0['proj_wfc'])
        kvectorproj=data0['kvectors']
                    
        check =  coffnumb==coefnumproj and np.array_equal(kvectorproj,kvectorout)
        
        return check




    def conv_proj_table(self, projdosout):
        
        
        """
        This function converts all data extracted by projdata to a tacle with following format \n
        
        
        |            | orb1         |  orb2         |     ..n is total orbitals of all atoms...|orb_n | \n 
        | band1,k1       | coeff(float) |  coeff(float) |     ......                               |      | \n
        | band2,k1
        | band3,k1             |               |                                                 | \n
        | ....       |              |               |                                          |      | \n                                                                                  
        
        
        
        but we need to resort above table in
        
        |band1,k1
        |band1,k2
        |band1,k3
        
        then we need a convertor 
        """
        data000=self.read_projout(projdosout)
        
        #read coefficients
        states=data000['statesofatoms']
                
        #number of quantum states
        nstate=states.shape[0]
        
        n_allstates=self.nkpoints*self.nband
        nk=self.nkpoints
        nbands=self.nband
        
        projband_mat=np.zeros((n_allstates, nstate))
#        for i in range(n_allstates):
            
        stateind=np.array(range(1,nstate+1))
        for i in range(n_allstates):
            temp=data000['proj_wfc'][i]
            listofstates=list(map(int, temp[1::2])) 
            listofstates=np.array(listofstates)
            values=list(map(float,temp[0::2]))    
            values=np.array(values)
            exsistate=np.intersect1d(listofstates,stateind)
            
            # this is list
            list1=list(map(int, exsistate))
            
            for k in list1:
                ind=np.where(listofstates==k)
                projband_mat[i,k-1]=values[ind]
#       

        resortprojcoeff=np.zeros_like(projband_mat)
        for iband in range(nbands):
            for ik in range(nk):
                res=projband_mat[ik*nbands+iband,:]
                resortprojcoeff[iband*nk+ik,:]=res
        np.savetxt("projban_mat0_0.txt", projband_mat, fmt="%8.4f")
        np.savetxt("projban_mat0_1.txt", resortprojcoeff, fmt="%8.4f")      
        return resortprojcoeff       
#projband is zeros for all states        
#        projband=np.zeros((n_allstates)
        
        
        
    def bandproject_group(self, projdosout, groupfile="groups.txt"):
        
        """
        This function reads gruops file and then calculate scaled projband for each group
        
        """
        
        

        with open(groupfile, "r") as f:
            statname=[]
            statinds=[]
            for line in f:
                statname.append(line.split()[0])
                statinds.append(line.split()[1:])
        f.close()
        
        
        numb=len(statname)
        
        
        #cring converted proj and orbital coefficients
        
        projband00=self.conv_proj_table(projdosout)
        n=np.shape(projband00)[0]
        
        banddata=self.eigenvalues
        kvector=self.kvectors
        nk=self.nkpoints
        x=np.zeros((nk,))
        for i in range(nk-1):
            x[i+1]=np.abs(np.linalg.norm(kvector[i+1]-kvector[i]))+x[i]
            
        eignumb=np.shape(banddata)[0]
        numball=np.shape(banddata)[0]*np.shape(banddata)[1]
        mybanddata=np.reshape(banddata, (numball))
        
        repeatedx=np.zeros_like(mybanddata)
        
        for i in range(eignumb):
            repeatedx[i*nk:(i+1)*nk]=x
        
        projected_group=np.zeros((n,numb+2))

        projected_group[:,0]=repeatedx
        projected_group[:,1]=mybanddata
        for i in range(numb):
            groupinds=list(map(int, statinds[i]))
            for j in groupinds:
                projected_group[:,i+2] += projband00[:,j-1]   
                
                
        
        
                    
     
        np.savetxt("proj.data", projected_group,fmt="%8.4f")
        
        

        with open('projgnu.txt', 'w') as f:
            for i in range(eignumb):
                for k in range (nk):
                    row=i*nk+k
                    tmp=np.array_str(projected_group[row,:], max_line_width=None, precision=6, suppress_small=None)
                    tmp=tmp.replace("[", " ")
                    tmp=tmp.replace("]", " ")
                    tmp0=tmp.split()
                    for item in tmp0:
                        a=float(item)
                        f.write('%8.4f ' % a)
                    f.write("\n")
                f.write("\n")
            f.close()
        
        
        
        return projected_group
    
    
             
    
    
    
    
    
    
    
    
    def plotfatband(self, fermi, projdosout, name="test.png", groupfile="groups.txt",  low=1.5, up=1.5):
        nk=self.nkpoints
        x=np.zeros((nk,1))
        kvector=self.kvectors
        #banddata=self.eigenvalues
        for i in range(nk-1):
            x[i+1]=np.abs(np.linalg.norm(kvector[i+1]-kvector[i]))+x[i]
        splist=self.read_splabel()
        ind0=splist[0]
        ind=x[ind0]
        xindex=splist[1]
        
#insert  fat coefficients        
#norm = plt.Normalize(dydx.min(), dydx.max())

        
        projbanddata=self.bandproject_group(projdosout, groupfile)
    #   first column in xk and second is eigen and 3th 4th ...are weights    
        lband=projbanddata[:,1]

        
        with open(groupfile, "r") as f:
            groupnames=[]
            for line in f:
                groupnames.append(line.split()[0])
        f.close()
        
       

        for iband in range(0, int(self.nband/self.spin)):
            spin=self.spin
            upn=iband*nk*spin
            yup=lband[upn:upn+nk]
            plt.plot(x,yup-fermi,  linewidth=1. ,  c= new_colors[9], alpha=0.4, zorder=1)
            dnn=(spin-1)*nk+upn
            ydn=lband[dnn:dnn+nk]
            plt.plot(x,ydn-fermi,  linewidth=1 ,  color=new_colors[4], alpha=0.4,zorder=2)            
        
        ngroups=len(groupnames)
        for q in range(ngroups):
            gname=groupnames[q]
            plt.scatter(x, -1000*x-100, c=new_colors[q+2], alpha=0.8-q*0.1,label=gname)
            size=projbanddata[:,q+2]
            for iband in range(0, int(self.nband/self.spin)):
                spin=self.spin
                upn=iband*nk*spin
                yup=lband[upn:upn+nk]
                usize=size[upn:upn+nk]
                plt.scatter(x,yup-fermi,  s=15*usize , c= new_colors[q+2],alpha=0.8-q*0.1, zorder=3)
                dnn=(spin-1)*nk+upn
                ydn=lband[dnn:dnn+nk]
                dsize=size[dnn:dnn+nk]
                plt.scatter(x,ydn-fermi,  s=15*dsize, color=new_colors[q+2], alpha=0.8-q*0.1,zorder=4)            
            
            
      
            

        plt.xticks(ind, xindex, fontsize =22)
        for xc in ind:
            plt.axvline(x=xc, color="gray")
        plt.yticks(fontsize =22)
        plt.axhline(y=0, linewidth=1.5, color="#4C3F64")
        plt.axhline(y=0, linewidth=0.75, color="#D0BDF2")
        y1=self.LUMO-fermi-low
        y2=self.LUMO-fermi+up
        plt.ylim(y1,y2)          
        plt.xlim(x[0],x[nk-1]) 
        plt.legend()
        plt.savefig(name,dpi=300)
        plt.show()
        
    
    
      
    
    
    
    
    
    
    
    
    
    
    
    def plotfatbandseg(self, fermi, projdosout,  groupfile="groups.txt", name="test.png", low=1.5, up=1.5):
        fig,a = plt.subplots()
        nk=self.nkpoints
        x=np.zeros((nk,1))
        kvector=self.kvectors
        banddata=self.eigenvalues
        shape1=np.shape(banddata)
        for i in range(nk-1):
            x[i+1]=np.abs(np.linalg.norm(kvector[i+1]-kvector[i]))+x[i]
        splist=self.read_splabel()
        ind0=splist[0]
        ind=x[ind0]
        xindex=splist[1]
        
#insert  fat coefficients        
#norm = plt.Normalize(dydx.min(), dydx.max())
        projband0=self.bandproject_group(projdosout)
        projband00=projband0[:,2:]
        
        fatnumbs=np.shape(projband00)[1]
        with open(groupfile, "r") as f:
            groupname=[]
            for line in f:
                groupname.append(line.split()[0])
        f.close()
        
        for q in range(fatnumbs):
            mat1d=projband00[:,q]
            scaledmat=np.zeros(shape1)
            for p in range(shape1[1]):
                temp=mat1d[p*shape1[0]:(p+1)*shape1[0]]
                scaledmat[:,p]=temp
                
            plt.plot(np.reshape(x,(self.nkpoints)), -100+x*0.00, label=groupname[q])
                

                
            for k in range(0, int(self.nband/self.spin)):
                x=np.reshape(x,(self.nkpoints))
                yup=np.array(banddata[self.spin*k,:])-fermi
                yupr=np.reshape(yup,(self.nkpoints))
                points = np.array([x, yupr]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                size1=100*scaledmat[self.spin*k,:]
                size1=size1[:-1]
                lc = LineCollection(segments, linewidths=0.05*size1,color=new_colors[q+3],antialiased=True)
                a.add_collection(lc)
#               plt.scatter(x,yup-fermi, s= size1, c=new_colors[q+1],  zorder=3,marker='|', edgecolors='none',label="qroup_"+str(q))
                ydn=banddata[self.spin*k+1*(self.spin-1),:]-fermi
                ydnr=np.reshape(yup,(self.nkpoints))
                points = np.array([x, ydnr]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                size2=100*scaledmat[self.spin*k+1*(self.spin-1),:]
                size2=size2[:-1]
                lc = LineCollection(segments, linewidths=0.05*size2, color=new_colors[q+3], antialiased=True)
                a.add_collection(lc)          
            
        for k in range(0, int(self.nband/self.spin)):
            yup=banddata[self.spin*k,:]
            plt.plot(x,yup-fermi,  linewidth=1. , alpha=0.4, c= new_colors[9],zorder=1)
            ydn=banddata[self.spin*k+1*(self.spin-1),:]
            plt.plot(x,ydn-fermi,  linewidth=1 , alpha=.4, color=new_colors[4], zorder=2)            
        plt.xticks(ind, xindex, fontsize =22)
        for xc in ind:
            plt.axvline(x=xc, color="gray")
        plt.yticks(fontsize =22)
        plt.axhline(y=0, linewidth=1.5, color="#4C3F64")
        plt.axhline(y=0, linewidth=0.75, color="#D0BDF2")
        y1=-low
        y2=self.LUMO-fermi+up
        plt.ylim(y1,y2)          
        plt.xlim(x[0],x[nk-1]) 
        fig.set_size_inches(8,8)
        plt.legend()
        plt.savefig(name,dpi=300)
    
    
# import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument("bandout", help="name of your band output file")
# parser.add_argument("pdosout", help="name of pdos.out file")
# parser.add_argument('-g', action='store', dest='groupname',
#                     help='groupname')
# parser.add_argument('-spin', action='store', dest='nspin', type=int, default=1,
#                     help='spin value')

# args = parser.parse_args()
# test=QElog(args.bandout, spin=args.nspin)

# if args.groupname is None:    
#     data=test.read_projout(args.pdosout)
#     test.bandproject_group(args.pdosout)
# #    fermi=test.HOMO
# #   test.plotfatband(fermi,args.pdosout)
# else:
#     test.bandproject_group(args.pdosout, groupfile=args.groupname)
    
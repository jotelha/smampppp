# coding: utf-8

# In[70]:


# reads united atoms structure from .pdb, positions from Bader analysis output ACF.dat 
# and different charges from custom made bader.csv,
# then writes extended .xyz file
import csv  
import ase
import numpy as np


# In[71]:


infile_pdb = 'system100.ortho.pdb'
ua_ase_struct = ase.io.read(infile_pdb)


# In[72]:


ua_ase_struct.center()


# In[73]:


X = np.genfromtxt('ACF.dat',skip_header=2,skip_footer=27,usecols=1)
Y = np.genfromtxt('ACF.dat',skip_header=2,skip_footer=27,usecols=2)
Z = np.genfromtxt('ACF.dat',skip_header=2,skip_footer=27,usecols=3)


# In[74]:


for i,a in enumerate(ua_ase_struct):
    a.position = np.array([X[i],Y[i],Z[i]])
    #print("Atom {}: element {}, position {}".format(i,a.symbol, a.position))


# In[75]:


ATB_charges = []
bader_charges = []
Horton_charges1=[]
Horton_charges2=[]
with open('bader.csv') as f:
        csvReader = csv.reader(f)
        next(csvReader)
        for row in csvReader:
             ATB_charges.append(float(row[2]))
             bader_charges.append(float(row[3]))
             Horton_charges1.append(float(row[4]))
             Horton_charges2.append(float(row[5]))


# In[76]:


ua_ase_struct.set_array("bader_charges",np.array(bader_charges))
ua_ase_struct.set_array("ATB_charges",np.array(ATB_charges))
ua_ase_struct.set_array("horton_charges",np.array(Horton_charges1))
ua_ase_struct.set_array("horton_charges_extra_benzene_constraints",np.array(Horton_charges2))


# In[77]:


ua_ase_struct.write("bader.xyz")

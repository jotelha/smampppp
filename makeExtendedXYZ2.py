# coding: utf-8

# In[70]:


# reads united atoms structure from .pdb, positions from Bader analysis output ACF.dat 
# and different charges from custom made bader.csv,
# then writes extended .xyz file
import csv  
import ase
import numpy as np
import ase.io

# In[71]:


infile_pdb = 'system100.pdb'
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
New_six = []
New_three = []
Bader_AA= []
Bader_UA= []
old_six = []
old_three = []
unconstrained_Horton = []
with open('chargesnh2.csv') as f:
        csvReader = csv.reader(f)
        next(csvReader)
        for row in csvReader:
             ATB_charges.append(float(row[2]))
             New_six.append(float(row[3]))
             New_three.append(float(row[4]))
             Bader_AA.append(float(row[5]))
             Bader_UA.append(float(row[6]))
             old_six.append(float(row[7]))
             old_three.append(float(row[8]))
             unconstrained_Horton.append(float(row[9]))

# In[76]:


ua_ase_struct.set_array("ATB_charges",np.array(ATB_charges))
ua_ase_struct.set_array("New_constrained_lnroh-6",np.array(New_six))
ua_ase_struct.set_array("New_constranined_lnroh-3",np.array(New_three))
ua_ase_struct.set_array("Bader_AA",np.array(Bader_AA))
ua_ase_struct.set_array("Bader_UA",np.array(Bader_UA))
ua_ase_struct.set_array("old_constrained_lnroh-6",np.array(old_six))
ua_ase_struct.set_array("old_constrained_lnroh-3",np.array(old_three))
ua_ase_struct.set_array("unconstrained_Horton",np.array(unconstrained_Horton))
# In[77]:


ua_ase_struct.write("NH2.xyz")

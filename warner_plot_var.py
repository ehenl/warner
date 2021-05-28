#!/usr/bin/env python
# coding: utf-8

# # 0. Preps

# ### 0.1 Load modules

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import netCDF4
import cmocean 
import matplotlib.colors as colors
from matplotlib import cm
import datetime


# ### 0.2 Pre-set parameters for plotting

# In[2]:


plt.rcParams.update({
    "font.weight": "normal",
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "lines.color": "k",
    "axes.titlesize": 18,
    "axes.titleweight": "normal",
    "axes.labelweight": "bold",
    "axes.labelsize": 18,
    "contour.linewidth": 0.8,
    "legend.fontsize": 14
})


# ### 0.3 Plotting functions

# In[3]:


def warner_plot(x_axis,y_axis,var,tmin,tmax,cmap,it,name,label):
    '''
    This function plots a chosen variable in a z-x-section.
    _____________________________________________________________________
    
    x_axis ---> distance along x-axis
    y_axis ---> depth
    var    ---> chosen variable (e.g. salinity)
    tmin   ---> start value of the chosen variable for its normalization
    tmax   ---> stop value of the chosen variable for its normalization
    cmap   ---> choose predefined colormap
    it     ---> choose time step
    label  ---> label of colormap
    _____________________________________________________________________
    '''
    fig, ax  = plt.subplots    (figsize=(9,6))                                                               # plot fig with one subplot
    norm_fld = colors.Normalize(vmin=tmin, vmax=tmax, clip=False)                                            # linearly normalizes data into [vmin,vmax] interval
    cs       = ax.pcolormesh   (x_axis,y_axis,var,norm=norm_fld,cmap=cmap)                                   # plot variable with pcolormesh
    cbar     = plt.colorbar    (cs,orientation='horizontal',shrink=0.75,pad=0.2)                             # add colorbar
    xticks   = np.linspace     (tmin,tmax,num=5)                                                             # specify location of colorbar ticks
    timestamp = datetime.datetime(year=2001,month=1,day=1,hour=0,minute=0,second=0)                          # start counting timesteps from 2001-01-01 00:00:00
    timestamp += datetime.timedelta(seconds=int(time.data[it]))                                              # add time steps to timestamp to make connection between real
                                                                                                             # time and time step
    cbar.set_label(label)                                                                                    # add colorbar label
    cbar.set_ticks(xticks)                                                                                   # set colorbar ticks
    ax.set(title=str(timestamp),xlabel='Distance along x [m]',ylabel='Depth [m]',
           xlim=(0, 100000),ylim=(-15,1.5))                                                                  # set axis labels and limits
    ax.tick_params(labeltop=False, labelright=False,direction = 'inout')                                     # specify axis ticks
    plt.axhline(0,color='k',lw=3,ls=':')                                                                     # add horizontal line at 0
    plt.savefig('/work/henell/tools/GETM/setups/warner/figs_warner/'+name+'_'+str(str(it).zfill(3))+'.png')  # save fig as png file with variable name and its timestep 
    plt.clf()                                                                                                # clear fig


# In[4]:


cmap_uu   = cm.get_cmap("RdBu_r")      # pre-defined colormap for zonal velocity
cmap_salt = cmocean.cm.delta           # pre-defined colormap for salinity


# # 1. Load file with resolved tidal data

# In[5]:


path_file = '/home/henell/WORK/tools/GETM/setups/warner/warner.3d.nc'    # location of file
id_file   = netCDF4.Dataset(path_file, 'r')                              # read file


# In[6]:


# find out dimensions of variables contained in file
list_var = id_file.variables.keys()
print('\n Variables in the files:')
for cvar in list_var:
    nb_dim = len(id_file.variables[cvar].dimensions)
    print(' *** '+cvar+' -> has '+str(nb_dim)+' dimensions')
print('')


# In[7]:


j = 0                                         # j is set to zero since in this case our variables do not depend on yc
xc    = id_file.variables['xc'][:]            # xc [m]
time  = id_file.variables['time'][1:]         # time [s], seconds since 2000-01-01 00:00:00
                                              # exclude initial tidal snapshot (avoid double-counting during tidal averaging)
sigma = id_file.variables['sigma'][:]         # depth in sigma coordinates
elev  = id_file.variables['elev'][:,j,:]      # elevation [m] (is equal to eta)
salt  = id_file.variables['salt'][:,:,j,:]    # salinity [PSU]
u     = id_file.variables['u'][1:,j,:]        # int. zonal vel. [m/s] 
uu    = id_file.variables['uu'][:,:,j,:]      # zonal vel. [m/s]
bathy = id_file.variables['bathymetry'][j,:]  # bathymetry [m] (is equal to H)
dx    = id_file.variables['dx'][:]            # grid spacing (x) [m]
dy    = id_file.variables['dy'][:]            # grid spacing (y) [m]
id_file.close()                               # close netCDF file


# In[8]:


(Nt,Nk,Ni) = np.shape(salt[:,:,:])
print('\n Shape of the domain Nt, Nk, Ni = ',Nt,Nk,Ni)


# # 2. Compute $dz$

# In[9]:


dz = np.empty( (Nt, Nk, Ni) )
for t in range(Nt):
    for k in range(1,Nk):
        for i in range(Ni):
            dz[t,k,i] = ( sigma[k] - sigma[k-1] ) * ( bathy[i] + elev[t,i] )


# In[10]:

# calculate zx based on dz
xx = ( xc[:-1] + xc[1:] ) / 2
zw = np.zeros_like(dz, subok=False)
zw[:,:,:] = np.cumsum(dz[:,:,:], axis=1)
zw = zw - bathy
zx = ( zw[:,:,:-1] + zw[:,:,1:] ) / 2
zx[:,:,0] = zw[:,:,1]
zx[:,:,-1] = zw[:,:,-2]


# In[11]:


xx2d = np.resize(xx, zx[0,:,:].shape)


# # 3. Plot variable(s) in warner testcase

# In[12]:


s = 0   ### start
f = 112 ### finish


# In[13]:


for i in range(s,f):
    warner_plot(xx2d,zx[i,:,:],salt[i,:,:-1],0,30,cmap_salt,i,'salt','S [gkg$^{-1}$]')
    warner_plot(xx2d,zx[i,:,:],uu[i,:,:-1],-0.5,0.5,cmap_uu,i,'uu','Zonal vel. [ms$^{-1}$]')


# In[ ]:





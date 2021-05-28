import numpy as np
import argparse
import netCDF4

length = 100e+3                         # length of estuary [m]
width  = 500.                           # width of estuary [m]
imax   = 1+200 ; dx = length / (imax-1) # number of cells (including 1 open boundary cell) and grid size
jmax   =     1 ; dy = width  / jmax     # number of cells and grid size (jmax=20 in Warner2005)
kmax   =    40                          # number of layers (20 in Warner2005)

bathy_riv =  5. # bathymetry at river end [m]
bathy_ocn = 15. # bathymetry at ocean end [m] (10. in Warner2005)
salt_riv  =  0. # salinity at river end [g kg-1]
salt_ocn  = 30. # salinity at ocean end [g kg-1]

Q_riv = 50. # river discharge [m3 s-1]

tide_period = 44800.  # tidal period [s] (correct 44714, but 44712 provides more integer divisions; 44800 is divisor of 14 days
tide_amp    =     0.6 # tidal elevation amplitude [m]

forcing_per_period = 56 ; dt_bdy = tide_period / forcing_per_period # resolution of open boundary forcing, should be a divisor of tide_period and a multiple of 4, 4*dt_bdy should be divisor of 43200,  (112 for 44800)
sim_cycles = 405 # length of simulation in number of tidal cycles (should be multiple of 27)

# tide_period |   44712 |  44800 |
# timestep    |      18 |     20 |
# nlast       | 1242000 | 907200 | = sim_cycles * meanout
# M           |      23 |     20 |
# meanout     |    2484 |   2240 | = tide_period / timestep
# hotout      |         |        | = meanout * 27
# first_3d    |         |        | = nlast - meanout


parser = argparse.ArgumentParser(description='Create input and forcing files for GETM Warner test case')
parser.add_argument("--bathy_ocn", dest="bathy_ocn", action="store", default=bathy_ocn, help='bathymetry at ocean end. Default is %f.' % (bathy_ocn) )
parser.add_argument("--tide_period", dest="tide_period", action="store", default=tide_period, help='tidal period in seconds. Default is %f.' % (tide_period) )
parser.add_argument("--sim_cycles", dest="sim_cycles", action="store", default=sim_cycles, help='simulation length in number of tidal cycles. Default is %d.' % (sim_cycles) )
args = parser.parse_args()

tide_period = float(args.tide_period)
sim_cycles = int(args.sim_cycles)

def func_bathy(y, x):
  return bathy_ocn + (bathy_riv - bathy_ocn) / length * x

def func_salt(z, y, x):
  x_salt  = 30e+3
  x_fresh = length - 20e+3

  if x < x_salt:
    salt = salt_ocn
  elif x > x_fresh:
    salt = salt_riv
  else:
    salt = salt_ocn + (salt_riv - salt_ocn) / (x_fresh - x_salt) * (x - x_salt)
  return salt

def func_elev(t):
  return tide_amp * np.sin( t / tide_period * 2 * np.pi )

nmax = 1 + sim_cycles*forcing_per_period

xx = np.linspace( -dx, length, imax+1 )
xc = ( xx[:-1] + xx[1:] ) / 2

yx = np.linspace( -width/2, width/2, jmax+1 )
yc = ( yx[:-1] + yx[1:] ) / 2

time = np.linspace( 0., sim_cycles*tide_period, nmax )

bathy = np.empty( (jmax, imax) )
for j in range(jmax):
  for i in range(imax):
    bathy[j,i] = func_bathy( yc[j], xc[i] )

# zax in initialisation must be postive downwards!
zax = bathy.max()

salt = np.empty( (1, 1, jmax, imax) )
for j in range(jmax):
  for i in range(imax):
    salt[0,0,j,i] = func_salt( -zax, yc[j], xc[i] )

elev = np.empty( (nmax, jmax) )
for n in range(nmax):
  elev[n,j] = func_elev( time[n] )

tname = 'time'
zname = 'zax'
yname = 'yc'
xname = 'xc'
gtype_name = 'grid_type'
dy_name = 'dy'
bathy_name = 'bathymetry'
salt_name = 'salt'
nbdyp_name = 'nbdyp'
elev_name = 'elev'

with open( 'warner.dim', 'w' ) as FILE:

  FILE.write('! specifies the dimensions - in case of static compilation\n')
  FILE.write('   integer, parameter :: iextr=%d,jextr=%d\n' % (imax, jmax) )
  FILE.write('   integer, parameter :: imin=1,imax=iextr,jmin=1,jmax=jextr\n')
  FILE.write('   integer, parameter :: iimin=1,iimax=imax,jjmin=1,jjmax=jmax\n')
  FILE.write('   integer, parameter :: kmax=%d\n' % (kmax) )


with open( 'bdyinfo.dat', 'w' ) as FILE:

  FILE.write('! Number of western boundaries (NWB)\n')
  FILE.write('1\n')
  FILE.write('! location in index space (wi, wfj:wlj) and boundary types\n')
  FILE.write('! wi wfj wlj bdy_2d_type bdy_3d_type\n')
  FILE.write('1 1 %d 3 0\n' % (jmax) )
  FILE.write('!\n')
  FILE.write('! Number of northern boundaries (NNB)\n')
  FILE.write('0\n')
  FILE.write('! location in index space (nj, nfi:nli) and boundary types\n')
  FILE.write('! ni nfi nli bdy_2d_type bdy_3d_type\n')
  FILE.write('!\n')
  FILE.write('! Number of eastern boundaries (NEB)\n')
  FILE.write('0\n')
  FILE.write('! location in index space (ei, efj:elj) and boundary types\n')
  FILE.write('! ei efj elj bdy_2d_type bdy_3d_type\n')
  FILE.write('!\n')
  FILE.write('! Number of southern boundaries (NSB)\n')
  FILE.write('0\n')
  FILE.write('! location in index space (sj, sfi:sli) and boundary types\n')
  FILE.write('! sj sfi sli bdy_2d_type bdy_3d_type\n')
  FILE.write('!\n')


with open( 'riverinfo.dat', 'w' ) as FILE:

  FILE.write('! Number of river entry points (nriver)\n')
  FILE.write('%d\n' % (jmax) )
  FILE.write('! location in index space (ir, jr) and river name [and optionally lower and upper river range rzl>rzu>0]\n')
  FILE.write('! ir jr river_name [rzl rzu]\n')
  for j in range(jmax):
    FILE.write('%d %d river\n' % (imax, j+1) )

with netCDF4.Dataset('warner.input.nc', mode='w') as NC:

  dims = {}
  vars = {}

  for name, dimlen in zip( (xname, yname, zname, tname), (imax, jmax, 1, 1) ):
    dims[name] = NC.createDimension( name, dimlen )
    vars[name] = NC.createVariable( name, 'f', (name,) )

  vars[xname][:] = xc
  vars[yname][:] = yc
  vars[zname][:] = zax
  vars[tname][:] = 0.0

  vars[dy_name] = NC.createVariable( dy_name, 'f')
  vars[dy_name][:] = dy

  vars[bathy_name] = NC.createVariable( bathy_name, 'f', (yname, xname) )
  vars[bathy_name][:] = bathy

  vars[salt_name ] = NC.createVariable( salt_name , 'f', (tname, zname, yname, xname) )
  vars[salt_name][:] = salt

  vars[gtype_name] = NC.createVariable( gtype_name, 'i' )
  vars[gtype_name][:] = 1


with netCDF4.Dataset('warner.bdy2d.nc',mode='w') as NC:

  dims = {}
  vars = {}

  # time dimension in bdy2d.nc must be unlimited
  for name, dimlen in zip( (tname, nbdyp_name), (None, jmax) ):
    dims[name] = NC.createDimension( name, dimlen )
    vars[name] = NC.createVariable( name, 'f', (name,) )

  vars[tname].units = 'seconds since 2000-01-01 00:00:00'
  vars[tname][:] = time

  vars[elev_name] = NC.createVariable( elev_name, 'f', (tname, nbdyp_name) )
  vars[elev_name][:] = elev


with netCDF4.Dataset('warner.bdy3d.nc',mode='w') as NC:

  dims = {}
  vars = {}

  for name, dimlen in zip( (zname, tname, nbdyp_name), (1, 2, jmax) ):
    dims[name] = NC.createDimension( name, dimlen )
    vars[name] = NC.createVariable( name, 'f', (name,) )

  vars[zname][:] = zax

  vars[tname].units = 'seconds since 2000-01-01 00:00:00'
  vars[tname][:] = ( 0., time[-1] )

  vars[salt_name] = NC.createVariable( salt_name, 'f', (tname, nbdyp_name, zname) )
  vars[salt_name][:] = salt_ocn
  vars['temp'] = NC.createVariable( 'temp', 'f', (tname, nbdyp_name, zname) )
  vars['temp'][:] = 0.


with netCDF4.Dataset('warner.rivers.nc',mode='w') as NC:

  dims = {}
  vars = {}
  # time in rivers.nc must be unlimited!
  dims[tname] = NC.createDimension( tname, None )
  vars[tname] = NC.createVariable( tname, 'f', (tname,) )
  vars[tname].units = 'seconds since 2000-01-01 00:00:00'
  vars[tname][:] = ( 0., time[-1] )

  riv_name = 'river'
  vars[riv_name] = NC.createVariable( riv_name, 'f', (tname,) )
  vars[riv_name][:] = Q_riv



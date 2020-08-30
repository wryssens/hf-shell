import numpy as np

def generate_level_densities(fname, s2=[]):

  #-----------------------------------------------------------------------------
  # The table should contain
  # 0 1 2  3   4 5 6  7   8   9   10      11     12      13  14  15  16  17  18    
  # i Z N beta E F S Q20 L20 Q22 Var(Q) ln Z_gc  ln Z_c  J2X J2Y J2Z IBX IBY IBZ  
  dat = np.genfromtxt(fname)  

  # Calculate the canonical energy and entropy
  E_c =  -findiff(dat[:,3], dat[:,11])
  S_c =           dat[:,3] * E_c + dat[:,11]

  # Heat capacity and state density
  C_c =   findiff(dat[:,3], E_c)
  rho_state = 1./(np.sqrt(2 * np.pi * np.abs(C_c))) * np.exp(S_c)

  # sigma^2
  if(len(s2) == 0):
    sigma = dat[:,16]
  else:
    sigma = s2  

  fac = 1./np.sqrt(2.0 * np.pi * sigma)
  rho_level = rho_state * fac

  return(E_c, S_c, rho_state, rho_level)

#-------------------------------------------------------------------------------
# First order finite difference on an array, not taking out points.
def findiff(x,y):
	if len(x) != len(y):
		sys.stderr.write('Error: Derivative variables have unequal length')
		sys.exit(1)
	f = np.zeros(len(y))
	for i in range(0,len(y)-1):
		if i == 0:
			dx = x[1]-x[0]
			dy = y[1]-y[0]
			f[i] = dy/dx
		else:
			dx = x[i+1] - x[i-1]
			dy = y[i+1]-y[i-1]
			f[i] = dy/dx
	f[-1] = (y[-1]-y[-2])/(x[-1]-x[-2])
	return f


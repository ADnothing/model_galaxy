from __future__ import print_function
from builtins import input
from builtins import range
import numpy as np
import matplotlib.pyplot as plt
import pypeg.pypegm as pp
import pypeg.populations as ppop
import os
import scipy
from scipy.interpolate import InterpolatedUnivariateSpline 
from scipy.interpolate import interp1d
from astropy.io import fits
from time import time
from importlib import reload 
from astropy import cosmology
from copy import deepcopy
from scipy import optimize
from scipy.optimize import minimize

import MCMC #some useful functions for the MCMC method

newcosmo = cosmology.LambdaCDM(70, 0.3, 0.7)
pp.define_cosmo(newcosmo)
filterb =  pp.Filter(filename = os.environ['ZPEG_ROOT']+'/data/filters/BJ.fil') # LF defined for this filter

def make_sdss_lf():
  sdss_lf = ppop.Lfunction(np.arange(-25.,-15.,0.1))
  M100 = -19.7 + filterb.ABVega # Vega to AB
  phi100 = 0.02
  alpha = -1.
  sdss_lf.zlf = 0.
  h = pp.cosmo_dict['cosmo'].H0.value/100.
  sdss_lf.Schechter(M100+5*np.log10(h), phi100*h**3,  alpha) 
  return sdss_lf
#sdss_lf =  make_sdss_lf()

def make_sdss_lf_Alhambra_2017(): # at z=0.5 (pivot of their model). Table 2 of López-Sanjuan et al 2017. Cosmo 0.3, 0.7, 70.
  sdss_lf = ppop.Lfunction(np.arange(-25.,-15.,0.1))
  sdss_lf_sf = ppop.Lfunction(np.arange(-25.,-15.,0.1))
  sdss_lf_q1 = ppop.Lfunction(np.arange(-25.,-15.,0.1))
  sdss_lf_q2 = ppop.Lfunction(np.arange(-25.,-15.,0.1))
  sdss_lf_q = ppop.Lfunction(np.arange(-25.,-15.,0.1))

  sdss_lf.zlf = 0.5
  sdss_lf_sf.zlf = 0.5
  sdss_lf_q.zlf = 0.5
  hlf = 0.7
  h = pp.cosmo_dict['cosmo'].H0.value/100.

  # SF
  phis = 10.**(-2.51)
  Ms = -21.00  #  AB
  alpha = -1.29
  sdss_lf_sf.Schechter(Ms-5*np.log10(hlf)+5*np.log10(h), phis*(h/hlf)**3,  alpha)

  # Q
  phis = 10.**(-2.76)
  Ms = -20.86  #  AB
  alpha = -0.53
  sdss_lf_q1.Schechter(Ms-5*np.log10(hlf)+5*np.log10(h), phis*(h/hlf)**3,  alpha)
  Mf = -17.00
  beta = -1.31
  sdss_lf_q2.Schechter(Mf-5*np.log10(hlf)+5*np.log10(h), phis*(h/hlf)**3,  beta)

  sdss_lf_q.N = sdss_lf_q1.N + sdss_lf_q2.N

  sdss_lf.N  = sdss_lf_sf.N + sdss_lf_q.N

  return sdss_lf

sdss_lf =  make_sdss_lf_Alhambra_2017()


def TBMF_from_z0LF_with_scen(model, z0LF = None): 
  """ Returns a total baryonic mass function for the model passed as a parameter
    required to match the z=0 LF """

  # LF
  if z0LF is None:
    z0LF = sdss_lf
  myTBMF = ppop.LF_to_MF(z0LF, model, filterb, z = z0LF.zlf, stellar = False, zfor = model.zfor)

  return(myTBMF)

def filename_encode(params, prefix = ''):
  return prefix+'ext_{:1d}_tinfall_{:08.2f}_tsf_{:08.2f}_twinds_{:08.2f}.fits'.\
        format(params['extinction'],params['tinfall'],params['tsf'],params['twinds'] )

def define_scen(params, ssp_prefix = None):
  """ Defines a scenario object and returns the object 

  Parameters
  ----------
  params : dictionary with the PEGASE scenario parameters: tsf, tinfall, twinds, zfor
  """

  model = pp.Model()

  if ssp_prefix is not None:
    model.SSPs_file = ssp_prefix+'_SSPs.dat'
  model.scen.stellib = 'stellibLCBcor.fits'
  model.scen.type_sf = model.scen.type_sf_dict['Schmidt']
  model.scen.p1 = 1.
  model.scen.p2 = params['tsf']
  model.scen.infall = True
  model.scen.tau_infall = params['tinfall']
  model.scen.tau_winds = params['twinds']
  model.scen.extinction = params['extinction'] # 0 = none, 1 = spheroidal, 2 = disk incl average
  model.scen.output_file = filename_encode(params, prefix = model.SSPs_file[:-4]+'_')
  model.zfor = params['zfor']
  
  return model

def plot_counts(mags, counts, *args, counts_are_in_log = False, normfunc = None, ax = None, **kwargs):
  """ normfunc is a function that normalizes the __linear__ counts"""

  if counts_are_in_log: 
    y = 10.**(counts)
  else:
    y = np.copy(counts)

  if normfunc is not None:
    y = normfunc(mags, y)

  # plot the counts, not the log(counts)
  if ax is None:
    plt.plot(mags, y, *args, **kwargs)
  else:
    ax.plot(mags, y, *args, **kwargs)
    
  return y

# DATA ------------------------------------
def define_obscounts():
  #read and calibrate filter
  prefix_filters = os.environ['ZPEG_ROOT']+'/data/filters/'
  fns = [prefix_filters + astr +'.fil' for astr in ['300','U3_BK78','JK']]+['i_subaru.res', 'K_uv.res']
  filters = [pp.Filter(filename = fn) for fn in fns]
  datafns = ['data/counts/FRV/'+ astr +'.res' for astr in ['comptF300W', 'comptU','comptB', 'comptI', 'comptK']]
  return prefix_filters, fns, datafns, filters
prefix_filters, fns_obs, datafns, filters = define_obscounts()
datas = []
for idata in range(len(filters)):
  data = np.loadtxt(datafns[idata]) # Vega magnitudes
  data[:,0] += filters[idata].ABVega # to AB
  data[:,1:4] += np.log10(2.) # to N/0.5mag/sqdeg to N/mag/sqdeg
  datas.append(data)
# DATA ------------------------------------

def reference_counts():
    # LF
  myLF_sch = sdss_lf

  # SED
  model_noevol_flat1 = pp.Model()
  model_noevol_flat1.read_from_p2file(os.environ['ZPEG_ROOT']+'/data/templates/Salp_200ages/Sb.dat', sigma = 10.)
  model_noevol_flat1.seds.fevol[:, :] = 5e35/3000./model_noevol_flat1.seds.w**1 # same SED at each timestep
  zfor = 10.
  myTBMF = ppop.LF_to_MF(myLF_sch, model_noevol_flat1, filterb, z=sdss_lf.zlf, 
                         stellar = False, zfor=zfor)
    
  filter_counts = filterb # why not?

  counts_generator = ppop.Counts()
  counts_values = counts_generator.counts(model_noevol_flat1, myTBMF, filter_counts, 
                                          zfor=zfor, extrapolate_MF = True, verbose = False)
  return counts_generator.marr, counts_values
ref_mags, ref_counts = reference_counts()

def norm_to_lambda1(mags, counts):
  global ref_mags, ref_counts
  try:
    if ref_mags is None:
      pass # just test if ref_mags is defined
  except:
    ref_mags, ref_counts = reference_counts()        
  return counts/10.**ppop.extrap(mags, ref_mags, np.log10(ref_counts))

def make_counts(myTBMF, model, filters, mags = None):
  counts_generator = ppop.Counts(marr = mags, zarr = None)
  mags = counts_generator.marr
  extrapolate_MF = True
  counts_values = []
  for ifilter in range(len(filters)):
    filter_counts = filters[ifilter]
    counts_values.append(counts_generator.counts(model, myTBMF, filter_counts,
                                                 zfor = model.zfor, extrapolate_MF = extrapolate_MF, 
                                                 verbose = False))
        
  return mags, counts_values

def compute_chi2(mod_mags, mod_counts, datas):
  """ Computes the chi2 array (1 value per filter) for the distance between the model counts and the counts observed in datafns files """

  assert len(mod_counts) == len(datafns)
  nfilters = len(mod_counts)
  chi2s = np.zeros(nfilters)
  for ifilter in range(nfilters):
    data = datas[ifilter]
    logyobs = data[:,1]
    finterp = interp1d(mod_mags, np.log10(mod_counts[ifilter]), assume_sorted=False)
    logymod = finterp(data[:,0])
    logsigmaobs = np.where(data[:,3] > data[:,1], data[:,3] - data[:,1], 0.2) # 0.2 dex by default...
    #logsigmaobs = np.where(logsigmaobs > 0.01 , logsigmaobs, 0.05) # 0.1 dex  is the minimum
    chi2s[ifilter] = np.sum(((logymod - logyobs)/logsigmaobs)**2)/len(logyobs)
  return chi2s

# define default SSPs
myssp = pp.SSP()
myssp.build()

def compute_counts(params):
  test_params = {'tinfall':10.**params[0], 'tsf': 10.**params[1], 'twinds': 10.**params[2], 'zfor': 7., 'extinction':2}
  model = define_scen(test_params, ssp_prefix = 'P2sp_'+myssp.ssp_prefix)
  pp.compute_scenarios([model.scen], skip_existing = True, verbose = False)
  model.read_from_fitsfile(model.scen.output_file, sigma = 10.)
  #model.read_from_p2file(os.environ['ZPEG_ROOT']+'/data/templates/Salp_200ages/Sb.dat', sigma = 10.)

  # fit its TMBF to SDSS z=0 LF
  myTBMF = TBMF_from_z0LF_with_scen(model, sdss_lf)

  # compute counts
  mags, counts_values = make_counts(myTBMF, model, filters)#, mags = np.arange(10,30,0.5))
  return mags, counts_values

def plot_counts_and_data(mags, counts_values):
  #normfunc = lambda x,y: y
  plt.ion()
  fig, axes = plt.subplots(2,3, figsize = (10,6))
  axes = axes.flatten()

  normfunc = norm_to_lambda1 # or to_euclidian or norm_to_lambda1
  for idata in range(len(fns_obs)):
    plt.sca(axes[idata])
    #plot data
    data = datas[idata]
    xmod = data[:,0]
    ymod = normfunc(data[:,0], 10**data[:,1])
    yerr = [normfunc(data[:,0], 10**data[:,1]) - normfunc(data[:,0], 10**data[:,2]), # lower
            normfunc(data[:,0], 10**data[:,3]) - normfunc(data[:,0], 10**data[:,1])] # upper
    
    #plt.plot(xmod, ymod, '.', c = 'red', ms = 5 ,label = 'obs : '+filter_counts.name)
    plt.errorbar(xmod, ymod, yerr, mfc = 'red', mec = 'black', ms = 7, marker = '.', ls = '', ecolor = 'green', fmt = '.', elinewidth = 0.4)

    plt.yscale('log')
    plt.ylim(np.median(ymod)*np.array([1e-1,1e1]))
    plt.xlim([10,30])
    plt.legend()

    # plot model prediction
    plot_counts(mags, normfunc(mags, counts_values[idata]), '-', linewidth = 2.)
    # plot ref model
    plt.plot(ref_mags, normfunc(ref_mags, ref_counts), '-', linewidth = 2.)
    plt.grid()
  a = input()

def plot_TBMF(myTBMF):
  plt.ion()
  plt.plot(myTBMF.logM, myTBMF.N)
  plt.xlabel('log M')
  plt.ylabel(u'N/dex/Mpc$^3$')
  plt.yscale('log')
  plt.ylim(1e-6,1)
  a = input()

def make_demo(x):    
  mags, counts_values = compute_counts(x)
  chi2s = compute_chi2(mags, counts_values, datas)
  print('CHI2 = ', chi2s)
  #plot_TBMF(myTBMF)
  plot_counts_and_data(mags, counts_values)


def make_uniform_prior(bxmin, bxmax):
    """
    Retourne une fonction qui renvoie une exception (erreur) 
    si la proposition est en dehors de [bmin, bmax]
    """
    def nchi2(parameters, *extra):
        if np.any(parameters > bxmax) or np.any(parameters < bxmin):
            raise MCMC.UndefinedChi2 #Exception
        return 0 # il suffit de retourner le chi2 à une constante près
    return nchi2

def combine_chi2(chi2_list):
    """
    Renvoie une fonction qui est la somme 
    des fonctions de log-vraisemblance fournies dans la liste chi2_list
    """
    def nchi2(prs):
        lkv = 0
        for chi2 in chi2_list:
            lkv += chi2(prs)
        return lkv
    return nchi2

def mcmc_step(chi2_func, xi, SigmaProp, extra=(), chi2_courant=None):
    """
    chi2_func calcule le chi2 de sorte que la fonction de vraisemblance soit exp(-.5*chi2)
    xi est le point actuel dans l'espace des paramètres
    SigmaProp est la matrice de covariance utilisée pour proposer un nouveau point
    """

    if chi2_courant==None:
        chi2_courant = chi2_func(xi, *extra) # calculer la vraisemblance initiale si nécessaire
 
    xf = np.random.multivariate_normal(xi, SigmaProp, 1)[0] # propose a new point

    try:
        chi2_1 = chi2_func(xf, *extra)
    except:
        # refuser proprement le point s'il est en dehors du domaine
        # de validité des paramètres
        return False, chi2_courant, xi

    accept = np.random.uniform() < np.exp(-.5 * (chi2_1-chi2_courant))

    if accept:    # on teste l'acceptance
        return True, chi2_1, xf # on accepte
    return False, chi2_courant, xi


def compute_resvec(mod_mags, mod_counts, datas):
  """ Computes the residual vector """
  liste=[]
  assert len(mod_counts) == len(datafns)
  nfilters = len(mod_counts)
  v = np.zeros(nfilters)
  for ifilter in range(nfilters):
    data = datas[ifilter]
    logyobs = data[:,1]
    finterp = interp1d(mod_mags, np.log10(mod_counts[ifilter]), assume_sorted=False)
    logymod = finterp(data[:,0])
    logsigmaobs = np.where(data[:,3] > data[:,1], data[:,3] - data[:,1], 0.2) # 0.2 dex by default...
    #logsigmaobs = np.where(logsigmaobs > 0.01 , logsigmaobs, 0.05) # 0.1 dex  is the minimum
    liste = liste+((logymod - logyobs)/logsigmaobs).tolist()
  return np.array(liste)

def resvec(params):
    """ Returns the residual vector """
    mags, counts_values = compute_counts(params)
    vec = compute_resvec(mags, counts_values, datas)
    return vec

resvec2D=lambda p2d:resvec([p2d[0], p2d[1], np.log10(20000)])

def extend(a):
    """ Converts the 2*2 input numpy array in a 3*3 numpy array 
    ex : [[1, 2], [3, 4]] => [[1, 2, 0], [3, 4, 0], [0, 0, 0]] """
    aext=np.append(a, [[0],[0]], axis=1)
    return np.append(aext, [[0, 0, 0]], axis=0)
    
def get_level(prob_arr,v):
   # Première étape, réordonner par ordre décroissant une version monodimensionnelle de la grille.
    sort_grid = prob_arr.flat[:]*1.
    sort_grid.sort() # ordre croissant
    sort_grid = sort_grid[::-1] # ordre décroissant

    # on calcule maintenant la cumulative
    cis=np.cumsum(sort_grid)
    nrm = cis[-1] # normalisation
    cis/=nrm # renormalisation
    # la fonction searchsorted() renvoie l'index ou se situe un nombre dans un tableau ordonné
    # cis est ordonné par construction
    a_pv = np.searchsorted(cis,v) 
    # a_pv est l'index pour lequel cis=v, 
    # c'est aussi l'index de la dernière valeur ajoutée pour estimer la cumulative
    levl = sort_grid[a_pv]
    return levl

def p_nsigma(n):
    return .5*(scipy.special.erf(n/np.sqrt(2))-scipy.special.erf(-n/np.sqrt(2)))

def MCMC_method(func_objective, resvec, x0, dparameters=2, steps=1500, plotting=True):
    """ MCMC_method function

	PARAMS :
	x0 : initial guess
        dparameters : number of parameters (2 for tSF and tinfall, 3 for tSF, tinfall and twinds)
        steps : number of steps used for the mcmc method
	plotting : True to visualize   """

    if dparameters==2:
        #Levenberg-Marquardt algorithm to compute Sigma
        x02D=x0[:2]
        result = scipy.optimize.least_squares(resvec2D, x02D, method="lm")
        H = result.jac.T @ result.jac
        #print(H)

        Sigma_2D = np.linalg.inv(H) #reduces H to a 2*2 array and inverse it
        Sigma=extend(Sigma_2D) #extend it 
    else :
        result = scipy.optimize.least_squares(resvec, x0, method="lm")
        H = result.jac.T @ result.jac
        Sigma = np.linalg.inv(H) #inverse H

    #numpy savez ??
    
    #MCMC method
    ch = MCMC.mcmc(mcmc_step,combine_chi2([func_objective,make_uniform_prior((0,0,0),(5,5,5))]), x0, Sigma* (2.38**2)*1./dparameters, nstep=steps)
    
    #resultats
    res=ch.mean()
    dev=ch.covariance().diagonal()**.5
    print("Acceptance",ch.acceptance())
    print("Moyenne",res)
    print("Covariance mesurée par MCMC",ch.covariance())
    print("Déviation std mesurée par MCMC",dev)

    if plotting==True:
        
        plt.figure(figsize=(5,5))
        ib,jb,gr = ch.regrid2d(0,1,20,20)
        MCMC.pretty_contour(ib,jb,gr.T,get_level,(p_nsigma(1),p_nsigma(2),p_nsigma(3)))
        plt.xlabel("$tinfall$", fontsize=20)
        plt.ylabel("$tSF$", fontsize=20)
        ch.plot_brown(0,1)
        
        if dparameters==3:
            plt.figure(figsize=(5,5))
            ib,jb,gr = ch.regrid2d(1,2,20,20)
            MCMC.pretty_contour(ib,jb,gr.T,get_level,(p_nsigma(1),p_nsigma(2),p_nsigma(3)))
            plt.xlabel("$tSF$", fontsize=20)
            plt.ylabel("$twind$", fontsize=20)
            ch.plot_brown(1,2)
            plt.figure(figsize=(5,5))
            ib,jb,gr = ch.regrid2d(0,2,20,20)
            MCMC.pretty_contour(ib,jb,gr.T,get_level,(p_nsigma(1),p_nsigma(2),p_nsigma(3)))
            plt.xlabel("$tinfall$", fontsize=20)
            plt.ylabel("$twind$", fontsize=20)
            ch.plot_brown(0,2)
        
    return res, dev

def main():

  # define a random model
  #x = np.array([0.30490286, 0.62787077, 4.30105171]) # 1pop, test_params = {'tinfall':10.**x[0], 'tsf': 10.**x[1], 'twinds': 10.**x[2], 'zfor': 7., 'extinction':2}
  #make_demo(x)
  
  #x = np.array([1.27867174, 3.61871842, 4.30105171])# 1pop, test_params = {'tinfall':10.**x[0], 'tsf': 10.**x[1], 'twinds': 10.**x[2], 'zfor': 7., 'extinction':2}
  #make_demo(x)
  
  x = np.array([3.45314666, 3.26105828, 4.30105171])# 1pop, test_params = {'tinfall':10.**x[0], 'tsf': 10.**x[1], 'twinds': 10.**x[2], 'zfor': 7., 'extinction':2}
  make_demo(x) # plot the predicted counts for this model

  def func_objective(params):
    mags, counts_values = compute_counts(params)
    chi2s = compute_chi2(mags, counts_values, datas)
    sumchi2 = np.sum(chi2s)
    return sumchi2


#r = optimize.fmin_l_bfgs_b(func_objective, x0 = np.log10([1000., 1000., 20000.]),
                             #bounds = ((0., 5.), (0., 5.), (np.log10(20001.),np.log10(20001.))), 
                             #epsilon = 0.1, approx_grad=True)
  r, ecarts = MCMC_method(func_objective, resvec,  np.log10([1000., 1000., 20000.]))

  make_demo(r) # plot the predicted counts for the best fit

  obj.writerow(r)
  for i in range(ecarts[0].size):
    obj2.writerow(ecarts[i][i])

if __name__ == "__main__":
  if os.path.isfile('data.csv') and os.path.isfile('data_err.csv'): #modifier le path enventuellement
    fichier = open('data.csv', 'a', newline='')
    fichier2 = open('data_err.csv', 'a', newline='')
    obj = csv.writer(fichier)
    obj2 = csv.writer(fichier2)
  else:
    fichier = open('data.csv', 'w', newline='')
    fichier2 = open('data_err.csv', 'w', newline='')
    obj = csv.writer(fichier)
    obj2 = csv.writer(fichier2)
    obj.writerow(np.array(['infall', 'sf', 'winds']))
    obj2.writerow(np.array(['infall', 'sf', 'winds']))

  main()

  fichier.close()
  fichier2.close()

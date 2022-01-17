from lib import *

fichier = open('param.csv', 'w', newline='')
obj = csv.writer(fichier)
obj.writerow(np.array(['T0_infall', 'T0_sf', 'T_infall', 'T_sf']))

a = np.array([[rand.uniform(0.1, 5.0), rand.uniform(0.0, 5.0), np.log10(2800), np.log10(1400)],
              [rand.uniform(0.1, 5.0), rand.uniform(0.0, 5.0), np.log10(3500), np.log10(2500)],
              [rand.uniform(0.1, 5.0), rand.uniform(0.0, 5.0), np.log10(6000), np.log10(5710)],
              [rand.uniform(0.1, 5.0), rand.uniform(0.0, 5.0), np.log10(8000), np.log10(10000)],
              [rand.uniform(0.1, 5.0), rand.uniform(0.0, 5.0), np.log10(7000), np.log10(14300)]])

for i in a:
    obj.writerow(i)
fichier.close()

for tours in range(5):
        
    os.system('python3 best_scen_MCMC.py')
        
    param = pd.read_csv('param.csv')
    os.remove('param.csv')
    param = param.drop([0])
        
    fichier = open('param.csv', 'w', newline='')
    obj = csv.writer(fichier)
    obj.writerow(np.array(['T0_infall', 'T0_sf', 'T_infall', 'T_sf']))
        
    a = np.array([param['T0_infall'].values, 
                  param['T0_sf'].values, 
                  param['T_infall'].values, 
                  param['T_sf'].values]).T
        
    for i in a:
        obj.writerow(i)
        
    fichier.close()

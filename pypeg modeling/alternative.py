from lib import *
from modeling import *
from models import *

# ---------main fonction-----------
if __name__ == "__main__":
    
# ---------Importing data and uncertainty---------
    datas = pd.read_csv('data.csv') #be sure it exist, rewrite the path just in case
    err = pd.read_csv('data_err.csv') #same for this file
    
    tab_tinfall = datas['infall'].values
    tab_tsf = datas['sf'].values
    
    tab_infall_err = err['infall'].values
    tab_sf_err = err['sf'].values
    
    
# ---------Creating Measurement array of the datas with their uncertainty---------
    tinfall = qexpy.MeasurementArray(tab_tinfall, unit="Gyr",
                                     error=tab_infall_err, name="Tinfall")
    tsf = qexpy.MeasurementArray(tab_tsf, unit="Gyr",
                                     error=tab_sf_err, name="Tsf")
    
    print('Select a model :\n')
    print('1 : linear\n')
    print('2 : polynomial\n')
    print('3 : quadratic\n')
    print('4 : exponential\n')
    n = input()
    
    if n == '1':
        mod = qexpy.FitModel.LINEAR
    elif n == '2':
        mod = qexpy.FitModel.POLYNOMIAL
    elif n == '3':
        mod =qexpy.FitModel.QUADRATIC
    elif n == '4':
        mod = qexpy.FitModel.EXPONENTIAL
    else:
        print("ERROR : No models of this format ready at the moment\n")
        
    plotQEXPY(mod, tinfall, tsf)
        
        
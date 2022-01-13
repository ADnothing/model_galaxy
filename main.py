from lib import *
from modeling import *
from models import *

# ---------main fonction-----------
if __name__ == "__main__":
    
    datas = pd.read_csv('data.csv') #be sure it exist, rewrite the path just in case
    
    tinfall = datas['infall'].values
    tsf = datas['sf'].values
    
    print('Select a model :\n')
    print('1 : linear\n')
    print('2 : parabolic\n')
    n = input()
    
    if n == '1':
        model = modlin
    elif n == '2':
        model = modparab
    else:
        print("ERROR : No models of this format ready at the moment\n")
    
    plots(curve_fit(model, tinfall, tsf), tinfall, tsf, model)
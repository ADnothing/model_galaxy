from lib import *
from models import *

def plots(ajust_args, X, Y, mod):
    """
    

    Parameters
    ----------
    ajust_args : Tuple containing : 
        - Optimized parameters in an array like shaped (n,)
        - Covariance matrice of the parameters shaped (n,n)
    X : Array like shaped (i,) X datas.
    Y : Array like shaped (i,) Y datas.
    mod : Model used to have the optimized parameters

    Returns
    -------
    None.
    Plots Y(X) datas and plot the model with the optimized parameters.
    Prints the uncertainty of all the parameters.

    """
    popt, pcov = ajust_args
    var = np.sqrt(pcov)
    x = np.arange(X.min(), X.max(), 0.1) #array used to plot the model used
    
    if mod == modlin:
        plt.scatter(X,Y)
        plt.plot(x, mod(x, popt[0], popt[1]))
        print("Parameters :")
        print("a = {:.3} (+/-) {:.3} ; b = {:.3} (+/-) {:.3}".format(popt[0], var[0][0], popt[1], var[1][1]))
    elif mod == modparab:
        plt.scatter(X,Y)
        plt.plot(x, mod(x, popt[0], popt[1], popt[2]))
        print("Parameters :")
        print("a = {:.3} (+/-) {:.3} ; b = {:.3} (+/-) {:.3} ; c = {:.3} (+/-) {:.3}".format(popt[0], var[0][0], popt[1], var[1][1], popt[2], var[2][2]))
    else:
        print("ERROR : No models of this format ready at the moment")
        

def plotQEXPY(MOD, X, Y):
    """
    

    Parameters
    ----------
    MOD : Model chosed.
    X : Measurement array of the X datas.
    Y : Measurement array of the Y datas.

    Returns
    -------
    None.
    Plots the figure with uncertainty and prints all the results of the fit

    """
    
    qplt.plot(X, Y)
    figure = qplt.get_plot()
    figure.error_bars
    resultats = figure.fit(model=MOD)
    figure.show()
from libraries import *
from binning import plot_binned
from models import *
from kde import *

def show_all():
    """
    Show the observation with the models made for each redshifts.

    """
    for i in range(10):
    
        a = np.load('z{}.npz'.format(i))
        
        plot_binned(a.f.arr_0[1][np.logical_and(a.f.arr_0[1]>-1, a.f.arr_0[1]<3)], a.f.arr_0[0][np.logical_and(a.f.arr_0[1]>-1, a.f.arr_0[1]<3)], bins = (20, 75))
        
        
        plt.plot(mod0[1][i], mod0[0][i],
                 linestyle='', marker='.', markersize=10, color='red',
                 label='Modèle MCMC')
        plt.plot(mod0bis[1][i], mod0bis[0][i],
                 linestyle='', marker='+', markersize=10, color='maroon',
                 label='Modèle MCMC bis')
        
        plt.plot(mod1[1][i], mod1[0][i],
                 linestyle='', marker='.', markersize=10, color='lawngreen',
                 label='Modèle E')
        plt.plot(mod1bis[1][i], mod1bis[0][i],
                 linestyle='', marker='+', markersize=10, color='darkolivegreen',
                 label='Modèle E bis')
        
        plt.plot(mod2[1][i], mod2[0][i],
                 linestyle='', marker='.', markersize=10, color='green',
                 label='Modèle S0')
        plt.plot(mod2bis[1][i], mod2bis[0][i],
                 linestyle='', marker='+', markersize=10, color='seagreen',
                 label='Modèle S0 bis')
        
        plt.plot(mod3[1][i], mod3[0][i],
                 linestyle='', marker='.', markersize=10, color='blue',
                 label='Modèle Sa')
        plt.plot(mod3bis[1][i], mod3bis[0][i],
                 linestyle='', marker='+', markersize=10, color='mediumpurple',
                 label='Modèle 3 bis')
        
        plt.plot(mod4[1][i], mod4[0][i],
                 linestyle='', marker='.', markersize=10, color='black',
                 label='Modèle Sb')
        plt.plot(mod4bis[1][i], mod4bis[0][i],
                 linestyle='', marker='+', markersize=10, color='grey',
                 label='Modèle Sb bis')
        
        plt.plot(mod5[1][i], mod5[0][i],
                 linestyle='', marker='.', markersize=10, color='crimson',
                 label='Modèle Sbc')
        plt.plot(mod6bis[1][i], mod6bis[0][i],
                 linestyle='', marker='+', markersize=10, color='pink',
                 label='Modèle Sbc bis')
        
        plt.plot(mod6[1][i], mod6[0][i],
                 linestyle='', marker='.', markersize=10, color='steelblue',
                 label='Modèle Sc')
        plt.plot(mod6bis[1][i], mod6bis[0][i],
                 linestyle='', marker='+', markersize=10, color='cyan',
                 label='Modèle Sc bis')
        
        plt.plot(mod7[1][i], mod7[0][i],
                 linestyle='', marker='.', markersize=10, color='orangered',
                 label='Modèle Sd')
        plt.plot(mod8bis[1][i], mod8bis[0][i],
                 linestyle='', marker='+', markersize=10, color='darkorange',
                 label='Modèle Sd bis')
        
        plt.plot(mod8[1][i], mod8[0][i],
                 linestyle='', marker='.', markersize=10, color='deeppink',
                 label='Modèle Im')
        plt.plot(mod8bis[1][i], mod8bis[0][i],
                 linestyle='', marker='+', markersize=10, color='purple',
                 label='Modèle Im bis')
        
        plt.title('Galaxies pour {} < z < {}'.format(i/2, (i+1)/2))
        plt.xlabel('B-V')
        plt.ylabel('V-I')
        plt.xlim(-1, 3)
        plt.ylim(-1, 3)
        plt.legend(bbox_to_anchor=(1.25, 1),
                             loc='upper left', borderaxespad=0.)
        plt.show()
        
def show_mod():
    """
    Show all our models for each redshifts
    """
    for i in range(10):
        plt.plot(mod0[1][i], mod0[0][i],
                 linestyle='', marker='.', markersize=10, color='red',
                 label='Modèle MCMC')
        plt.plot(mod0bis[1][i], mod0bis[0][i],
                 linestyle='', marker='+', markersize=10, color='maroon',
                 label='Modèle MCMC bis')
        
        plt.plot(mod1[1][i], mod1[0][i],
                 linestyle='', marker='.', markersize=10, color='lawngreen',
                 label='Modèle E')
        plt.plot(mod1bis[1][i], mod1bis[0][i],
                 linestyle='', marker='+', markersize=10, color='darkolivegreen',
                 label='Modèle E bis')
        
        plt.plot(mod2[1][i], mod2[0][i],
                 linestyle='', marker='.', markersize=10, color='green',
                 label='Modèle S0')
        plt.plot(mod2bis[1][i], mod2bis[0][i],
                 linestyle='', marker='+', markersize=10, color='seagreen',
                 label='Modèle S0 bis')
        
        plt.plot(mod3[1][i], mod3[0][i],
                 linestyle='', marker='.', markersize=10, color='blue',
                 label='Modèle Sa')
        plt.plot(mod3bis[1][i], mod3bis[0][i],
                 linestyle='', marker='+', markersize=10, color='mediumpurple',
                 label='Modèle 3 bis')
        
        plt.plot(mod4[1][i], mod4[0][i],
                 linestyle='', marker='.', markersize=10, color='black',
                 label='Modèle Sb')
        plt.plot(mod4bis[1][i], mod4bis[0][i],
                 linestyle='', marker='+', markersize=10, color='grey',
                 label='Modèle Sb bis')
        
        plt.plot(mod5[1][i], mod5[0][i],
                 linestyle='', marker='.', markersize=10, color='crimson',
                 label='Modèle Sbc')
        plt.plot(mod6bis[1][i], mod6bis[0][i],
                 linestyle='', marker='+', markersize=10, color='pink',
                 label='Modèle Sbc bis')
        
        plt.plot(mod6[1][i], mod6[0][i],
                 linestyle='', marker='.', markersize=10, color='steelblue',
                 label='Modèle Sc')
        plt.plot(mod6bis[1][i], mod6bis[0][i],
                 linestyle='', marker='+', markersize=10, color='cyan',
                 label='Modèle Sc bis')
        
        plt.plot(mod7[1][i], mod7[0][i],
                 linestyle='', marker='.', markersize=10, color='orangered',
                 label='Modèle Sd')
        plt.plot(mod8bis[1][i], mod8bis[0][i],
                 linestyle='', marker='+', markersize=10, color='darkorange',
                 label='Modèle Sd bis')
        
        plt.plot(mod8[1][i], mod8[0][i],
                 linestyle='', marker='.', markersize=10, color='deeppink',
                 label='Modèle Im')
        plt.plot(mod8bis[1][i], mod8bis[0][i],
                 linestyle='', marker='+', markersize=10, color='purple',
                 label='Modèle Im bis')
        
        plt.title('Galaxies pour {} < z < {}'.format(i/2, (i+1)/2))
        plt.xlabel('B-V')
        plt.ylabel('V-I')
        plt.legend(bbox_to_anchor=(1.25, 1),
                             loc='upper left', borderaxespad=0.)
        plt.show()


def show_obsv():
    """
    Display the observation for each redshifts
    """
    
    for i in range(10):
        a = np.load('z{}.npz'.format(i))
        obj = sns.jointplot(x=a.f.arr_0[1][np.logical_and(a.f.arr_0[1]>-1, a.f.arr_0[1]<3)],
                      y=a.f.arr_0[0][np.logical_and(a.f.arr_0[1]>-1, a.f.arr_0[1]<3)],
                      kind='hex')
        obj.set_axis_labels('B-V', 'V-I')
        
def show_all2():
    """
    """
    for i in range(10):
        a = np.load('z{}.npz'.format(i))
        xs = a.f.arr_0[1]
        ys = a.f.arr_0[0]
        
        extent = [np.min(xs), np.max(xs), np.min(ys), np.max(ys)]
        xv = data_coord2view_coord(xs, extent[0], extent[1])
        yv = data_coord2view_coord(ys, extent[2], extent[3])
        
        neighbours = 64
        
        ima = kNN2DDens(xv, yv, neighbours)
    
        plt.imshow(ima, origin='lower', extent=extent, cmap=cm.Blues)
        plt.title('Galaxies pour {} < z < {}'.format(i/2, (i+1)/2))
        
        plt.plot(mod0[1][i], mod0[0][i],
                 linestyle='', marker='.', markersize=10, color='red',
                 label='Modèle MCMC')
        plt.plot(mod0bis[1][i], mod0bis[0][i],
                 linestyle='', marker='+', markersize=10, color='maroon',
                 label='Modèle MCMC bis')
        
        plt.plot(mod1[1][i], mod1[0][i],
                 linestyle='', marker='.', markersize=10, color='lawngreen',
                 label='Modèle E')
        plt.plot(mod1bis[1][i], mod1bis[0][i],
                 linestyle='', marker='+', markersize=10, color='darkolivegreen',
                 label='Modèle E bis')
        
        plt.plot(mod2[1][i], mod2[0][i],
                 linestyle='', marker='.', markersize=10, color='green',
                 label='Modèle S0')
        plt.plot(mod2bis[1][i], mod2bis[0][i],
                 linestyle='', marker='+', markersize=10, color='seagreen',
                 label='Modèle S0 bis')
        
        plt.plot(mod3[1][i], mod3[0][i],
                 linestyle='', marker='.', markersize=10, color='blue',
                 label='Modèle Sa')
        plt.plot(mod3bis[1][i], mod3bis[0][i],
                 linestyle='', marker='+', markersize=10, color='mediumpurple',
                 label='Modèle 3 bis')
        
        plt.plot(mod4[1][i], mod4[0][i],
                 linestyle='', marker='.', markersize=10, color='black',
                 label='Modèle Sb')
        plt.plot(mod4bis[1][i], mod4bis[0][i],
                 linestyle='', marker='+', markersize=10, color='grey',
                 label='Modèle Sb bis')
        
        plt.plot(mod5[1][i], mod5[0][i],
                 linestyle='', marker='.', markersize=10, color='crimson',
                 label='Modèle Sbc')
        plt.plot(mod6bis[1][i], mod6bis[0][i],
                 linestyle='', marker='+', markersize=10, color='pink',
                 label='Modèle Sbc bis')
        
        plt.plot(mod6[1][i], mod6[0][i],
                 linestyle='', marker='.', markersize=10, color='steelblue',
                 label='Modèle Sc')
        plt.plot(mod6bis[1][i], mod6bis[0][i],
                 linestyle='', marker='+', markersize=10, color='cyan',
                 label='Modèle Sc bis')
        
        plt.plot(mod7[1][i], mod7[0][i],
                 linestyle='', marker='.', markersize=10, color='orangered',
                 label='Modèle Sd')
        plt.plot(mod8bis[1][i], mod8bis[0][i],
                 linestyle='', marker='+', markersize=10, color='darkorange',
                 label='Modèle Sd bis')
        
        plt.plot(mod8[1][i], mod8[0][i],
                 linestyle='', marker='.', markersize=10, color='deeppink',
                 label='Modèle Im')
        plt.plot(mod8bis[1][i], mod8bis[0][i],
                 linestyle='', marker='+', markersize=10, color='purple',
                 label='Modèle Im bis')
        
        plt.plot(mod02pop1[1][i], mod02pop1[0][i],
                 linestyle='', marker='.', markersize=10, color='gold',
                 label='Modèle MCMC 2 pop (pop1)')
        plt.plot(mod02pop1bis[1][i], mod02pop1bis[0][i],
                 linestyle='', marker='+', markersize=10, color='y',
                 label='Modèle MCMC 2 pop (pop1) bis')
        plt.plot(mod02pop2[1][i], mod02pop2[0][i],
                 linestyle='', marker='.', markersize=10, color='fuchsia',
                 label='Modèle MCMC 2 pop (pop2)')
        plt.plot(mod02pop2bis[1][i], mod02pop2bis[0][i],
                 linestyle='', marker='+', markersize=10, color='magenta',
                 label='Modèle MCMC 2 pop (pop2) bis')
        
        plt.colorbar()
        plt.xlim(-1,3)
        plt.ylim(-1,3)
        plt.xlabel('B-V')
        plt.ylabel('V-I')
        plt.legend(bbox_to_anchor=(1.25, 1),
                             loc='upper left', borderaxespad=0.)
        plt.show()

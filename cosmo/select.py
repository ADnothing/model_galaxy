from libraries import *

mask = np.array(['SC_V_MAG_AUTO', 'SC_B_MAG_AUTO', 'SC_IB827_MAG_AUTO', 'ez_z_phot'])
#ref filters :
#SC_V -> 5470.20 A
#SC_B -> 4448.00 A
#SC_IB827 -> 8243.85 A
#http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=Subaru/Suprime.IB827&&mode=browse&gname=Subaru&gname2=Suprime

hdul = fits.open('COSMOS2020_CLASSIC_R1_v2.0.fits') 
data = hdul[1].data
hdul.close()

V = data.field(mask[0])
B = data.field(mask[1])
I = data.field(mask[2])
Z = data.field(mask[3])

VI = V-I
BV = B-V

zs = np.array([Z<0.5,
               np.logical_and(Z>0.5, Z<1),
               np.logical_and(Z>1, Z<1.5),
               np.logical_and(Z>1.5, Z<2),
               np.logical_and(Z>2, Z<2.5),
               np.logical_and(Z>2.5, Z<3),
               np.logical_and(Z>3, Z<3.5),
               np.logical_and(Z>3.5, Z<4),
               np.logical_and(Z>4, Z<4.5),
               np.logical_and(Z>4.5, Z<5)])

for  i in range(10):
    vi = VI[zs[i]]
    bv = BV[zs[i]]
    maskf = np.logical_and(np.logical_not(np.isnan(vi)), np.logical_not(np.isnan(bv)))
    VItrue = vi[maskf]
    BVtrue = bv[maskf]
    np.savez('z{}.npz'.format(i), np.array([VItrue, BVtrue]))
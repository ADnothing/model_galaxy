from lib import *

#asking about the data suppression
if os.path.isfile('data.csv') and os.path.isfile('data_err.csv'):
    print('supress data.csv and data_err.csv ? (y/n)')
    a = input()
else:
    pass

if a == 'y':
    os.remove('data.csv')
    os.remove('data_err.csv')
    

print('Number of points')
n = input()
for i in range(int(n)):
    os.system('python3 best_scen_MCMC')
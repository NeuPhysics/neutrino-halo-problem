from numpy import genfromtxt, savetxt

filename = 'halo_sim_osc_100_1000_1.000000_0.001000.csv'

# print('cmake-build-debug/' + filename)

data = genfromtxt('cmake-build-debug/' + filename, delimiter=', ')
# data = genfromtxt('gxx/1/' + filename, delimiter=', ')

savetxt('data/transposed'+ filename, data.T, delimiter=', ')

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
plt.style.use('seaborn-poster')
lines = ["-","--","-.",":", "o", "v", "d", "^"]

# path = sys.argv

lbls = ['incline-alpha-0.0','original-code']
# folder = str(path[1])
folder = '../gxx/neutrino-headon-incline/'
# filename = str(path[2])
# filename = 'damping-method-mu-0.1-refl-0.01-'+lbls[0] + lbls[1]
# filename = 'damping-method-mu-0.1-refl-0.1-'+lbls[0] + lbls[1]
# filename = 'damping-method-mu-1-refl-0.1-'+lbls[0] + '-' + lbls[1]
# filename = 'damping-method-mu-1-refl-0.01-'+lbls[0] + '-' + lbls[1]
filename = 'incline-method-range-1-mu-5-refl-0.2-'+lbls[0] + '-' + lbls[1]

# datarawlist = [np.genfromtxt('../gxx/neutrino-headon-avg/halo_parallel_avg0.5_mu_0.1_REFL0.010000_ITER100000_STEPS1000_RANGE1.000000_TH_20_t2017-11-7-11-54-42.csv', delimiter = ", "), np.genfromtxt('../gxx/neutrino-headon-avg/halo_parallel_noavg_mu_0.1_REFL0.010000_ITER100000_STEPS1000_RANGE1.000000_TH_20_t2017-11-7-11-58-26.csv', delimiter = ", ")]
# datarawlist = [np.genfromtxt('../gxx/neutrino-headon-avg/halo_parallel_avg0.5_mu_0.1_REFL0.100000_ITER100000_STEPS1000_RANGE1.000000_TH_20_t2017-11-7-11-55-54.csv', delimiter = ", "), np.genfromtxt('../gxx/neutrino-headon-avg/halo_parallel_noavg_mu_0.1_REFL0.100000_ITER100000_STEPS1000_RANGE1.000000_TH_20_t2017-11-7-11-58-47.csv', delimiter = ", ")]
# datarawlist = [np.genfromtxt('../gxx/neutrino-headon-avg/halo_parallel_avg0.5_mu_1.0_REFL0.100000_ITER100000_STEPS1000_RANGE1.000000_TH_20_t2017-11-7-11-56-8.csv', delimiter = ", "), np.genfromtxt('../gxx/neutrino-headon-avg/halo_parallel_noavg_mu_1.0_REFL0.100000_ITER100000_STEPS1000_RANGE1.000000_TH_20_t2017-11-7-11-58-57.csv', delimiter = ", ")]
# datarawlist = [np.genfromtxt('../gxx/neutrino-headon-avg/halo_parallel_avg0.5_mu_1.0_REFL0.010000_ITER100000_STEPS1000_RANGE1.000000_TH_20_t2017-11-7-11-54-55.csv', delimiter = ", "), np.genfromtxt('../gxx/neutrino-headon-avg/halo_parallel_noavg_mu_1.0_REFL0.010000_ITER100000_STEPS1000_RANGE1.000000_TH_20_t2017-11-7-11-58-36.csv', delimiter = ", ")]
#datarawlist = [np.genfromtxt('../gxx/neutrino-headon-avg/neutrino-headon-avg_MU_1.000000_REFL0.200000_ITER10000000_STEPS50000_RANGE5.000000_TH_20_t2017-11-7-17-31-19.csv', delimiter = ", "), np.genfromtxt('../gxx/neutrino-headon/halo_parallel_REFL0.200000_ITER10000000_STEPS50000_RANGE5.000000_TH_20_t2017-11-5-17-59-36.csv', delimiter = ", ")]
#datarawlist = [np.genfromtxt('../gxx/neutrino-headon-avg/varyalpha-bahcall/neutrino-headon-avg-0.1_MU_1.000000_REFL0.200000_ITER10000000_STEPS50000_RANGE5.000000_TH_22_t2017-11-9-2-39-2.csv', delimiter = ", "), np.genfromtxt('../gxx/neutrino-headon-avg/varyalpha-bahcall/neutrino-headon-avg-0.0_MU_1.000000_REFL0.200000_ITER10000000_STEPS50000_RANGE5.000000_TH_22_t2017-11-8-18-29-56.csv', delimiter = ", ")]
datarawlist = [np.genfromtxt('../gxx/neutrino-headon-incline/neutrino-headon-incline-0.00_MU_5.000000_REFL0.100000_ITER100000_STEPS1000_RANGE1.000000_TH_20_t2017-11-11-15-58-38.csv', delimiter = ", "), np.genfromtxt('../gxx/neutrino-headon-incline/neutrino-headon_MU_5.000000_REFL0.100000_ITER100000_STEPS1000_RANGE1.000000_TH_20_t2017-11-11-15-59-1.csv', delimiter = ", ")]
#data = np.genfromtxt('gxx/1/halo_sim_osc_20000_100000_1.000000_0.000010.csv', delimiter = ", ")
# dataraw = np.genfromtxt( folder + filename , delimiter = ", ")
print("Data Loaded")





#for i in np.arange(1, 1000):
#    plt.plot( data[0], data[-2*i+1][::-1], label="backward "+str(i) +"th" )
#    plt.plot( data[0], data[-2*i], label="forward "+str(i)+"th" )
#plt.legend()
#plt.show()

# Clean up data: grab the first element of state

def extract(dataraw):

    datafull = dataraw
    datafull = np.delete(datafull,0, axis=0)
    data = datafull[::3]
    data = np.insert(data, 0 ,dataraw[0], axis=0)

    return data

datalist = np.array([extract(datarawlist[0]),extract(datarawlist[1])])

# data = dataraw

totlen = len(datalist[0])
length = int(len(datalist[0][0])/2)
totlenref = len(datalist[-1])

fig, ax = plt.subplots()

x = datalist[0][0][:length+1]
line1, = ax.plot(x, datalist[0][1][:length+1]/2+0.5, '.r', label=lbls[0]+' forward')
line2, = ax.plot(x[1:length+1], datalist[0][1][length:2*length][::-1]/2+0.5, '.k', label=lbls[0]+' backward')
line3, = ax.plot(x, datalist[1][1][:length+1]/2+0.5, '-r', label=lbls[1] + ' forward')
line4, = ax.plot(x[1:length+1], datalist[1][1][length:2*length][::-1]/2+0.5, '-k', label=lbls[1] + ' backward')
time_text = ax.text(0.90, 0.95, '', transform=ax.transAxes)
ax.set_title('Survival Probability')


def animate(i):
    i = int(i)
    line1.set_ydata(datalist[0][i][:length+1]/2+0.5)  # update the data for forward
    line2.set_ydata(datalist[0][i][length:2*length][::-1]/2+0.5)  # update the data for backward
    if i < totlenref:
        line3.set_ydata(datalist[1][i][:length+1]/2+0.5)  # update the data for forward
        line4.set_ydata(datalist[1][i][length:2*length][::-1]/2+0.5)  # update the data for backward
    else:
        line3.set_ydata(datalist[1][-1][:length+1]/2+0.5)  # update the data for forward
        line4.set_ydata(datalist[1][-1][length:2*length][::-1]/2+0.5)  # update the data for backward
    time_text.set_text('time = %.1f' % i)
    return line1, line2, line3, line4, time_text

def init():
    #line.set_ydata(np.ma.array(x, mask=True))
    line1.set_ydata(datalist[0][1][:length+1]/2+0.5)  # update the data for backward
    # line2.set_ydata(data[2]/2+0.5)  # update the data for backward
    time_text.set_text('')
    return line1, line2, line3, line4, time_text


ani = animation.FuncAnimation(fig, animate, np.arange(1, totlen-2, 2), init_func=init, interval=200, blit=True)



plt.ylim(0,1)
plt.xlabel('$x/\omega$')
plt.ylabel('Probability')
plt.legend(loc=8, ncol=2) #, bbox_to_anchor=(1.05, -0.1))



Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
print("Export MP4 Started")

ani.save(folder + filename + ".mp4", writer=writer)

print("Export MP4 Finished")

plt.show()


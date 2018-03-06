import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
plt.style.use('seaborn-poster')
lines = ["-","--","-.",":", "o", "v", "d", "^"]

# path = sys.argv

# folder = str(path[1])
folder = '../gxx/neutrino-headon-avg/'
# filename = str(path[2])
filename = 'Damping Method'

datarawlist = [np.genfromtxt('../gxx/omp/halo_parallel_REFL0.100000_ITER20000000_STEPS200000_RANGE20.000000_TH_20_t2017-10-12-16-8-6.csv', delimiter = ", "), np.genfromtxt('../gxx/omp/halo_parallel_REFL0.100000_ITER2000000_STEPS20000_RANGE20.000000_TH_20_t2017-10-12-15-8-40.csv', delimiter = ", ")]
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

fig, ax = plt.subplots()

x = datalist[0][0][:length+1]
line1, = ax.plot(x, datalist[0][1][:length+1]/2+0.5, '--r', label='1e-4 forward')
line2, = ax.plot(x[1:length+1], datalist[0][1][length:2*length][::-1]/2+0.5, '-k', label='1e-4 backward')
line3, = ax.plot(x, datalist[1][1][:length+1]/2+0.5, '-.r', label='1e-3 forward')
line4, = ax.plot(x[1:length+1], datalist[1][1][length:2*length][::-1]/2+0.5, ':k', label='1e-3 backward')
time_text = ax.text(0.90, 0.95, '', transform=ax.transAxes)
ax.set_title('Survival Probability')


def animate(i):
    i = int(i)
    line1.set_ydata(datalist[0][i][:length+1]/2+0.5)  # update the data for forward
    line2.set_ydata(datalist[0][i][length:2*length][::-1]/2+0.5)  # update the data for backward
    line3.set_ydata(datalist[1][i][:length+1]/2+0.5)  # update the data for forward
    line4.set_ydata(datalist[1][i][length:2*length][::-1]/2+0.5)  # update the data for backward
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
plt.legend(loc=0, ncol=4, bbox_to_anchor=(1.05, -0.1))



# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
# print("Export MP4 Started")

# ani.save(folder + filename + ".mp4", writer=writer)

# print("Export MP4 Finished")

plt.show()


import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
plt.style.use('seaborn-poster')
lines = ["-","--","-.",":", "o", "v", "d", "^"]

path = sys.argv

folder = str(path[1])
filename = str(path[2])

#data = np.genfromtxt('halo_osc_10000000.csv', delimiter = ", ")
#data = np.genfromtxt('halo_osc.csv', delimiter = ", ")
#data = np.genfromtxt('gxx/halo_osc_5.000000_10000000_2000000.csv', delimiter = ", ")
#data = np.genfromtxt('gxx/halo_osc_1.000000_10000000_200000.csv', delimiter = ", ")
#data = np.genfromtxt('gxx/1/halo_sim_osc_20000_100000_1.000000_0.000010.csv', delimiter = ", ")
dataraw = np.genfromtxt( folder + filename , delimiter = ", ")
print("Data Loaded")





#for i in np.arange(1, 1000):
#    plt.plot( data[0], data[-2*i+1][::-1], label="backward "+str(i) +"th" )
#    plt.plot( data[0], data[-2*i], label="forward "+str(i)+"th" )
#plt.legend()
#plt.show()

# # Clean up data: grab the first element of state
# datafull = dataraw
# datafull = np.delete(datafull,0, axis=0)
# data = datafull[::3]
# data = np.insert(data, 0 ,dataraw[0], axis=0)

data = dataraw


totlen = len(data)
length = int(len(data[0])/2)

fig, ax = plt.subplots()

x = data[0][:length+1]
line1, = ax.plot(x, data[1][:length+1]/2+0.5, '--r', label='forward')
line2, = ax.plot(x[1:length+1], data[1][length:2*length][::-1]/2+0.5, '-k', label='backward')
time_text = ax.text(0.90, 0.95, '', transform=ax.transAxes)
ax.set_title('Forward Beam Survival Probability')


def animate(i):
    i = int(i)
    line1.set_ydata(data[i][:length+1]/2+0.5)  # update the data for forward
    line2.set_ydata(data[i][length:2*length][::-1]/2+0.5)  # update the data for backward
    time_text.set_text('time = %.1f' % i)
    return line1, line2, time_text

def init():
    #line.set_ydata(np.ma.array(x, mask=True))
    line1.set_ydata(data[1][:length+1]/2+0.5)  # update the data for backward
    # line2.set_ydata(data[2]/2+0.5)  # update the data for backward
    time_text.set_text('')
    return line1, line2, time_text


ani = animation.FuncAnimation(fig, animate, np.arange(1, totlen-2, 2), init_func=init, interval=200, blit=True)



plt.ylim([0,1])
plt.legend(loc=8, ncol=2)



Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
print("Export MP4 Started")

ani.save(folder + filename + ".mp4", writer=writer)

print("Export MP4 Finished")

plt.show()


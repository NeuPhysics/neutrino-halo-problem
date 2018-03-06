import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
plt.style.use('seaborn-poster')
lines = ["-","--","-.",":", "o", "v", "d", "^"]

#path = sys.argv

#folder = str(path[1])
#filename = str(path[2])

#data = np.genfromtxt('halo_osc_10000000.csv', delimiter = ", ")
#data = np.genfromtxt('halo_osc.csv', delimiter = ", ")
#data = np.genfromtxt('gxx/halo_osc_5.000000_10000000_2000000.csv', delimiter = ", ")
#data = np.genfromtxt('gxx/halo_osc_1.000000_10000000_200000.csv', delimiter = ", ")
#data = np.genfromtxt('gxx/1/halo_sim_osc_20000_100000_1.000000_0.000010.csv', delimiter = ", ")
#data = np.genfromtxt( folder + filename , delimiter = ", ")

#data1 = np.genfromtxt('gxx/1e4/halo_sim_osc_20000_10000_1.000000_0.000100.csv', delimiter = ", ")
data1 = np.genfromtxt('gxx/1e6/halo_sim_osc_20000_1000000_1.000000_0.000001.csv', delimiter = ", ")
data2 = np.genfromtxt('gxx/1e5/halo_sim_osc_20000_100000_1.000000_0.000010.csv', delimiter = ", ")
print("Data Loaded")

length1 = len(data1)
length2 = len(data2)
#for i in np.arange(1, 1000):
#    plt.plot( data[0], data[-2*i+1][::-1], label="backward "+str(i) +"th" )
#    plt.plot( data[0], data[-2*i], label="forward "+str(i)+"th" )
#plt.legend()
#plt.show()


fig, ax = plt.subplots()

x1 = data1[0]
x2 = data2[0]
line11, = ax.plot(x1, data1[1], '--r', label='forward 1e6')
line12, = ax.plot(x1, data1[2], '-k', label='backward 1e6')
line21, = ax.plot(x2, data2[1], '-.r',  markersize=3, label='forward 1e5')
line22, = ax.plot(x2, data2[2], ':k',  markersize=3, label='backward 1e5')
time_text = ax.text(0.90, 0.95, '', transform=ax.transAxes)
ax.set_title('Forward Beam Survival Probability')


def animate(i):
    i = int(i)
    line11.set_ydata(data1[i]/2+0.5)  # update the data for forward
    line12.set_ydata(data1[i+1]/2+0.5)  # update the data for backward
    line21.set_ydata(data2[i]/2+0.5)  # update the data for forward
    line22.set_ydata(data2[i+1]/2+0.5)  # update the data for backward
    timetrack = int(i/2);
    time_text.set_text('time = %.1f' % timetrack)
    return line11, line12, line21, line22, time_text

def init():
    #line.set_ydata(np.ma.array(x, mask=True))
    line11.set_ydata(data1[1]/2+0.5)  # update the data for backward
    line21.set_ydata(data2[1]/2+0.5)  # update the data for backward
    # line2.set_ydata(data[2]/2+0.5)  # update the data for backward
    time_text.set_text('')
    return line11, line12, line21, line22, time_text


ani = animation.FuncAnimation(fig, animate, np.arange(1, length1-2, 2), init_func=init, interval=200, blit=True)



plt.ylim([0,1])
plt.legend(loc=8, ncol=2)



Writer = animation.writers['ffmpeg']
writer = Writer(fps=1, metadata=dict(artist='Me'), bitrate=1800)
print("Export MP4 Started")

# ani.save(folder + filename + ".mp4", writer=writer)
ani.save("gxx/convergence-1e4-5" + ".mp4", writer=writer)

print("Export MP4 Finished")

plt.show()


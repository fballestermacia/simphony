import numpy as np
import matplotlib.pyplot as plt


def readbulkekgnu(prefix):
    filelines = [line for line in open(prefix+'.gnu') if line.strip()]
    for line in filelines:
        if line.split()[0] == 'set' and line.split()[1] == 'xtics':
            labeltics = line.split('(')[1].split(')')[0].split(',')
            labels = []
            ticks = []
            for lt in labeltics:
                labels.append(lt.split('    ')[0].replace('"','').strip())
                ticks.append(float(lt.split('    ')[1]))
            
            for i in range(len(labels)):
                if labels[i] == 'G':
                    labels[i] = '$\\Gamma$'
    data = np.transpose(np.loadtxt(prefix+'.dat'))
    
    splits = np.flatnonzero(np.diff(data[0])<0)
    qpoints = np.split(data[0],splits+1)
    bands = data[1].reshape(np.shape(qpoints))
    return qpoints, bands, labels, ticks
    

factor = 0.123983#*0.24180
thztomev = 1#4.15665538536

qpoints, tbbandsperline, labels, ticks = readbulkekgnu('examples/PbTe/bulkek')

tbbandsperline *= thztomev



plt.figure()
topocolors = 'k'*8 + 'g' +'r'+'g'+'r'+'k'*2+'k'+'b'*3+'b'*3

#topocolors = 'g'*10+'g'*4+'k'*12+'g'*6+'k'*4 # For AgP2
ax1 = plt.axes()

for i, tbands in enumerate(tbbandsperline):
    ax1.plot(qpoints[i], tbands,linewidth = 3, alpha = 1, color = topocolors[i])

plt.ylabel("Frequency (ThZ)", fontsize=20 )#(cm$^{-1}$)")
plt.yticks(fontsize=20 )
plt.xlim(np.min(qpoints), np.max(qpoints))
plt.xticks(ticks=ticks, labels=labels, fontsize=20 )
for pos in ticks:
    ax1.axvline(x=pos, linewidth=0.5, color='k' )

ax1.axhline(y=0, linewidth=0.5, color='b', linestyle='--')
plt.title('PbTe with long-range', fontsize=20 )

#plt.ylim(0, )
plt.show()
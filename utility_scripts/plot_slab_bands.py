import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns


def readslabek(prefix):
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
    bottom = data[2].reshape(np.shape(qpoints))
    top = data[3].reshape(np.shape(qpoints))
    return qpoints, bands, labels, ticks, bottom, top
   

factor = 0.123983#*0.24180
thztomev = 4.15665538536

qpoints, tbbandsperline, labels, ticks, bottom, top = readslabek('examples/PbTe/slabek')

tbbandsperline *= thztomev



plt.figure(figsize=(10,5))

colors = top - bottom

ax1 = plt.axes()
hue_neg, hue_pos = 250, 15
colormap = sns.diverging_palette(hue_neg, hue_pos, center="dark", as_cmap=True,s=100,sep=100)

sc = ax1.scatter(np.array(qpoints).flatten(), np.array(tbbandsperline).flatten(), c = colors.flatten(), label='Phonnier',s=2,cmap = colormap, vmin=-1, vmax=1)

cb = plt.colorbar(sc,ticks=[-1, 0, 1])
cb.set_ticklabels(['Bottom','Bulk','Top'], fontsize = 20)

plt.yticks(fontsize=20)
plt.ylabel("Frequency (meV)",fontsize=20 )#(cm$^{-1}$)")
plt.xlim(np.min(qpoints), np.max(qpoints))
plt.xticks(ticks=ticks, labels=labels,fontsize=20)
for pos in ticks:
    ax1.axvline(x=pos, linewidth=0.5, color='k')

ax1.axhline(y=0, linewidth=0.5, color='b', linestyle='--')
plt.title('PbTe, (001) surface', fontsize=20)

#plt.ylim(17, 23.2)
#plt.savefig('DFTpycodes/plotsWT/testfig.jpeg')
plt.show()
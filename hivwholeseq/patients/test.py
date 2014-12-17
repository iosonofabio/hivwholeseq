# vim: fdm=marker
'''
author:     Fabio Zanini
date:       16/12/14
content:    Test script for quick things.
'''
import matplotlib.pyplot as plt

# Legend: 0 = conserved, 1 = polymorphic, 2 = sweeping
site_class = np.zeros((afts.shape[0], afts.shape[1]), int)
for i, aft in enumerate(afts):
    for j, aftpos in enumerate(aft):
        aftpos = np.vstack(aftpos)
        ia0 = aftpos[:, 0].argmax()
        if set((aftpos > 0.5).nonzero()[0]) - set([ia0]):
            site_class[i, j] = 2

        # FIXME: make this dependent on depth
        elif set((aftpos > 0.01).nonzero()[0]) - set([ia0]):
            site_class[i, j] = 1


# Plot as image
from scipy.cluster.hierarchy import linkage, dendrogram
fig, axs = plt.subplots(1, 2, figsize=(16, 5), gridspec_kw={'width_ratios': [1, 4]})
d = np.zeros((site_class.shape[0], site_class.shape[0]), float)
for i in xrange(d.shape[0]):
    for j in xrange(d.shape[0]):
        d[i, j] = ((site_class[i] != 0) != (site_class[j] != 0)).mean()
Z = linkage(d)
dg = dendrogram(Z, orientation='right', ax=axs[0], leaf_label_func=lambda x: '')
axs[0].set_xticklabels('')
ax = axs[1]
ax.imshow(site_class[dg['leaves']][::-1], interpolation='nearest', aspect='auto')
ax.set_yticks(np.arange(len(patients)))
ax.set_yticklabels(patients.index[dg['leaves']][::-1].tolist())
fig.suptitle(region)

plt.tight_layout(rect=(0, 0, 1, 0.96))
plt.ion()
plt.show()


is_sweeping = site_class == 2
h = np.bincount(is_sweeping.sum(axis=0))
fig, ax = plt.subplots()
plt.plot(np.arange(len(h)), h, lw=2, label='Data, '+region)
from scipy.stats import poisson
mu = curve_fit(poisson.pmf, np.arange(len(h)), 1.0 * h / h.sum(), p0=1)[0][0]
hth = poisson.pmf(np.arange(len(h)), mu)
plt.plot(np.arange(len(hth)), hth * h.sum(), lw=2, c='k',
         label='mu = '+'{:1.1e}'.format(mu)+', '+region)

plt.grid()
plt.ylabel('# sweeping sites')
plt.xlabel('# patients')
plt.yscale('log')
plt.legend(loc=1)


plt.ion()
plt.show()

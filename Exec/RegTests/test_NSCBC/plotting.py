import yt
import os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm, Normalize
from matplotlib.cm import ScalarMappable
import numpy as np

case_folder = './Exec/RegTests/test_NSCBC/'
opath = os.path.join(case_folder, "0images", )
if not os.path.exists(opath):
    os.makedirs(opath)

res = [128, 128, 4]


num = []
for name in os.listdir(case_folder):
    path = os.path.join(case_folder, name)
    if os.path.isdir(path) and name[:3]=='plt' and (not 'old' in name):
        num.append(int(name[3:]))
num = np.sort(num)
files = [case_folder + f'plt{n:05}' for n in num]
print(len(files))

# file = os.path.join(case_folder, 'plt00100')
for i, file in enumerate(files):
    ds = yt.load(file)
    t = ds.current_time
    L = (ds.domain_right_edge - ds.domain_left_edge).d

    s = yt.SlicePlot(ds, "z", "pressure")
    s.annotate_contour("pressure")
    s.save(os.path.join(opath, f"pressure{i:05}.png"))

# xy = np.meshgrid(np.array(ds.domain_right_edge[0], )
# slc = ds.slice("z", 0)
# data = slc.to_frb(width=L[0], height=L[1], resolution=res[:-1])
# array = np.array(data['pressure'])

# fig, ax = plt.subplots(1, 1, sharex=True, figsize=(5, 5))
# plot = ax.imshow(array, origin='lower', extent = (-0.065, 0.065, -0.065, 0.065))
# ax.text(4.2, 2.6, f't = {float(t)*1e3:.04} ms')
# divider = make_axes_locatable(ax)
# cax = divider.append_axes('right', size='10%', pad=0.45)
# norm = Normalize(vmin=np.min(array), vmax=np.max(array))
# plt.colorbar(ScalarMappable(norm=norm,cmap=plt.cm.viridis), cax=cax, orientation='vertical')
# fig.savefig(os.path.join(opath, 'pressure'))

from dynspec_tools import dedisperse as dd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

with open("../1413381294-I.yaml", 'r') as yaml_file:
    params = dd.parse_yaml(yaml_file)
    params['input'] = "../" + params['input']

dynspec = dd.Dynspec(**params)

fig, ax = plt.subplots(figsize=(12,5))

dynspec.plot(ax, fscrunch=64)
ax.set_xlim([1413381450, 1413381720])
plt.tight_layout()
plt.savefig('dm_detail.png')

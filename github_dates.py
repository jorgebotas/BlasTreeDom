import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
from datetime import datetime

# In case the above fails, e.g. because of missing internet connection
# use the following lists as fallback.
names = ['v2.2.4', 'v3.0.3', 'v3.0.2', 'v3.0.1', 'v3.0.0', 'v2.2.3',
            'v2.2.2', 'v2.2.1', 'v2.2.0', 'v2.1.2', 'v2.1.1', 'v2.1.0',
            'v2.0.2', 'v2.0.1', 'v2.0.0', 'v1.5.3', 'v1.5.2', 'v1.5.1',
            'v1.5.0', 'v1.4.3', 'v1.4.2', 'v1.4.1', 'v1.4.0']

dates = ['2019-02-26', '2019-02-26', '2018-11-10', '2018-11-10',
            '2018-09-18', '2018-08-10', '2018-03-17', '2018-03-16',
            '2018-03-06', '2018-01-18', '2017-12-10', '2017-10-07',
            '2017-05-10', '2017-05-02', '2017-01-17', '2016-09-09',
            '2016-07-03', '2016-01-10', '2015-10-29', '2015-02-16',
            '2014-10-26', '2014-10-18', '2014-08-26']

names = ['ASN_GLYCOSYLATION', 'ASN_GLYCOSYLATION', 'ASN_GLYCOSYLATION', 'ASN_GLYCOSYLATION', 'CAMP_PHOSPHO_SITE', 'CAMP_PHOSPHO_SITE', 'CAMP_PHOSPHO_SITE', 'CAMP_PHOSPHO_SITE', 'CAMP_PHOSPHO_SITE', 'CAMP_PHOSPHO_SITE', 'PKC_PHOSPHO_SITE', 'PKC_PHOSPHO_SITE', 'PKC_PHOSPHO_SITE', 'PKC_PHOSPHO_SITE', 'PKC_PHOSPHO_SITE', 'PKC_PHOSPHO_SITE', 'PKC_PHOSPHO_SITE', 'PKC_PHOSPHO_SITE', 'PKC_PHOSPHO_SITE', 'PKC_PHOSPHO_SITE', 'PKC_PHOSPHO_SITE', 'PKC_PHOSPHO_SITE', 'PKC_PHOSPHO_SITE', 'PKC_PHOSPHO_SITE', 'PKC_PHOSPHO_SITE', 'PKC_PHOSPHO_SITE', 'PKC_PHOSPHO_SITE', 'CK2_PHOSPHO_SITE', 'CK2_PHOSPHO_SITE', 'CK2_PHOSPHO_SITE', 'CK2_PHOSPHO_SITE', 'CK2_PHOSPHO_SITE', 'CK2_PHOSPHO_SITE', 'CK2_PHOSPHO_SITE', 'CK2_PHOSPHO_SITE', 'CK2_PHOSPHO_SITE', 'CK2_PHOSPHO_SITE', 'CK2_PHOSPHO_SITE', 'CK2_PHOSPHO_SITE', 'CK2_PHOSPHO_SITE', 'TYR_PHOSPHO_SITE_1', 'MYRISTYL', 'MYRISTYL', 'MYRISTYL', 'MYRISTYL', 'MYRISTYL', 'MYRISTYL', 'MYRISTYL', 'MYRISTYL', 'MYRISTYL', 'MYRISTYL', 'AMIDATION', 'AMIDATION', 'RGD', 'ATP_GTP_A']

dates = [75.0, 547.0, 607.0, 650.0, 22.0, 37.0, 42.0, 57.0, 71.0, 662.0, 27.5, 39.5, 59.5, 97.5, 143.5, 158.5, 214.5, 223.5, 245.5, 327.5, 391.5, 425.5, 511.5, 548.5, 567.5, 632.5, 668.5, 46.0, 77.0, 159.0, 224.0, 318.0, 521.0, 526.0, 585.0, 590.0, 598.0, 633.0, 657.0, 667.0, 209.0, 74.0, 106.0, 141.0, 252.0, 277.0, 332.0, 356.0, 486.0, 501.0, 630.0, 69.0, 228.0, 166.5, 440.0]

# Convert date strings (e.g. 2014-10-18) to datetime
# dates = [datetime.strptime(d, "%Y-%m-%d") for d in dates]


# Choose some nice levels
levels = np.tile(np.arange(-9, 9 , 2),
                 int(np.ceil(len(dates)/9)))[:len(dates)]

# Create figure and plot a stem plot with the date
fig, ax = plt.subplots(figsize=(15, 7), constrained_layout=True)
ax.set(title="ProSite domains of sseq")



### STEM PLOT

plt.axhline(y=0, color='grey', linestyle='-')

markerline, stemline, baseline = ax.stem(dates, levels,
                                         linefmt="C3-", basefmt="k-",
                                         use_line_collection=True)

plt.setp(markerline, mec="k", mfc="w", zorder=3)

# Shift the markers to the baseline by replacing the y-data by zeros.
markerline.set_ydata(np.zeros(len(dates)))

# annotate lines
vert = np.array(['top', 'bottom'])[(levels > 0).astype(int)]
for d, l, r, va in zip(dates, levels, names, vert):
    ax.annotate(r, xy=(d, l), xytext=(3, np.sign(l)*3),
                textcoords="offset points", va=va, ha="right", **{'fontsize':6, 'rotation':'vertical'}) # xytext=(-3, np.sign(l)*3), textcoords="offset points"




#### BARPLOT

# plt.bar(dates, np.full(len(dates), 1))



# format xaxis with 4 month intervals
plt.xticks(list(range(0,800, 100)))
plt.ylim(-15, 13)

# remove y axis and spines
ax.get_yaxis().set_visible(False)
# for spine in ["left", "top", "right"]:
#     ax.spines[spine].set_visible(False)

ax.margins(y=0.1)
plt.show()
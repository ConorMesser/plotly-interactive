from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import pandas as pd
import numpy as np
import plotly.graph_objs as go
import plotly.express as px
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
from ipywidgets import widgets

# params
jitter = 0.05
filter_disc_num = 10

# read in maf files todo from gcloud
maf_files = [pd.read_csv("~/Documents/Data/plotly/RP-1886_CM53575_v1_Exome_OnPrem_pair.validated.maf.txt", sep='\t'),
             pd.read_csv("~/Documents/Data/plotly/RP-1886_CM66515_v1_Exome_OnPrem_pair.validated.maf.txt", sep='\t'),
             pd.read_csv("~/Documents/Data/plotly/RP-1886_CM68786_v1_Exome_OnPrem_pair.validated.maf.txt", sep='\t')]

# sort samples by purity
# combine samples/patients into single dataframe
maf_keys = [df['Tumor_Sample_Barcode'][0] for df in maf_files]
samples = pd.concat(maf_files, keys=maf_keys)

try:
    half_purity = samples['purity'][0] / 2  # gets purity from first mutation of each sample todo****
except KeyError:
    half_purity = input("Purity call failed.\nEnter purity values for each sample (separated by commas): ")
    half_purity = half_purity.split(', ')  # [0.94, 0.15, 0.55]
    half_purity[:] = [float(x) / 2 for x in half_purity]

# sort samples by purity value to determine x_values  todo better add purity to df and then sort
helper_arr = [(pur, order) for order, pur in enumerate(half_purity)]
helper_arr = sorted(helper_arr)
helper_arr = [(order, pur, x_val + 1) for x_val, (pur, order) in enumerate(helper_arr)]
helper_arr = sorted(helper_arr)

x_vals = [[x_val] * len(maf_files[order]) for (order, _, x_val) in helper_arr]
x_vals = [item for sublist in x_vals for item in sublist]

# add random jitter to x_values
x_cat = [int(level) + (np.random.random() - 0.5) * jitter for level in x_vals]
samples['x_jitter'] = x_cat

# select variables to plot
y_axis = 'tumor_f'
gnomad_freq = 'gnomADg_AF'

# variables for hover text todo clean up NaNs
scatter_hover_1 = samples['Hugo_Symbol'].values  # todo or tolist()?
scatter_hover_2 = samples['Protein_Change'].values
scatter_hover_3 = samples['t_ref_count'].values + samples['t_alt_count'].values

# clean up fields with multiple values
temp_gnomad = samples[gnomad_freq]
for i, t in enumerate(temp_gnomad):
    if type(t) == str:
        temp_gnomad[i] = float(t.split(',')[0])
    else:
        temp_gnomad[i] = t
# samples[gnomad_freq] = [s.split(',')[0] for s in samples[gnomad_freq]]
# samples[gnomad_freq] = pd.to_numeric(samples[gnomad_freq])

fig = make_subplots(rows=2, cols=1, row_heights=[0.1, 0.9], specs=[[{"type": "heatmap"}],
                                                                   [{"type": "scatter"}]])

# add Scatter plot of sample data
filter_partitions = np.linspace(0, max(samples[gnomad_freq]), num=filter_disc_num) # todo max or set number for filter?
for step in filter_partitions:
    filtered_df = samples[samples[gnomad_freq] > step]
    fig.add_trace(go.Scatter(visible=False,
                             y=filtered_df[y_axis], x=filtered_df['x_jitter'],
                             marker_color=samples[gnomad_freq], mode='markers',
                             customdata=np.stack((scatter_hover_1,
                                                  scatter_hover_2,
                                                  scatter_hover_3), axis=-1),
                             hovertemplate='<extra></extra>' +
                                           'allelic_fraction: %{y:.3f}<br>' +
                                           'Gene: %{customdata[0]} <br>' +
                                           'Protein Change: %{customdata[1]} <br>' +
                                           'Read Depth: %{customdata[2]:d}'),
                  row=2, col=1)

# make all data points (marker_color > 0) visible initially
fig.data[0].visible = True

# create slider for filter
steps = []
for i in range(len(fig.data)):
    step = dict(
        method="update",
        args=[{'visible': [False] * (len(fig.data) + 1)}],
        label=f"{filter_partitions[i]:.2f}")
    step["args"][0]["visible"][i] = True  # toggle ith trace to visible
    step["args"][0]["visible"][-1] = True  # toggle last trace (heatmap) to visible todo
    steps.append(step)

sliders = [dict(
    active=0,
    currentvalue_visible=False,
    pad={"t": 100},
    steps=steps
)]

fig.update_layout(sliders=sliders)

# add heatmap plot of sum data
fig.add_trace(go.Heatmap(z=samples.sort_values(by='x_jitter').groupby(level=[0], sort=False).size().to_list(),
                         x=samples.sort_values(by='x_jitter').index.get_level_values(0).unique(), y=[''] * len(samples.groupby(level=0)),
                         colorscale='Viridis', showscale=False,
                         hovertemplate='<extra></extra>' +
                                       'Number of Mutations: %{z:d}'),
              row=1, col=1)

# add lines (purity / 2) to scatter plot
line_shape_dicts = []
for _, purity, x_val in helper_arr:
    x0 = x_val - 0.4
    x1 = x_val + 0.44
    line_shape_dicts.append(dict(type='line',
                                 yref='y2', y0=purity, y1=purity,
                                 xref='x2', x0=x0, x1=x1,
                                 line_color='red', line_width=2.5))
fig.update_layout(shapes=line_shape_dicts)
# todo add numbers (using add_annotation)


##################

# layout updates
fig.update_layout(title='Test Plot')

# heatmap layout
fig.update_xaxes(showticklabels=False, ticks="outside", tick0=1, dtick=1, tickwidth=2, row=1)
fig.update_yaxes(visible=False, row=1)

# scatter updates
fig.update_xaxes(showgrid=False, range=[0.5, 3.5], ticks="outside", tickwidth=2, tickmode='array', ticktext=['Sample1', 'Sample2', 'Sample3'], tickvals=[1, 2, 3], tickangle=270, row=2)
fig.update_yaxes(showgrid=False, range=[-0.1, 1.1], ticks="outside", tickwidth=2, tick0=0, dtick=0.2, zeroline=False, row=2)
fig.update_layout(
    hoverlabel=dict(
        bgcolor="white",
        font_size=16,
        font_family="Rockwell"
    )
)


# use selectedpoints in scatter plot to specify which mutations are silent
#    (can be made opaque, or very small using button update_layout)  ???

# Use jupyterNotebook to add more interactivity

# figure out how to use dash??


# look into widgets**

# plot fig
fig.show()


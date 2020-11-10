import pandas as pd
import numpy as np
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import csv
import dalmatian


def get_static_plot():
    # params
    jitter = 0.3
    line_length = 0.8
    workspace = 'broad-firecloud-ibmwatson/Getz_Ebert_IBM_13-583_Exomes_Analysis'

    # get list of maf files
    maf_list_filename = input("Filename for list of maf files: ")  # can be either local files or from google cloud urls
    maf_files = []
    with open(maf_list_filename) as file:
        filename_reader = csv.reader(file)
        for row in filename_reader:
            this_maf = pd.read_csv(row[0], sep='\t')
            maf_files.append(this_maf)

    # combine samples/patients into single dataframe
    maf_keys = [df['Tumor_Sample_Barcode'][0] for df in maf_files]
    samples = pd.concat(maf_files, keys=maf_keys)

    try:
        # get purity values from Terra
        wm = dalmatian.WorkspaceManager(workspace)
        pairs = wm.get_pairs()
        purities = pairs.loc[[key + '_pair' for key in maf_keys], 'wxs_purity']
        half_purity = purities.apply(lambda x: float(x) / 2)
        half_purity.index = [maf_key[:-5] for maf_key in half_purity.index]
    except Exception as e:
        print(e)
        purities = input("\nPurity call failed.\nEnter purity values for each sample (separated by commas): ")
        purities = purities.split(', ')
        half_purity = pd.Series(maf_keys, [float(x) / 2 for x in purities])

    # sort samples df by purity
    full_pur_list = [half_purity[s[0]] for s in samples.index]
    samples['half_purity'] = full_pur_list
    samples.sort_values('half_purity', inplace=True)  # doesn't matter the order within each sample
    half_purity.sort_values(inplace=True)

    # get x_values for plotting scatter
    sum_data = samples.groupby(level=[0], sort=False).size().to_list()
    x_vals = [[x_val] * sum_data[x_val-1] for x_val in range(1, len(sum_data) + 1)]
    x_vals = [item for sublist in x_vals for item in sublist]

    # add random jitter to x_values
    x_cat = [int(level) + (np.random.random() - 0.5) * jitter for level in x_vals]
    samples['x_jitter'] = x_cat

    # select variables to plot
    gnomad_freq = 'gnomADg_AF'

    # clean up multiple frequencies (take first value)
    samples[gnomad_freq] = samples[gnomad_freq].apply(lambda x: float(x.split(',')[0] if type(x) == str else x))

    fig = make_subplots(rows=2, cols=1,
                        row_heights=[0.1, 0.9],
                        shared_xaxes=True,
                        specs=[[{"type": "heatmap"}],
                               [{"type": "scatter"}]])

    # add Scatter plot of sample data
    gnomad_min = samples[gnomad_freq].min()
    gnomad_max = samples[gnomad_freq].max()
    fig.add_trace(go.Scatter(visible=True,
                             y=samples['tumor_f'], x=samples['x_jitter'],
                             marker_color=samples[gnomad_freq], mode='markers',
                             marker_reversescale=True,
                             # marker_showscale=True,
                             marker_cmin=gnomad_min,
                             marker_cmax=gnomad_max,
                             customdata=np.stack((samples['Hugo_Symbol'].tolist(),
                                                  samples['Protein_Change'].tolist(),
                                                  samples['t_ref_count'].values + samples['t_alt_count'].values),
                                                 axis=-1),
                             hovertemplate='<extra></extra>' +
                                           'allelic_fraction: %{y:.3f}<br>' +
                                           'Gene: %{customdata[0]} <br>' +
                                           'Protein Change: %{customdata[1]} <br>' +
                                           'Read Depth: %{customdata[2]:d}'),
                  row=2, col=1)

    # add heatmap plot of sum data
    fig.add_trace(go.Heatmap(z=sum_data,
                             x=samples.sort_values(by='x_jitter').index.get_level_values(0).unique(),
                             y=[''] * len(samples.groupby(level=0)),
                             colorscale='Viridis', showscale=False,
                             zmin=0, zmax=max(sum_data),
                             hovertemplate='<extra></extra>' +
                                           'Num Mutations: %{z:d}'),
                  row=1, col=1)

    # add annotation for number of mutations on heat map
    for x_val, sum in enumerate(sum_data):
        fig.add_annotation(x=x_val, y=0, text=sum, showarrow=False, font_color="#ffffff")

    # add lines (purity / 2) to scatter plot
    line_shape_dicts = []
    for x_val, purity in enumerate(half_purity.tolist()):

        x0 = x_val + (1 - line_length / 2)  # x from enumerate is zero-indexed
        x1 = x_val + (1 + line_length / 2)
        line_shape_dicts.append(dict(type='line',
                                     yref='y2', y0=purity, y1=purity,
                                     xref='x2', x0=x0, x1=x1,
                                     line_color='red', line_width=2.5))
    fig.update_layout(shapes=line_shape_dicts)

    # layout updates
    fig.update_layout(title='Allelic Frequency Plot')

    # heatmap layout
    fig.update_xaxes(showticklabels=False, ticks="outside",
                     tick0=1, dtick=1, tickwidth=2, row=1)
    fig.update_yaxes(visible=False, row=1)

    # scatter updates
    sample_names = [name[8:15] for name in half_purity.index.tolist()]
    tick_x = np.linspace(1, len(sample_names), len(sample_names))
    fig.update_xaxes(showgrid=False,
                     range=[0.5, len(sample_names) + 0.5],
                     ticks="outside", tickwidth=2, tickmode='array',
                     ticktext=sample_names, tickvals=tick_x,
                     tickangle=270, row=2)
    fig.update_yaxes(showgrid=False,
                     range=[-0.1, 1.1],
                     ticks="outside", tickwidth=2,
                     tick0=0, dtick=0.2, zeroline=False,
                     row=2)
    fig.update_layout(
        hoverlabel=dict(
            bgcolor="white",
            font_size=16,
            font_family="Rockwell"
        )
    )

    return fig, samples, gnomad_min, gnomad_max, gnomad_freq


def get_filenames(workspace, output_file='maf_filenames.txt'):
    # Import Workspace from Firecloud
    wm = dalmatian.WorkspaceManager(workspace)
    pairs = wm.get_pairs()
    # pairs = wm.get_pairs_in_pair_set('DESIRED_PAIR_SET')

    pairs = pairs[pairs['mutation_validator_validated_maf'].notnull()]

    ###################
    # Perform Filtering

    ###################

    desired_files = pairs['mutation_validator_validated_maf'].tolist()

    with open(output_file, 'w') as o_file:
        o_file.writelines(f'{f_name}\n' for f_name in desired_files)

# plotly-interactive

Allows interactive plotting using plotly backend. Implements a plot showing allelic frequencies of mutations for many samples, with filtering for mutation type and GnomAD frequency. 

Currently, interactivity lives in a Jupyter Notebook but it can be extended to work with Plotly Dash.

# Usage

Modify get_filenames() to put together a txt file of maf filenames (can be local files or on gcloud), each on their own line.

Open the Jupyter Notebook and run the cells, inputting the created maf_file txt file. The plot will be output after the last cell.

import dalmatian
import numpy as np
import csv
import os
import argparse


def collect_filenames(filepath=None, directory=None,
                      workspace=None, set_type=None,
                      set_name=None, column_name=None):
    """Get a list of filenames/urls from a text file, directory, or Terra workspace.

    Filenames can be retrieved from a specified filepath with each filename (/path) delimited by newline, comma, or tab.
    A specified directory will return the paths for every file inside it. A specified workspace will return the fields
    in the given column_name of a particular set (set_name of set_type: 'pair', 'sample', 'participant'). If multiple
    options are given, priority is filepath, directory, then workspace.

    :param filepath: path and name of file to parse for desired filenames
    :param directory: path and name of directory with desired files
    :param workspace: full name of Terra workspace
    :param set_type: one of [pair, sample, participant]
    :param set_name: name of desired set
    :param column_name: column name holding desired file/data
    :return: list of filenames
    """
    if filepath:
        with open(filepath) as f:
            reader = csv.reader(f, delimiter="\t")  # split by tab
            flat_list = [x for sublist in list(reader) for x in sublist]  # split by row and gives single-nested list
            lines = [row.split(',') for row in flat_list]  # split by comma
            filename_list = [name.strip() for sublist in lines for name in sublist]  # strip whitespace + give flat list
    elif directory:
        filename_list = []
        for path in os.listdir(directory):
            full_path = os.path.join(directory, path)
            if os.path.isfile(full_path):  # check to ensure this is a file (not a directory)
                filename_list.append(full_path)
    elif workspace:
        final_set = df_from_workspace_set(workspace, set_type, set_name, column_name)
        filename_list = final_set[column_name].tolist()
    else:
        raise ValueError("One of filepath, directory, or workspace options must be specified.")

    return filename_list


def df_from_workspace_set(workspace, set_type, set_name, column_name):
    """Get the dataframe specified by set_name (of set_type) in your Terra workspace.

    The dataframe is filtered by indices that have non-null values in the column specified by column_name. The whole
    dataframe is returned so that further filtering can be performed beyond what is in the given set_name.

    :param workspace: full name of Terra workspace
    :param set_type: one of [pair, sample, participant]
    :param set_name: name of desired set
    :param column_name: column name holding desired file/data
    :return: pandas.Dataframe of the desired set
    """
    if not set_type or not set_name or not column_name:
        raise ValueError("If calling from Terra workspace, the set_type, "
                         "set_name, and column_name must be specified.")

    # Import Workspace from Terra/Firecloud
    wm = dalmatian.WorkspaceManager(workspace)
    if set_type == 'pair':
        set_df = wm.get_pairs_in_pair_set(set_name)
    elif set_type == 'sample':
        set_df = wm.get_samples()
        set_df = set_df[np.in1d(set_df.index.values, wm.get_sample_sets().loc[set_name]['samples'])]
    elif set_type == 'participant':
        set_df = wm.get_participants()
        set_df = set_df[np.in1d(set_df.index.values, wm.get_participant_sets().loc[set_name]['participants'])]
    else:
        raise ValueError(f"set_type must be one of pair, sample, participant, not {set_type}.")

    return set_df[set_df[column_name].notnull()]


def main():
    parser = argparse.ArgumentParser(description='get a list of filenames from a file, directory, or Terra workspace')

    parser.add_argument("-f", "--filepath", help='file with desired filenames separated by linebreaks', nargs='*')
    parser.add_argument("-d", "--directory", help='directory containing all desired files', nargs='*')
    parser.add_argument("-w", "--workspace", help='Terra/firecloud workspace name. If given in addition to -f/-d, '
                                                  'workspace name will be output from method', nargs='*')
    parser.add_argument("-st", "--set_type", choices=['pair', 'sample', 'participant'])
    parser.add_argument("-sn", "--set_name")
    parser.add_argument("-cn", "--column_name")

    parser.add_argument('--output_path', help='desired path for output file')

    args = parser.parse_args()

    # saves the parser if there is a space (in a path name or workspace, God forbid)
    #   just a convenience tool, so spaces don't have to be \ (back-slashed)
    if args.filepath:
        args.filepath = ' '.join(args.filepath)
    if args.directory:
        args.directory = ' '.join(args.directory)
    if args.workspace:
        args.workspace = ' '.join(args.workspace)

    filenames_list = collect_filenames(filepath=args.filepath, directory=args.directory,
                                       workspace=args.workspace, set_type=args.set_type,
                                       set_name=args.set_name, column_name=args.column_name)

    print("Filenames: ")
    print(*filenames_list, sep='\n')


if __name__ == '__main__':
    main()

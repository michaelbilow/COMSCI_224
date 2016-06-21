import shutil
from os.path import join, split, splitext, exists
from os import listdir, makedirs, walk
from distutils.dir_util import copy_tree


if __name__ == "__main__":
    input_folder = '../Final Projects'
    data_folder = '../Final_Project_Data'
    solutions_folder = '../Final_Project_Solutions'
    if not exists(data_folder):
        makedirs(data_folder)
    if not exists(solutions_folder):
        makedirs(solutions_folder)

    for root, dirs, files in walk(input_folder):
        if all(_ in dirs for _ in ('input', 'output')):
            data_output_folder = join(data_folder, '/'.join(root.split('/')[2:]))
            solutions_output_folder = join(solutions_folder, '/'.join(root.split('/')[2:]))
            copy_tree(join(root, 'input'), join(data_output_folder, 'input'))

            this_output_folder = solutions_output_folder if split(root)[1] == 'test' else data_output_folder
            copy_tree(join(root, 'output'), join(this_output_folder, 'output'))
            # if split(root) == 'test':

            print root
            print dirs
            print files
            print data_output_folder, solutions_output_folder
            print '\n\n'
    # for final_prj_folder in [join(input_folder, _) for _ in listdir(input_folder)]:
    #     difficulties = [join(final_prj_folder, _) for _ in listdir(final_prj_folder)]
    #     for diff in difficulties:
    #         datasets = [join(diff, _) for _ in listdir(diff)]
    #         for dataset in datasets:
    #             in_outs = [join(dataset, _) for _ in listdir(dataset)]
    #             for in_out in in_out:
    #                 if dataset == 'test' and in_out = 'output':
    #                     shutil.copy()
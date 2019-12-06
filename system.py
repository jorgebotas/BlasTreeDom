#!/usr/bin/env python

from subprocess import call, PIPE, Popen

def list_files(dir):
    # Create iterable list with all the files in input directory
    dir_ls = Popen(['ls', os.path.abspath(dir)], stdout = PIPE, stderr = PIPE)
    dir_list = dir_ls.stdout.read().decode('utf-8').splitlines()
    dir_ls.stdout.close()
    dir_ls.stderr.close()
    return dir_list

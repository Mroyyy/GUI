import os

def get_path():
    filepath = '/home/gallegolab/Desktop/Maria/Saccharomyces/Pymol/SSO1'
    fpath = os.path.realpath(filepath)

    return fpath

print(get_path)
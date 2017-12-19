import os


def save_mkdir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

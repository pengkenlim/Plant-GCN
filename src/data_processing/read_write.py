#setting sys.path for importing modules
import os
import sys

if __name__ == "__main__":
         abspath= __file__
         parent_module= "/".join(abspath.split("/")[:-2])
         print(parent_module)
         sys.path.insert(0, parent_module)

import pickle

# func definitions
def load_pickle(read_path):
    with open(read_path, 'rb') as f:
        data = pickle.load(f)
    return data

def to_pickle( data , write_path):
    with open(write_path, 'wb') as f:
        pickle.dump(data, f)


import sys

import numpy as np
import matplotlib.pyplot as plt


def prepare_file():
    
    try:
        y_file_name = sys.argv[1]
    except:
        print("Usage: " + sys.argv[0] + " data_Transport_file_name"); sys.exit(1)

    y_spt  = y_file_name.split("_")
    x_file_name = ''
    
    print(y_spt)
    
    for y in y_spt[:5]:
        x_file_name += y + '_'

    print("First loop: ",x_file_name)

    for y in y_spt[7:]:
        x_file_name += y + '_'
    
    print("Second loop: ",x_file_name)

    x_file_name = x_file_name[:-1]

    x_data = np.load(x_file_name)
    y_data = np.load(y_file_name)

    # x_data = 0
    # y_data = 0

    return x_data, y_data



def main():
    x, y = prepare_file()

    plt.plot(x, y)
    plt.show()


if __name__ == "__main__":
    main()
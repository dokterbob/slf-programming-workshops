import numpy as np


def process_data(data):
    """
    Process data. Expects a list of lists with numbers.

    * Convert data to 2D numpy array
    * Rotate the 2D array
    """
    # Convert to Numpy float32
    data = np.array(data, np.float32)

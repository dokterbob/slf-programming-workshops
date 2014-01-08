from input import import_data
from process import process_data


# Read the data
data = import_data('data.txt')
data = process_data(data)

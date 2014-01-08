def read_file(filename):
    # Open the file
    file_obj = open(filename)

    # Iterate over lines in the file
    for line in file_obj:

        # Split line by spaces (creates a list)
        # Alternatives: split(',')
        numbers = line.split()

        if len(numbers) != 2:
            # Convert strings to numbers

            # map() calls the first argument on every item in the second
            # argument and returns a list of results.
            numbers = map(float, numbers)
        else:
            # We're processing a header
            print 'Skipping header line'

    return contents

# Just for debugging purposes
read_file('data.txt')

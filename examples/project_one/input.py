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

            numbers2 = []
            for number in numbers:
                # Convert number to float
                number = float(number)

                # Append to temperary list
                numbers2.append(number)

            # Replace numbers by numbers2
            numbers = numbers2

        else:
            # We're processing a header
            print 'Skipping header line'

    return contents

# Just for debugging purposes
read_file('data.txt')

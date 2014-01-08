def import_data(filename):
    """
    Read filename, skip headers, parse numbers.
    Returns list with lists.
    """

    # Open the file
    file_obj = open(filename)

    contents = []
    imported_lines = 0

    # Iterate over lines in the file
    for line in file_obj:

        # Split line by spaces (creates a list)
        # Alternatives: split(',')
        numbers = line.split()

        if len(numbers) != 2:
            # Convert strings to numbers

            # Assert that the column number corresponds to our
            # expectations.
            assert len(numbers) == 103, 'Unkown column number.'

            # map() calls the first argument on every item in the second
            # argument and returns a list of results.
            numbers = map(float, numbers)

            contents.append(numbers)
            imported_lines += 1
        else:
            # We're processing a header
            print 'Skipping header line'

    assert imported_lines == len(contents)
    print 'Imported %d lines' % imported_lines

    return contents

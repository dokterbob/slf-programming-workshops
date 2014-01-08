def read_file(filename):
    file_obj = open(filename)

    contents = file_obj.read()

    return contents

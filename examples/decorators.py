def my_function():
    # <DO SOMETHING>
    pass

my_function = other_function(my_function)


@other_function
def my_function():
    """ Equivalent to previous example. """
    pass

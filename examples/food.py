class Food(object):
    """ Something you can eat. """

    # Class properties are shared across instances
    destination = 'mouth'

    def __init__(self, name):
        """
        The initialization method is called whenever
        instances are created.
        """
        # Store the name property on the object
        self.name = name

        # Create the object unused
        self.used = False

    def use(self):
        """ Use this particular item, marking it as used. """

        print 'Using %s in %s' % (self.name, self.destination)

        # Instance properties are specific to one instance
        self.used = True


lunch = Food(name='banana')
lunch.used
lunch.use()
lunch.used

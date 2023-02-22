class VCException(Exception):
    """

    """
    min_queue_size = 1
    max_queue_size = 10

    def __init__(self, queue_size, *args):
        """

        @param queue_size:
        @param args:
        """
        super().__init__(args)
        self.queue_size = queue_size

    def __str__(self):
        """

        @return:
        """
        return f'The provided queue size ({self.queue_size}) is not in a valid range {self.min_queue_size, self.max_queue_size}'

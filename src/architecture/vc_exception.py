import src.config_util.config_io as cio


class VCException(Exception):
    """
    Exception for VCQueue
    """

    def __init__(self, queue_size: int, *args):
        """
        Constructor
        @param queue_size: queue size so exception can be checked and handled accordingly
        """
        super().__init__(args)
        self.queue_size = queue_size
        self.min_queue_size = cio.get_min_queue_size()
        self.max_queue_size = cio.get_max_queue_size()

    def __str__(self) -> str:
        """
        Function to stringify exception
        @return: str representation of exception
        """
        return (
            f"The provided queue size ({self.queue_size}) is not in a valid range "
            f"{self.min_queue_size, self.max_queue_size} "
        )

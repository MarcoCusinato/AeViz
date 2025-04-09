from . import wraps
import matplotlib.pyplot as plt

def fig_window_open(func):
    """
    Check if the window figure has been closed. If it has it destroys it.
    """
    @wraps(func)
    def check_opens(*args, **kwargs):
        if args[0].fig is not None:
            if not plt.fignum_exists(args[0].fig.number):
                args[0].Close()
        return func(*args, **kwargs)
    return check_opens

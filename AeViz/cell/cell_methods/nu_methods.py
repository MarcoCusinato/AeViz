def E_nu_left(self):
    """
    Method that returns an array with the left neutrino energy grid.
    Results:
        E_nu_left: (numpy array)
    """
    if self.path_grid is None:
        raise TypeError("You shouldn't be here.")
    return self._cell__nu_grid_file[:, 0]

def E_nu_right(self):
    """
    Method that returns an array with the right neutrino energy grid.
    Results:
        E_nu_right: (numpy array)
    """
    if self.path_grid is None:
        raise TypeError("You shouldn't be here.")
    return self._cell__nu_grid_file[:, 2]

def E_nu(self):
    """
    Method that returns an array with the neutrino energy grid.
    Results:
        E_nu: (numpy array)
    """
    if self.path_grid is None:
        raise TypeError("You shouldn't be here.")
    return self._cell__nu_grid_file[:, 1]

def dE_nu(self):
    """
    Method that returns an array with the neutrino energy integration element.
    Results:
        dE_nu: (numpy array)
    """
    return self.E_nu_right() - self.E_nu_left()
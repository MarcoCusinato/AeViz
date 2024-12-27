import numpy as np

## X
def x_left(self, ghost):
    """
    Method that returns an array with the 'left' coordinates of the
    x, with the selected number of ghost cells
    Parameters:
        ghost: (object) ghost
    Results:
        left x coordinates: (numpy array)
    """
    return ghost.remove_ghost_cells(self._cell__radius_file[:, 0], self.dim,
                                    'radius')

def x_right(self, ghost):
    """
    Method that returns an array with the 'right' coordinates of the 
    x, with the selected number of ghost cells
    Parameters:
        ghost: (object) ghost
    Results:
        right x coordinates: (numpy array)
    """
    return ghost.remove_ghost_cells(self._cell__radius_file[:, 2], self.dim,
                                    'radius')

def x(self, ghost):
    """
    Method that returns an array with the 'center' coordinates of the
    x, with the selected number of ghost cells
    Parameters:
        ghost: (object) ghost
    Results:
        center x coordinates: (numpy array)
    """
    return ghost.remove_ghost_cells(self._cell__radius_file[:, 1], self.dim,
                                    'radius')

def dx(self, ghost):
    """
    Method that returns an array with the lenght of each radial cell
    Parameters:
        ghost: (object) ghost
    Results:
        radial lenght of a cell: (numpy array)
    """
    return self.x_right(ghost) - self.x_left(ghost)

## Y
def y_left(self, ghost):
    """
    Method that returns an array with the 'left' coordinates of the
    y, with the selected number of ghost cells
    Parameters:
        ghost: (object) ghost
    Results:
        left y coordinates: (numpy array)
    """
    try:
        return ghost.remove_ghost_cells(self._cell__theta_file[:, 0], self.dim,
                                    'theta')
    except:
        return self._cell__theta_file[0]
    
def y_right(self, ghost):
    """
    Method that returns an array with the 'right' coordinates of the 
    y, with the selected number of ghost cells
    Parameters:
        ghost: (object) ghost
    Results:
        right y coordinates: (numpy array)
    """
    try:
        return ghost.remove_ghost_cells(self._cell__theta_file[:, 2], self.dim,
                                    'theta')
    except:
        return self._cell__theta_file[2]

def y(self, ghost):
    """
    Method that returns an array with the 'center' coordinates of the
    y, with the selected number of ghost cells
    Parameters:
        ghost: (object) ghost
    Results:
        center y coordinates: (numpy array)
    """
    try:
        return ghost.remove_ghost_cells(self._cell__theta_file[:, 1], self.dim,
                                    'theta')
    except:
        return self._cell__theta_file[1]

def dy(self, ghost):
    """
    Method that returns an array with the lenght of each radial cell
    Parameters:
        ghost: (object) ghost
    Results:
        radial lenght of a cell: (numpy array)
    """
    if self.dim < 2:
        return 1
    return self.y_right(ghost) - self.y_left(ghost)

## Z
def z_left(self, ghost):
    """
    Method that returns an array with the 'left' coordinates of the
    z, with the selected number of ghost cells
    Parameters:
        ghost: (object) ghost
    Results:
        left z coordinates: (numpy array)
    """
    try:
        return ghost.remove_ghost_cells(self._cell__phi_file[:, 0], self.dim,
                                    'phi')
    except:
        return self._cell__phi_file[0]

def z_right(self, ghost):
    """
    Method that returns an array with the 'right' coordinates of the 
    z, with the selected number of ghost cells
    Parameters:
        ghost: (object) ghost
    Results:
        right z coordinates: (numpy array)
    """
    try:
        return ghost.remove_ghost_cells(self._cell__phi_file[:, 2], self.dim,
                                    'phi')
    except:
        return self._cell__phi_file[2]
    
def z(self, ghost):
    """
    Method that returns an array with the 'center' coordinates of the
    z, with the selected number of ghost cells
    Parameters:
        ghost: (object) ghost
    Results:
        center z coordinates: (numpy array)
    """
    try:
        return ghost.remove_ghost_cells(self._cell__phi_file[:, 1], self.dim,
                                    'phi')
    except:
        return self._cell__phi_file[1]
    
def dz(self, ghost):
    """
    Method that returns an array with the lenght of each radial cell
    Parameters:
        ghost: (object) ghost
    Results:
        radial lenght of a cell: (numpy array)
    """
    if self.dim < 3:
        return 1
    return self.z_right(ghost) - self.z_left(ghost)

def dVolume(self, ghost):
    """
    Method that returns an array with the volume of each cell
    Parameters:
        ghost: (object) ghost
    Results:
        volume of a cell: (numpy array)
    """
    if self.dim == 1:
        return self.dx(ghost)
    elif self.dim == 2:
        return self.dx(ghost)[None, :] * self.dy(ghost)[:, None] * \
            self.dz(ghost)
    else:
        return self.dx(ghost)[None, None, :] * \
            self.dy(ghost)[None, :, None] * self.dz(ghost)[:, None, None]
    

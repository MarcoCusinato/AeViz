import numpy as np
from scipy.special import factorial, binom
import numpy as np

class AssociatedLegendrePolynomials:
    def __init__(self):
        pass
    
    def change_m_sign(self, m, l, P):
        """
        Renormaqlizes the polynomial order from m to -m, using the
        relation:
        P^m_l(x) = (-1)^m (l-m)!/(l+m)! P^-m_l(x)
        Takes as input:
            l: original degree of the polynomial
            m: original order of the polynomial
            P: original polynomial
        """
        assert m >= 0, "m must be a positive integer"
        assert l >= m, "l must be greater or equal than m"
        return (-1)**m * factorial(l-m) / factorial(l+m) * P
    
    def P(self, m, l, x):
        """
        Calculates the associated Legendre polynomial of order m and
        degree l. Uses the polynomials closed form expression:
        P^m_l(x) = (-)^m 2^l (1-x^2)^(m/2) Σ_k=m^l k!/(k-m)! binom(
            (l+k-1)/2, l)
        Takes as input:
            m: order of the polynomial
            l: degree of the polynomial
            x: value at which the polynomial is evaluated, array or
                scalar
        Returns:
            P^m_l(x): value of the polynomial at x, array
        """
        assert np.abs(m) <= l, "m must be smaller or equal than l"
        x = np.asarray(x)
        if l == 0:
            return 1
        if m<0:
            return self.change_m_sign(-m, l, self.P(-m, l, x))
        else:
            summ = np.zeros(x.shape)
            for k in range(m, l+1):
                summ += factorial(k) / factorial(k-m) * \
                    binom((l + k - 1) / 2, l) * x ** (k - m)
            return (-1)**m * 2**l * (1 - x ** 2) ** (m / 2) * summ
    

class SphericalHarmonics(AssociatedLegendrePolynomials):
    def __init__(self):
        pass
    
    def Ylm(self, m, l, theta, phi):
        """
        Calculates the spherical harmonic of order m and degree l.
        Uses the definition:
        Y^m_l(theta, phi) = (-)^m (2l+1)/(4π) (l-m)!/(l+m)! P^m_l(
            cos(theta)) exp(i m phi)
        Takes as input:
            m: order of the polynomial
            l: degree of the polynomial
            theta: azimutal angle, array or scalar
            phi: polar angle, array or scalar
        Returns:
            Y^m_l(theta, phi): value of the polynomial at (theta, phi),
                array
        """
        assert np.abs(m) <= l, "m must be smaller or equal than l"
        theta = np.cos(np.asarray(theta))
        phi = np.asarray(phi)
        if phi.size > 1:
            phi = phi[..., None]
            theta = theta[None, ...]
        norm = np.sqrt((2 * l + 1) / (4 * np.pi) * factorial(l - m) / \
            factorial(l + m))
        Plm = self.P(m, l, theta)
        return norm * Plm * np.exp(1j * m * phi)
    
    def Ylm_conj(self, m, l, theta, phi):
        """
        Calculates the complex conjugate of the spherical harmonic of
        order m and degree l.
        Takes as input:
            m: order of the polynomial
            l: degree of the polynomial
            theta: azimutal angle, array or scalar
            phi: polar angle, array or scalar
        Returns:
            Y^m_l(theta, phi): value of the polynomial at (theta, phi),
                array
        """
        return (-1) ** m * self.Ylm(-m, l, theta, phi)
    
    def Ylm_norm(self, m, l, theta, phi):
        if m < 0:
            return np.sqrt(2) * (-1) ** m * self.Ylm(m, l, theta, phi).imag
        elif m == 0:
            return self.Ylm(m, l, theta, phi).real
        else:
            return np.sqrt(2) * (-1) ** m * self.Ylm(m, l, theta, phi).real

    def spin_weighted_Ylm(self, s, m, l, theta, phi):
        """
        Calculates the spin weighted spherical harmonic of order m and
        degree l, with spin s.

        """
        assert np.abs(m) <= l, "m must be smaller or equal than l"
        if l < np.abs(s):
            return 0
        theta = np.asarray(theta)
        phi = np.asarray(phi)
        if phi.size > 1:
            phi = phi[..., None]
            theta = theta[None, ...]
        norm = (-1) ** (l + m - s) * \
            (((2 * l +1) / (4 * np.pi)) * (factorial(l + m) / 
                                             factorial(l + s)) * \
                (factorial(l - m) / factorial(l - s))) ** 0.5 * \
                (np.sin(0.5 * theta)) ** (2 * l) * np.exp(1j * m * phi)
        summ = 0
        for n in range(0, l - s + 1):
            summ += np.nan_to_num((-1) ** (n) * binom(l - s, n) * \
                binom(l + s, n + s - m) * (np.cos(0.5 * theta) / \
                                     np.sin(0.5 * theta)) ** (2 * n + s - m), 
                                       False,  0)
        return norm * summ
        
import unittest
import numpy as np
import matplotlib.pyplot as plt

class Orbit(object):

    """
    this class evaluates orbit parameters according to Algorithm 4.1 and 4.2 Curtis
    mu - gravitational parameter (kmˆ3/sˆ2); mu=G(m1+m2)
    r - position vector
    v - velocity vector
    r_mod, v_mod - the magnitudes of r and v
    vr - radial velocity component (km/s)
    h - the angular momentum vector (kmˆ2/s)
    h_mod - the magnitude of H (kmˆ2/s)
    incl - inclination of the orbit (rad)
    N - the node line vector (kmˆ2/s)
    N_mod - the magnitude of N
    RA - right ascension of the ascending node (grad)
    e_vec - eccentricity vector
    e - eccentricity (magnitude of e_vec)
    w - argument of perigee (grad)
    TA - true anomaly (grad)
    T - period, only for an ellipse (e<1) otherwise is 0
    pi - 3.1415926...
    """

    def __init__(self, coords, type, mu=398600):

        """
        :param coords: array of 6 elements, cartesian coordinates must be in km and km/s
            cartesian coordinates must be given as [x,y,z,vx,vy,vz]
            keplerian coordinates must be given as [h_mod, incl, RA, e, w, TA]
        :param type: string, can be 'cartesian' or 'keplerian'
        :param mu: gravitational parameter (kmˆ3/sˆ2); mu=G(m1+m2); default value is mu=398600
        """
        self.mu = mu
        self.type = type
        if type == 'cartesian':
            self.cart = {'x':coords[0],
                         'y':coords[1],
                         'z':coords[2],
                         'vx':coords[3],
                         'vy':coords[4],
                         'vz':coords[5]}
            # convert coordinates
            self.toKep()
        elif type == 'keplerian':
            self.kep = {'h':coords[0],
                         'incl':coords[1],
                         'RA':coords[2],
                         'e':coords[3],
                         'w':coords[4],
                         'TA':coords[5]}
            # convert coordinates
            self.toCart()
        else:
            raise ValueError('type must be "cartesian" or "keplerian"')

    def getKep(self):
        """
        :return: keplerian coordinates
        """
        return np.array(list(self.kep.values()))

    def getKepDict(self):
        """
        :return: keplerian coordinates dictionary
        """
        return self.kep

    def getCartDict(self):
        """
        :return: cartesian coordinates dictionary
        """
        return self.cart


    def getCart(self):
        """
        :return: cartesian coordinates
        """
        return np.array(list(self.cart.values()))

    def getSemiMajorAxes(self):

        """
        :return:  semimajor axis of ellipse in km
        """

        if self.kep['e'] >= 1:
            raise ValueError("The orbit is not an ellipse.")
        else:
            h_mod = self.kep['h']
            incl = self.kep['incl'] * np.pi / 180
            RA = self.kep['RA'] * np.pi / 180
            e = self.kep['e']
            w = self.kep['w'] * np.pi / 180
            TA = self.kep['TA'] * np.pi / 180

            rp = (h_mod**2/self.mu)/(1+e)
            ra = (h_mod**2/self.mu)/(1-e)

            return (rp+ra)/2

    def getRadius(self):

        """
        :return: orbit radius in km
        """

        r = [self.cart['x'], self.cart['y'], self.cart['z']]
        return np.linalg.norm(r)

    def getPeriod(self):

        """
        :return: orbit period in seconds
        """

        if self.kep['e'] >=1:
            raise ValueError("The orbit is not an ellipse.")
        else:
            return 2*np.pi*np.sqrt(self.getSemiMajorAxes()**3)/np.sqrt(self.mu)


    def getEnergy(self):
        pass

    def toKep(self):

        """
        Convert cartesian coordinates to kelerian
        :return: None
        """

        if self.cart is not None:
            r = np.array([self.cart['x'], self.cart['y'], self.cart['z']])
            v = np.array([self.cart['vx'], self.cart['vy'], self.cart['vz']])
        else:
            raise ValueError("cartesian coordinates are not defined")

        # 1. calculate the distance
        r_mod = np.linalg.norm(r)

        # 2. calculate the speed
        v_mod = np.linalg.norm(v)

        # 3. calculate the radial velocity
        vr = np.dot(r, v)/r_mod

        # 4. calculate the specific angular momentum
        h = np.cross(r, v)

        # 5. calculate the magnitude of h
        h_mod = np.linalg.norm(h)

        # 6. calculate the inclination
        incl = np.arccos(h[2]/h_mod)

        # 7. calculate N
        N = np.cross([0, 0, 1], h)

        # 8. calculate the magnitude of N
        N_mod = np.linalg.norm(N)

        # 9. Calculate RA of ascending node
        if N[1]>=0:
            RA = np.arccos(N[0]/N_mod)
        else:
            RA = 2*np.pi - np.arccos(N[0]/N_mod)

        # 10. calculate the eccentricity vector
        e_vec = ((v_mod**2 - self.mu/r_mod)*r - r_mod*vr*v)/self.mu

        # 11. calculate the eccentricity
        e = np.linalg.norm(e_vec)

        # 12. calculate the argument of perigee
        if e_vec[2] >= 0:
            w = np.arccos(np.dot(N, e_vec)/(N_mod*e))
        else:
            w = 2*np.pi - np.arccos(np.dot(N, e_vec)/(N_mod*e))

        # 13. calculate the true anomaly
        if vr >= 0:
            TA = np.arccos(np.dot(e_vec, r)/(e*r_mod))
        else:
            TA = 2*np.pi - np.arccos(np.dot(e_vec, r)/(e*r_mod))

        # define Keplerian coordinates
        self.kep = {'h':h_mod,
                    'incl':incl*180/np.pi,
                    'RA':RA*180/np.pi,
                    'e':e,
                    'w':w*180/np.pi,
                    'TA':TA*180/np.pi}

    def toCart(self):

        """
        Convert keplerian coordinates to cartesian
        :return: None
        """

        if self.kep is not None:
            h_mod = self.kep['h']
            incl = self.kep['incl']*np.pi/180
            RA = self.kep['RA']*np.pi/180
            e = self.kep['e']
            w = self.kep['w']*np.pi/180
            TA = self.kep['TA']*np.pi/180
        else:
            raise ValueError("keplerian coordinates are not defined")

        # 1. calculate r perifocal
        r_perif = (h_mod**2/self.mu)/(1 + e*np.cos(TA))*np.array([np.cos(TA), np.sin(TA), 0.0])


        # 2. calculate v perifocal
        v_perif = (self.mu/h_mod)*np.array([-np.sin(TA), e + np.cos(TA), 0.0])

        # 3. calulate QXx
        Q = np.array([[np.cos(RA)*np.cos(w) - np.sin(RA)*np.sin(w)*np.cos(incl),    -np.cos(RA)*np.sin(w) - np.sin(RA)*np.cos(incl)*np.cos(w),    np.sin(RA)*np.sin(incl)],
                       [np.sin(RA)*np.cos(w) + np.cos(RA)*np.sin(w)*np.cos(incl),    -np.sin(RA)*np.sin(w) + np.cos(RA)*np.cos(incl)*np.cos(w),    -np.cos(RA)*np.sin(incl)],
                       [np.sin(incl)*np.sin(w),                                       np.sin(incl)*np.cos(w),                                       np.cos(incl)]])

        # 4 calculate r and v
        r = Q.dot(r_perif)
        v = Q.dot(v_perif)


        # define cartesian coordinates
        self.cart = {'x':r[0],
                     'y':r[1],
                     'z':r[2],
                     'vx':v[0],
                     'vy':v[1],
                     'vz':v[2]}

    def draw(self):
        TA_span = np.linspace(0, 360, 500)
        R = []

        for t in TA_span:
            coord = [self.kep['h'], self.kep['incl'], self.kep['RA'], self.kep['e'], self.kep['w'], t]
            R.append(Orbit(coord, "keplerian", mu=self.mu).getCart())
            #print(Orbit(coord, "keplerian", mu=self.mu).getKep())
            #print(coord)
        R = np.array(R).T

        ax = plt.axes(projection='3d')
        ax.plot(R[0, :], R[1, :], R[2, :])

        return ax

    def __str__(self):

        if self.type == 'keplerian':
            return f'Keplerian coordinates: \n' \
                   f'[h: {self.kep["h"]}, incl: {self.kep["incl"]}, ' \
                   f'RA: {self.kep["RA"]}, e: {self.kep["e"]}, ' \
                   f'w: {self.kep["w"]}, TA: {self.kep["TA"]}]'
        else:
            return f'Cartesian coordinates: \n' \
                   f'[x: {self.cart["x"]}, y: {self.cart["y"]}, ' \
                   f'z: {self.cart["z"]}, vx: {self.cart["vx"]}, ' \
                   f'vy: {self.cart["vy"]}, vz: {self.cart["vz"]}]'


# =============================================================================
# Begin tests
# =============================================================================
class TestConversion(unittest.TestCase):

    def setUp(self):
        self.cart_coords_1 = [-6045, -3490, 2500, -3.457, 6.618, 2.533]
        self.kep_coords_1 = [58311, 153.2, 255.3, 0.17, 20.07, 28.45]

        self.cart_coords_2 = [-4040, 4815, 3629, -10.39, -4.772, 1.744]
        self.kep_coords_2 = [80000, 30, 40, 1.4, 60, 30]


    def test_conversion_to_kep(self):
        self.orbit = Orbit(self.cart_coords_1, "cartesian")
        np.testing.assert_almost_equal(self.orbit.getKep(), self.kep_coords_1, 0)

    def test_conversion_to_kep(self):
        self.orbit = Orbit(self.kep_coords_2, "keplerian")
        #print(self.orbit.getCart())
        np.testing.assert_almost_equal(self.orbit.getCart(), self.cart_coords_2, 0)

    def testSMA(self):
        self.orbit = Orbit(self.cart_coords_1, 'cartesian')
        np.testing.assert_almost_equal(self.orbit.getSemiMajorAxes(), 8788, 1)

    def testPeriod(self):
        self.orbit = Orbit(self.cart_coords_1, 'cartesian')
        np.testing.assert_almost_equal(self.orbit.getPeriod(), 2.278*3600, -1)

    def testString(self):
        self.orbit1 = Orbit(self.kep_coords_1, 'keplerian')
        self.orbit2 = Orbit(self.cart_coords_1, 'cartesian')
        print(self.orbit1)
        print(self.orbit2)


if __name__ == '__main__':
    unittest.main()

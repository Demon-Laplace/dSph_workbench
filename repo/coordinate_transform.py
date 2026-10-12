import numpy as np
from astropy import units as u
from astropy.coordinates import Galactocentric, Galactic, SkyCoord


class CoordinateTransform:
    """Handle coordinate transformations."""

    @staticmethod
    def to_heliocentric(x, y, z, vx=None, vy=None, vz=None):
        galactocentric_frame = Galactocentric()
        x_sun = galactocentric_frame.galcen_distance.to(u.kpc).value
        z_sun = galactocentric_frame.z_sun.to(u.kpc).value

        xh = x - x_sun
        yh = y
        zh = z - z_sun

        if all(v is not None for v in [vx, vy, vz]):
            solar_velocity = galactocentric_frame.galcen_v_sun
            # Astropy <=7 exposes this vector as a CartesianDifferential,
            # while Astropy 8 exposes an equivalent CartesianRepresentation.
            components = solar_velocity.d_xyz if hasattr(solar_velocity, "d_xyz") else solar_velocity.xyz
            vsun = components.to(u.km / u.s).value
            vxh = vx - vsun[0]
            vyh = vy - vsun[1]
            vzh = vz - vsun[2]
            return xh, yh, zh, vxh, vyh, vzh

        return xh, yh, zh

    @staticmethod
    def to_galactic(x, y, z, chunk_size=100000):
        if len(x) == 0:
            return np.array([], dtype=float), np.array([], dtype=float)

        l_list = []
        b_list = []

        for i in range(0, len(x), chunk_size):
            chunk_x = x[i:i + chunk_size]
            chunk_y = y[i:i + chunk_size]
            chunk_z = z[i:i + chunk_size]

            galactic_coords = Galactocentric(
                x=chunk_x * u.kpc,
                y=chunk_y * u.kpc,
                z=chunk_z * u.kpc,
            ).transform_to(Galactic())

            l_list.append(galactic_coords.l.deg)
            b_list.append(galactic_coords.b.deg)

        return np.concatenate(l_list), np.concatenate(b_list)

    @staticmethod
    def calculate_distances(x, y, z):
        return np.sqrt(
            x.astype(np.float32) ** 2 +
            y.astype(np.float32) ** 2 +
            z.astype(np.float32) ** 2
        )

    @staticmethod
    def galactic_to_equatorial(l, b, chunk_size=100000):
        if len(l) == 0:
            return np.array([], dtype=float), np.array([], dtype=float)

        ra_list = []
        dec_list = []

        for i in range(0, len(l), chunk_size):
            chunk_l = l[i:i + chunk_size]
            chunk_b = b[i:i + chunk_size]

            coords = SkyCoord(l=chunk_l * u.degree, b=chunk_b * u.degree, frame='galactic')
            eq_coords = coords.transform_to('icrs')

            ra_list.append(eq_coords.ra.deg)
            dec_list.append(eq_coords.dec.deg)

        return np.concatenate(ra_list), np.concatenate(dec_list)

import os
import numpy as np
import scipy as sp
from scipy import signal
import warnings
import csv
import random
from astropy.io import fits
from astropy.wcs import WCS

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


class ClusterSimulator:

    def __init__(self):

        self.n_clusters = 3
        self.image_size = 2048
        self.pixel_scale = 0.2
        self.n_cluster_galaxies = 20
        self.n_background_galaxies = 200
        self.n_stars = 30
        self.background_noise = 1.1e-12
        self.gain = 52062244760000.0
        self.simulation_dir = "simulations"
        self.parameters_dir = "simulations/parameters"
        self.setup_directories()


    def setup_directories(self):

        """Create necessary directories"""

        os.makedirs(self.simulation_dir, exist_ok=True)
        os.makedirs(self.parameters_dir, exist_ok=True)


    def create_empty_fits(self, size=2048):

        """Create an empty FITS file with proper WCS headers"""

        data = np.zeros((size, size), dtype=np.float32)
        hdu = fits.PrimaryHDU(data)

        # Add standard astronomical headers
        hdu.header['BITPIX'] = -32
        hdu.header['NAXIS'] = 2
        hdu.header['NAXIS1'] = size
        hdu.header['NAXIS2'] = size
        hdu.header['BSCALE'] = 1.0
        hdu.header['BZERO'] = 0.0
        hdu.header['EXPTIME'] = 1000.0
        hdu.header['GAIN'] = self.gain
        hdu.header['RDNOISE'] = 5.0
        hdu.header['TELESCOP'] = 'SIMULATED'
        hdu.header['INSTRUME'] = 'SIMCAM'
        hdu.header['FILTER'] = 'r'

        # Add WCS information
        hdu.header['CTYPE1'] = 'RA---TAN'
        hdu.header['CTYPE2'] = 'DEC--TAN'
        hdu.header['CRVAL1'] = 150.0  # RA center
        hdu.header['CRVAL2'] = 2.0  # DEC center
        hdu.header['CRPIX1'] = size / 2
        hdu.header['CRPIX2'] = size / 2
        hdu.header['CD1_1'] = -self.pixel_scale / 3600.0  # deg/pixel
        hdu.header['CD1_2'] = 0.0
        hdu.header['CD2_1'] = 0.0
        hdu.header['CD2_2'] = self.pixel_scale / 3600.0  # deg/pixel

        return hdu


    def get_plane_distance(self, data, xc, yc):

        nx, ny = len(data[0]), len(data)
        x, y = np.ogrid[-yc:float(ny) - yc, -xc:float(nx) - xc]

        return np.sqrt(x ** 2. + y ** 2.)


    def gaussian(self, r, sigma, mu0):

        return np.sqrt(mu0 ** 2.) * np.exp(-0.5 * (r / sigma) ** 2.)


    def exponential(self, r, h, mu0):

        return mu0 * np.exp(-r / h)


    def moffat(self, r, I0, alfa, beta):

        mof = np.sqrt(I0 ** 2) * (1 + (r / alfa) ** 2.) ** (-beta)

        return mof


    def analytic_ocam_psf(self, r, mu0, mag=None):

        r = r * 5.
        coef = np.array([2.10164715e-02, 2.30883062e-11])
        p = np.poly1d(coef)

        if mag:

            mu0 = p(10. ** (mag / (-2.5)))

        else:

            mu0 = 1.

        gauss_mu0, gauss_sigma = 0.061, 0.757 * 5.
        mof_I0, mof_alpha, mof_beta = 0.938, 0.635 * 5., 1.60
        exponential_h, exponential_I0 = 74.26 * 5., 6.022e-06

        gaussian_mod = self.gaussian(r, gauss_mu0, gauss_sigma) * mu0
        moffat_mod = self.moffat(r, mof_I0, mof_alpha, mof_beta) * mu0
        exponential_mod = self.exponential(r, exponential_h, exponential_I0) * mu0
        psf = gaussian_mod + moffat_mod + exponential_mod

        return psf


    def surface_brightness(self, flux, pix_scale=0.2):

        return -2.5 * np.log10(flux) + 2.5 * np.log10(pix_scale ** 2.)


    def find_ocam_psf_sfb_rad(self, mag, sfb):

        r = 0

        while r <= 8000:

            r += 1.
            sb = self.analytic_ocam_psf(r / 5., 0., mag=mag)

            if self.surface_brightness(sb) > sfb:
                return r

        return 0


    def sersic(self, mue, n, r, reff, flux=False, mag=None):

        if n > 0.36:

            k = (2.0 * n - 1.0 / 3. + 4.0 / 405. / n + 46.0 / 25515. / n ** 2.
                 + 131.0 / 1148175. / n ** 3. + 2194697. / 30690717750.0 / n ** 4.)

        else:

            k = 0.01945 - 0.8902 * n + 10.95 * n ** 2. - 19.67 * n ** 3. + 13.43 * n ** 4.

        if not flux:

            profile = mue + 2.5 * k * (r / reff) ** (1. / n)

        else:

            if mag is None:

                profile = mue * np.exp(-k * ((r / reff) ** (1. / n) - 1.))

            else:

                mue = (10. ** (-mag / 2.5) /
                       (2. * np.pi * reff ** 2. * np.exp(k) * n * k ** (-2. * n) *
                        sp.special.gamma(2. * n)))
                profile = mue * np.exp(-k * ((r / reff) ** (1. / n) - 1.))

        return profile


    def make_galaxy(self, mue, reff, n, nx, arat, pa):

        yc, xc = nx // 2, nx // 2
        ny = nx
        x, y = np.ogrid[-yc:float(ny) - yc, -xc:float(nx) - xc]
        angles = np.arctan2(y, x)
        distance = np.sqrt(x ** 2. + y ** 2.) / reff

        A = 1.
        B = arat
        edist = distance / (B / np.sqrt((B * np.cos(angles - pa)) ** 2. + (A * np.sin(angles - pa)) ** 2.))
        edist[yc, xc] = 0.

        r = edist.reshape(nx * nx)
        flux = self.sersic(mue, n, r, 1., flux=True)
        flux.resize(nx, nx)

        return flux


    def generate_objects(self):

        """Generate cluster galaxies, background galaxies, and stars"""

        maxx = self.image_size - 1
        maxy = self.image_size - 1

        bg_data = []

        for i in range(self.n_background_galaxies):

            x = int(np.random.random() * maxx)
            y = int(np.random.random() * maxy)
            mue = 10. ** (-11. - np.random.random() * 3.)
            reff = (0.5 + np.random.random() * 3.) * 5.
            n = 2. + np.random.random() * 2.
            nx = int(10 * reff * n)
            arat = 0.3 + np.random.random() * 0.7
            pa = np.random.random() * np.pi
            bg_data.append((x, y, mue, reff, n, nx, arat, pa))

        cg_data = []

        for i in range(self.n_cluster_galaxies):

            x = int(np.random.random() * maxx)
            y = int(np.random.random() * maxy)
            mue = 10. ** (-11. - np.random.random() * 4.)
            reff = (1. + np.random.random() * 99.) * 5.
            n = 0.5 + np.random.random() * 1.5
            nx = int(10 * reff * n)
            arat = 0.3 + np.random.random() * 0.7
            pa = np.random.random() * np.pi
            cg_data.append((x, y, mue, reff, n, nx, arat, pa))

        st_data = []

        for i in range(self.n_stars):

            mag = 26 - np.random.exponential(scale=3)
            x = int(np.random.random() * maxx)
            y = int(np.random.random() * maxy)
            st_data.append((x, y, mag))

        return (sorted(bg_data, key=lambda t: t[2]),
                sorted(cg_data, key=lambda t: t[2]),
                sorted(st_data, key=lambda t: t[2]))

    def run_simulation(self):

        """Run the complete simulation"""

        print("Starting cluster simulation...")

        bg_data, cg_data, st_data = self.generate_objects()

        for cluster_idx in range(self.n_clusters):

            print(f"\n=== Creating Cluster {cluster_idx + 1} ===")

            hdu = self.create_empty_fits(self.image_size)
            data = np.zeros((self.image_size, self.image_size), dtype=np.float32)
            maxx = self.image_size - 1
            maxy = self.image_size - 1

            cluster_params = []
            background_params = []
            star_params = []

            rand_indexes = random.sample(range(len(bg_data)), min(1400, len(bg_data)))
            skip_first = rand_indexes[0:149]
            skip_second = rand_indexes[150:549]
            skip_third = rand_indexes[550:min(1149, len(rand_indexes))]

            print("Adding cluster galaxies...")

            for gal in cg_data:

                x, y, base_mue, base_reff, n, nx, arat, pa = gal


                if cluster_idx == 1:

                    mu = base_mue * random.uniform(0.75, 0.85)
                    reff = base_reff * random.uniform(0.65, 0.75)

                elif cluster_idx == 2:

                    mu = base_mue * random.uniform(0.75, 0.85)
                    reff = base_reff * random.uniform(0.65, 0.75)

                else:

                    mu = base_mue
                    reff = base_reff

                galaxy = self.make_galaxy(mu, reff, n, nx, arat, pa)

                xe, ye = galaxy.shape
                x0 = max(0, x - xe // 2)
                x1 = min(maxx, x + xe // 2)
                y0 = max(0, y - ye // 2)
                y1 = min(maxy, y + ye // 2)

                gal_slice = galaxy[(y0 - y + ye // 2):(y1 - y + ye // 2),
                (x0 - x + xe // 2):(x1 - x + xe // 2)]

                data[y0:y1, x0:x1] += gal_slice

                cluster_params.append([x, y, mu, reff, arat, n, pa])

            print("Adding background galaxies...")

            for idx, gal in enumerate(bg_data):

                x, y, base_mue, base_reff, n, nx, arat, pa = gal

                if cluster_idx == 1 and idx in skip_second:

                    continue

                elif cluster_idx == 2 and idx in skip_third:

                    continue

                elif cluster_idx == 0 and idx in skip_first:

                    continue

                if cluster_idx == 1:

                    mu = base_mue * 0.94
                    reff = base_reff * 0.92

                elif cluster_idx == 2:

                    mu = base_mue * 0.80
                    reff = base_reff * 0.64

                else:

                    mu = base_mue
                    reff = base_reff

                galaxy = self.make_galaxy(mu, reff, n, nx, arat, pa)

                xe, ye = galaxy.shape
                x0 = max(0, x - xe // 2)
                x1 = min(maxx, x + xe // 2)
                y0 = max(0, y - ye // 2)
                y1 = min(maxy, y + ye // 2)

                gal_slice = galaxy[(y0 - y + ye // 2):(y1 - y + ye // 2),
                (x0 - x + xe // 2):(x1 - x + xe // 2)]

                data[y0:y1, x0:x1] += gal_slice

                background_params.append([x, y, mu, reff, arat, n, pa])

            print("Convolving with PSF...")

            h = self.get_plane_distance(np.zeros([81, 81]), 40, 40) * 0.2
            h = self.analytic_ocam_psf(h, 1, mag=None)
            h = h / np.sum(h)
            data = signal.fftconvolve(data, h, mode='same')
            data[data < 0] = 0.

            print("Adding stars...")

            for star in st_data:

                x, y, mag = star
                rad = int(self.find_ocam_psf_sfb_rad(mag, self.surface_brightness(5e-14)))

                if rad <= 0:

                    continue

                star_profile = np.zeros([2 * rad, 2 * rad])
                star_profile = self.get_plane_distance(star_profile, rad, rad) * 0.2
                star_profile = self.analytic_ocam_psf(star_profile, 1, mag=mag)

                xe, ye = star_profile.shape
                x0 = max(0, x - xe // 2)
                x1 = min(maxx, x + xe // 2)
                y0 = max(0, y - ye // 2)
                y1 = min(maxy, y + ye // 2)

                star_slice = star_profile[(y0 - y + ye // 2):(y1 - y + ye // 2),
                (x0 - x + xe // 2):(x1 - x + xe // 2)]

                data[y0:y1, x0:x1] += star_slice
                star_params.append([x, y, mag])

            print("Adding noise...")

            data = np.random.poisson(data * self.gain) / self.gain
            data += np.random.normal(loc=0., scale=self.background_noise, size=data.shape)

            hdu.data = data.astype('float32')
            output_path = os.path.join(self.simulation_dir, f'cluster{cluster_idx + 1}.fits')
            hdu.writeto(output_path, overwrite=True)
            print(f"Saved: {output_path}")

            self.save_parameters(cluster_idx + 1, cluster_params, background_params, star_params)

        print(f"\nSimulation completed! Check the '{self.simulation_dir}' directory.")

    def save_parameters(self, cluster_id, cluster_params, background_params, star_params):

        """Save object parameters to CSV files"""

        if cluster_params:

            with open(os.path.join(self.parameters_dir, f'cluster{cluster_id}_galaxies.csv'), 'w') as f:

                writer = csv.writer(f)
                writer.writerow(['x', 'y', 'mu', 'reff', 'axis_ratio', 'sersic_n', 'pa'])
                writer.writerows(cluster_params)

        if background_params:

            with open(os.path.join(self.parameters_dir, f'cluster{cluster_id}_background.csv'), 'w') as f:

                writer = csv.writer(f)
                writer.writerow(['x', 'y', 'mu', 'reff', 'axis_ratio', 'sersic_n', 'pa'])
                writer.writerows(background_params)

        if star_params:

            with open(os.path.join(self.parameters_dir, f'cluster{cluster_id}_stars.csv'), 'w') as f:

                writer = csv.writer(f)
                writer.writerow(['x', 'y', 'mag'])
                writer.writerows(star_params)


def main():

    simulator = ClusterSimulator()
    simulator.run_simulation()


if __name__ == "__main__":
    main()


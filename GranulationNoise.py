import numpy as np

def convective_blueshift(T, g):
    """
    Compute the convective blueshift of the star with respect to solar convective blueshift.

    Parameters
    ----------
    T : array_like
        Temperature of the star in Kelvin.
    g : array_like
        Surface gravity of the star in cm/s2.

    Returns
    -------
    Convective Blueshift : array_like
        Convective Blueshift of the star with respect to solar convective blueshift.
    """
    Tsun = 5772
    gsun = 275.4

    star_velo = (T ** (32 / 9)) * (g ** (-2 / 9))
    sun_velo = (Tsun ** (32 / 9)) * (gsun ** (-2 / 9))
    v = star_velo / sun_velo
    star_i = (54.98 * np.log10(T)) - (4.80 * np.log10(g)) - 169.0
    sun_i = (54.98 * np.log10(Tsun)) - (4.80 * np.log10(gsun)) - 169.0
    irms = star_i / sun_i
    conv_blue = (irms ** 2) * (v)
    
    return conv_blue

def convective_blueshift_error(T, g, errT, errg):
    ''' This function calculates the error in convective blueshift of the star with respect to solar convective blueshift'''
    Tsun = 5772
    gsun = 275.4

    star_velo = (T ** (32 / 9)) * (g ** (-2 / 9))
    sun_velo = (Tsun ** (32 / 9)) * (gsun ** (-2 / 9))
    star_i = (54.98 * np.log10(T)) - (4.80 * np.log10(g)) - 169.0
    sun_i = (54.98 * np.log10(Tsun)) - (4.80 * np.log10(gsun)) - 169.0
    d = (sun_i ** 2) * (sun_velo)
    conv_blue = ((star_i ** 2) * (star_velo)) / d

    sigma_v = star_velo * (np.sqrt((((32 * errT) / (9 * T)) ** 2) + (((2 * errg) / (9 * g)) ** 2)))
    sigma_i = np.sqrt(((54.98 * errT) / (np.log(10) * T)) ** 2 + ((4.80 * errg) / (np.log(10) * g)) ** 2)
    sigma_cb = (conv_blue) * (np.sqrt((sigma_v / star_velo) ** 2 + ((2 * sigma_i) / (star_i)) ** 2))
    err_CB_shift = sigma_cb

    return err_CB_shift

def radial_velocity_dispersion(CB, T, M, R):
    ''' This function calculates the dispersion in radial velocity due to granulation of the star with respect to solar value'''
    ngsun = (1 / (5772 * 1) ** 2)
    sigmasun = (1) / (np.sqrt(ngsun))
    ng = (M / (T * R)) ** 2
    sigmarv = (CB) / (np.sqrt(ng))
    sigma_rv = sigmarv / sigmasun

    return sigma_rv

def radial_velocity_dispersion_error(CB, T, M, R, errCB, errT, errM, errR):
    ''' This function calculates the error in dispersion in radial velocity due to granulation of the star with respect to solar value'''
    ngsun = (1 / (5772 * 1) ** 2)
    sigmasun = (1) / (np.sqrt(ngsun))
    ng = (M / (T * R)) ** 2
    sigmarv = (CB) / (np.sqrt(ng))
    sigmarvsigma = sigmarv * (np.sqrt((errCB / CB) ** 2 + (errT / T) ** 2 + (errM / M) ** 2 + (errR / R) ** 2))
    err_sigma_rv = sigmarvsigma / sigmasun

    return err_sigma_rv

# Constants
STAR_DATA = {
    "Alpha Cen B": {"M": 0.91, "R": 0.859, "Teff": 5248, "ge": 4.55},
    # "HD 166620": {"M": 0.76, "R": 0.77, "Teff": 4989, "ge": 4.65},
    # "Alpha Cen A": {"M": 1.11, "R": 1.22, "Teff": 5790, "ge": 4.34},
    # "Tau Ceti": {"M": 1.2, "R": 1.22, "Teff": 6220, "ge": 4.496}
}

def main():
    for star, data in STAR_DATA.items():
        M, R, Teff, ge = data["M"], data["R"], data["Teff"], data["ge"]
        g = (10 ** ge) / 100

        CB = convective_blueshift(Teff, g)
        CB_error = convective_blueshift_error(Teff, g, 48, 0.001)
        rv_rms = radial_velocity_dispersion(CB, Teff, M, R)
        rv_rms_error = radial_velocity_dispersion_error(CB, Teff, M, R, CB_error, 64, 0.02, 0.04)

        print(f"The value of the convective blueshift for {star} is: {CB * 350} m/s")
        print(f"The error on the convective blueshift for {star} is: {CB_error * 350} m/s")
        print(f"The value of the sigma RV for {star} is: {rv_rms * 0.40} m/s")
        print(f"The error on the sigma RV for {star} is: {rv_rms_error * 0.40} m/s")
        print()

if __name__ == "__main__":
    main()

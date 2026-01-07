import astropy.units as u
import numpy as np
from pyirf.spectral import (
    CRAB_HEGRA,
    CRAB_MAGIC_JHEAP2015,
    DIFFUSE_FLUX_UNIT,
    POINT_SOURCE_FLUX_UNIT,
    DAMPE_P_He_SPECTRUM,
    LogParabola,
    PowerLaw,
)

__all__ = [
    "CRAB_HEGRA",
    "CRAB_MAGIC_JHEAP2015",
    "DAMPE_P_He_SPECTRUM",
    "MRK421_VERITAS_VERY_LOW",
    "MRK421_VERITAS_LOW",
    "MRK421_VERITAS_HIGH",
    "MRK501_TACTIC",
    "MRK501_VERITAS_LOW",
    "MRK501_VERITAS_HIGH",
    "LHAASO_J2226_ASg",
    "LHAASO_J2226_HAWC",
    "LHAASO_J2108",
    "LHAASO_J1908",
    "ES1959_HEGRA_LOW",
    "ES1959_HEGRA_HIGH",
    "POINT_SOURCE_FLUX_UNIT",
    "DIFFUSE_FLUX_UNIT",
    "CRAB_LHAASO_2022",
    "CRAB_MAGIC_2020",
    "CRAB_HAWC_2019",
    "CRAB_HAWC_2019_NN",
    "CRAB_VERITAS_2015",
    "HESS_J1702",
    "CRAB_HESS_2024",
    "CRAB_LST_2023"
]

# Surprisingly, ECPL class is not defined in pyirf, so we define it here
class PowerLawExpCutoff:

    @u.quantity_input(
        normalization=[DIFFUSE_FLUX_UNIT, POINT_SOURCE_FLUX_UNIT], e_ref=u.TeV
    )
    def __init__(self, normalization, index, e_cutoff, e_ref=1 * u.TeV):
        self.normalization = normalization
        self.index = index
        self.e_cutoff = e_cutoff
        self.e_ref = e_ref

    @u.quantity_input(energy=u.TeV)
    def __call__(self, energy):
        e = (energy / self.e_ref).to_value(u.one)
        ec = (energy / self.e_cutoff).to_value(u.one)
        return self.normalization * e ** self.index * np.exp(-ec)

    def __repr__(self):
        return f"{self.__class__.__name__}({self.normalization} * (E / {self.e_ref})**{self.index} * exp (-E / {self.e_cutoff}))"

# LogParabola in pyirf is with log_10, while in mane papers, it is defined as log_e
class LogParabola_natural:
    r"""
    A log parabola flux parameterization with natural logaritm

    .. math::

        \Phi(E, \Phi_0, \alpha, \beta, E_\text{ref}) =
        \Phi_0 \left(
            \frac{E}{E_\text{ref}}
        \right)^{\alpha + \beta \cdot \log_{10}(E / E_\text{ref})}

    Attributes
    ----------
    normalization: astropy.units.Quantity[flux]
        :math:`\Phi_0`,
    a: float
        :math:`\alpha`
    b: float
        :math:`\beta`
    e_ref: astropy.units.Quantity[energy]
        :math:`E_\text{ref}`
    """

    @u.quantity_input(
        normalization=[DIFFUSE_FLUX_UNIT, POINT_SOURCE_FLUX_UNIT], e_ref=u.TeV
    )
    def __init__(self, normalization, a, b, e_ref=1 * u.TeV):
        """Create a new LogParabola spectrum"""
        self.normalization = normalization
        self.a = a
        self.b = b
        self.e_ref = e_ref


    @u.quantity_input(energy=u.TeV)
    def __call__(self, energy):
        e = (energy / self.e_ref).to_value(u.one)
        return self.normalization * e ** (self.a + self.b * np.log(e))



    def __repr__(self):
        return f"{self.__class__.__name__}({self.normalization} * (E / {self.e_ref})**({self.a} + {self.b} * log10(E / {self.e_ref}))"





# ECPL parametrization of the Mrk 421
# VERITAS 2011, very low, low and high-B state: https://iopscience.iop.org/article/10.1088/0004-637X/738/1/25/pdf
# very-low state is in the best agreement with other experiments, including HAWC longterm SED
# HAWC 2022, https://iopscience.iop.org/article/10.3847/1538-4357/ac58f6/pdf
MRK421_VERITAS_VERY_LOW = PowerLawExpCutoff(
    normalization=4.48e-11 / (u.TeV * u.cm ** 2 * u.s), index=-2.29, e_cutoff=1.59 * u.TeV, e_ref = 1 * u.TeV
)
MRK421_VERITAS_LOW = PowerLawExpCutoff(
    normalization=7.60e-11 / (u.TeV * u.cm ** 2 * u.s), index=-2.285, e_cutoff=2.95 * u.TeV, e_ref = 1 * u.TeV
)
MRK421_VERITAS_HIGH = PowerLawExpCutoff(
    normalization=22.23e-11 / (u.TeV * u.cm ** 2 * u.s), index=-1.88, e_cutoff=3.06 * u.TeV, e_ref = 1 * u.TeV
)

# PL parametrization of the Mrk 501
# TACTIC 2008, representing probably low or mid state: https://iopscience.iop.org/article/10.1088/0954-3899/35/6/065202
# High energy cutoff not seen by other experiments as well (HAWC, MAGIC, VERITAS, ARGO)
MRK501_TACTIC = PowerLaw(
    normalization=1.66e-11 / (u.TeV * u.cm ** 2 * u.s), index=-2.80, e_ref = 1 * u.TeV
)

# VERITAS_LOW, 2009 March 24 (54914)
# low state, https://iopscience.iop.org/article/10.1088/0004-637X/729/1/2/pdf
# Fermi Flux from LCR (0.1-100 GeV): 4.21e-8 GeV/cm2/s2
MRK501_VERITAS_LOW = PowerLaw(
    normalization=5.78e-12 / (u.TeV * u.cm ** 2 * u.s), index=-2.72, e_ref = 1 * u.TeV
)

# VERITAS_HIGH
# high state (non flaring), https://iopscience.iop.org/article/10.1088/0004-637X/727/2/129/pdf
MRK501_VERITAS_HIGH = PowerLaw(
    normalization=0.88e-11 / (u.TeV * u.cm ** 2 * u.s), index=-2.26, e_ref = 1 * u.TeV
)

# Log-parabola parametrization of LHAASO J1908+0621 (eHWC J1907+063)
# HAWC flux points: https://arxiv.org/abs/1909.08609
LHAASO_J1908 = LogParabola(
    normalization=0.95e-13 / (u.TeV * u.cm ** 2 * u.s), a=-2.46, b=-0.11, e_ref = 10 * u.TeV
)

# ECPL parametrization of LHAASO J2108+5157
# LHAASO flux points + LST-1 ULs: LST col 2023, https://arxiv.org/abs/2210.00775
LHAASO_J2108 = PowerLawExpCutoff(
    normalization=7.6e-14 / (u.TeV * u.cm ** 2 * u.s), index=-1.37, e_cutoff=50 * u.TeV, e_ref = 1 * u.TeV
)

# LHAASO J2226+6057 (Boomerang, G106.3+2.7)
# Tibet ASg flux points above 6 TeV can be described with PL https://www.nature.com/articles/s41550-020-01294-9
LHAASO_J2226_ASg = PowerLaw(
    normalization=9.5e-16 / (u.TeV * u.cm ** 2 * u.s), index=-2.95, e_ref = 40 * u.TeV
)
# VERITAS+HAWC joint fit PL https://iopscience.iop.org/article/10.3847/2041-8213/ab96cc/pdf
LHAASO_J2226_HAWC = PowerLaw(
    normalization=2.46e-15 / (u.TeV * u.cm ** 2 * u.s), index=-2.29, e_ref = 20 * u.TeV
)

# ECPL parametrization of 1ES 1959+650
# HEGRA 2003, high and low state: https://www.aanda.org/articles/aa/pdf/2003/28/aafe161.pdf
ES1959_HEGRA_LOW = PowerLawExpCutoff(
    normalization=6.0e-12 / (u.TeV * u.cm ** 2 * u.s), index=-1.8, e_cutoff=2.7 * u.TeV, e_ref = 1 * u.TeV
)
ES1959_HEGRA_HIGH = PowerLawExpCutoff(
    normalization=5.6e-11 / (u.TeV * u.cm ** 2 * u.s), index=-1.83, e_cutoff=4.2 * u.TeV, e_ref = 1 * u.TeV
)

# Log-parabola parametrization of CRAB (WCDA+KM2A of LHAASO)
# https://www.science.org/doi/epdf/10.1126/science.abg5137
CRAB_LHAASO_2022 = LogParabola(
    normalization=8.2e-14 / (u.TeV * u.cm ** 2 * u.s), a=-2.9, b=-0.19, e_ref = 10 * u.TeV
)

# Log-parabola parametrization of CRAB (MAGIC LZA, up to 100 TeV)
# https://www.aanda.org/articles/aa/pdf/2020/03/aa36899-19.pdf 
CRAB_MAGIC_2020 = LogParabola(
    normalization=2.95e-23 / (u.eV * u.cm ** 2 * u.s), a=-2.48, b=-0.23, e_ref = 1 * u.TeV
)

# Log-parabola parametrization of CRAB (HAWC, "ground parameter")
# https://iopscience.iop.org/article/10.3847/1538-4357/ab2f7d/pdf 
CRAB_HAWC_2019 = LogParabola(
    normalization=2.35e-13 / (u.TeV * u.cm ** 2 * u.s), a=-2.79, b=-0.10, e_ref = 7 * u.TeV
)

# Log-parabola parametrization of CRAB (HAWC, "neural network")
# https://iopscience.iop.org/article/10.3847/1538-4357/ab2f7d/pdf 
CRAB_HAWC_2019_NN = LogParabola(
    normalization=2.31e-13 / (u.TeV * u.cm ** 2 * u.s), a=-2.73, b=-0.06, e_ref = 7 * u.TeV
)

# Log-parabola parametrization of CRAB (VERITAS)
# https://arxiv.org/pdf/1508.06442.pdf 
CRAB_VERITAS_2015 = LogParabola(
    normalization=3.75e-11 / (u.TeV * u.cm ** 2 * u.s), a=-2.467, b=-0.16, e_ref = 1 * u.TeV
)

# Log-parabola parametrization of CRAB (HESS 2024)
# https://arxiv.org/pdf/2403.12608
# DONT USE IT, there is a typo in the paper, the parametrization is wrong
CRAB_HESS_2024 = LogParabola_natural(
    normalization=4.06e-11 / (u.TeV * u.cm ** 2 * u.s), a=-2.530, b=-0.086, e_ref = 1 * u.TeV
)

# Log-parabola parametrization of CRAB (LST 2023)
# https://arxiv.org/pdf/2306.12960 
CRAB_LST_2023 = LogParabola_natural(
    normalization=3.05e-10 / (u.TeV * u.cm ** 2 * u.s), a=-2.25, b=-0.114, e_ref = 400 * u.GeV
)

# HESS J1702−420A
# https://arxiv.org/pdf/2402.03511.pdf
# https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.124.021102 
HESS_J1702 = PowerLaw(
    normalization=1.6e-13 / (u.TeV * u.cm ** 2 * u.s), index=-1.53, e_ref = 1 * u.TeV
)
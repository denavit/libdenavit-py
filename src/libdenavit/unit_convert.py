from typing import Union

import pint

ureg = pint.UnitRegistry()
units_for_pint = {
    ## Length
    "in": ureg.inch, # inch
    "ft": ureg.ft,   # foot
    "m" : ureg.m,    # meter
    "cm": ureg.cm,   # centimeter
    "mm": ureg.mm,   # milimeter

    ## Area
    "sqin": ureg.inch ** 2, # inch^2
    "sqft": ureg.ft ** 2,   # foot^2
    "sqm": ureg.m ** 2,    # meter^2
    "sqcm": ureg.cm ** 2,  # centimeter^2
    "sqmm": ureg.mm ** 2,  # milimeter^2

    ## Curvature
    "1/in": 1/ureg.inch, # 1/inch
    "1/mm": 1/ureg.mm,   # 1/millimeter
    "1/m" : 1/ureg.m,    # 1/meter

    ## Time
    "s"     : ureg.second, # second
    "second": ureg.second, # second
    "min"   : ureg.minute, # minute
    "h"     : ureg.hour,   # hour
    "hr"    : ureg.hour,   # hour

    ## Force
    "kip"      : ureg.kip,              # kip
    "kips"     : ureg.kip,              # kip
    "lbf"      : ureg.lbf,              # pound-force
    "kn"       : ureg.knewton,          # kilonewton
    "mn"       : ureg.meganewton,       # meganewton
    "n"        : ureg.newton,           # kilonewton
    "metricton": ureg.force_metric_ton, # tonne
    "tonne"    : ureg.force_metric_ton, # tonne
    "longton"  : ureg.force_long_ton,   # long ton
    "shortton" : ureg.force_ton,        # short ton

    ## Moment
    'kin'   : ureg.kforce_pound*ureg.inch,      # kilo pound-inch
    'k-in'  : ureg.kforce_pound*ureg.inch,      # kilo pound-inch
    'kip-in': ureg.kforce_pound*ureg.inch,      # kilo pound-inch
    'kft'   : ureg.kforce_pound*ureg.foot,      # kilo pound-foot
    'k-ft'  : ureg.kforce_pound*ureg.foot,      # kilo pound-foot
    'kip-ft': ureg.kforce_pound*ureg.foot,      # kilo pound-foot
    'lbft'  : ureg.foot_pound,                  # pound-foot
    'lb-ft' : ureg.foot_pound,                  # pound-foot
    'nmm'   : ureg.newton*ureg.millimeter,      # newton-milimeter
    'n-mm'  : ureg.newton*ureg.millimeter,      # newton-milimeter
    'knm'   : ureg.knewton*ureg.meter,          # kilo newton-milimeter
    'kn-m'  : ureg.knewton*ureg.meter,          # kilo newton-milimeter
    'mnm'   : ureg.mnewton*ureg.meter,          # millinewton meter
    'mn-m'  : ureg.mnewton*ureg.meter,          # millinewton meter
    'tfm'   : ureg.force_metric_ton*ureg.meter, # metric ton-force meter
    'tf-m'  : ureg.force_metric_ton*ureg.meter, # metric ton-force meter

    ## Pressure
    'psi'           : ureg.psi,                         # pound per inch^2
    'psf'           : ureg.force_pound/ureg.foot**2,    # pound per foot^2
    'ksi'           : ureg.kpsi,                        # killo pound per foot^2
    'mpa'           : ureg.megaPa,                      # megapascal
    'n/mm^2'        : ureg.megaPa,                      # newton per millimeter^2
    'kpa'           : ureg.kPa,                         # killopascal
    'gpa'           : ureg.gigaPa,                      # gigapascal
    'kn/cm^2'       : ureg.knewton/ureg.cm**2,          # killo-newton/centimeter^2
    'kgscm'         : ureg.kilogram_force/ureg.cm**2,   # Kilograms Per Square Centimeter
    'tscm'          : ureg.force_metric_ton/ureg.cm**2, # metric tonne per cm^2
    'metricton/cm^2': ureg.force_metric_ton/ureg.cm**2, # metric tonne per cm^2
    'longton/in^2'  : ureg.force_long_ton/ureg.inch**2, # long ton per inch^2
    'shortton/in^2' : ureg.force_ton/ureg.inch**2,      # short ton per inch^2

    ## Density
    'pcf' : ureg.pound/ureg.ft**3,      # pounds per cubic foot
    'kgcm': ureg.metric_ton/ureg.cm**3, # metric tonne per cm^2

    ## Angle
    'rad': ureg.radian, # radian
    'deg': ureg.degree, # degree

    ## Volume
    'cbin': ureg.inch**3,  # cubic inch
    'cbft': ureg.ft**3,    # cubic foot
    'cbyd': ureg.yard**3,  # cubic yard
    'cbm' : ureg.meter**3, # cubic meter
    'cbmm': ureg.mmeter**3, # cubic millimeter

    ## Area
    'in^2': ureg.inch ** 2,  # square inch
    'ft^2': ureg.ft ** 2,  # square foot
    'yd^2': ureg.yard ** 2,  # square yard
    'm^2' : ureg.meter ** 2,  # square meter
    'mm^2': ureg.mmeter ** 2,  # square millimeter
}


def unit_conversion_factor(old_unit: str, new_unit: str) -> float:
    return unit_convert(1, old_unit, new_unit)


def unit_convert(old_value: Union[float, int], old_unit: str, new_unit: str) -> float:
    """
    **Defined units to use:**
        Length Units: in, ft, m, cm, mm

        Curvature units: 1/in, 1/mm, 1/m

        Time Units: s, second, min, h, hr

        Force Units: kip, kips, lbf, kn, mn, n, metricton, tonne, longton, shortton

        Moment Units: kin, k-in, kip-in, kft, k-ft, kip-ft, lbft, lb-ft, nmm, n-mm,
        knm, kn-m, mnm, mn-m, tfm, tf-m

        Pressure Units: psi, psf, ksi, mpa, n/mm^2 kpa, gpa, kn/cm^2, kgcm, tscm, metricton/cm^2,
        longton/in^2, shortton/in^2

        Density Units: pcf, kgcm

        Angle Units: rad, deg

        Volume Units: cbin, cbft, cbyd, cbm, cbmm

        Area Units: in^2, ft^2, yd^2, m^2, mm^2
    """

    old_unit = old_unit.strip()
    new_unit = new_unit.strip()

    value_with_unit = old_value*units_for_pint[old_unit.lower()]
    new_value = value_with_unit.to(units_for_pint[new_unit.lower()]).magnitude

    return new_value

module MoistAirProp

import Roots
"""
    psat_water(T)

Calculates the saturation pressure of the water in ``kPa`` at a given temperature in ``Kelvin``.
`T` should range between 173.15 to 473.15.
Ref -> ASHRAE Handbook - Fundamentals - Chapter-1

"""
function psat_water(T)
    C1 = -5.6745359E+03
    C2 = 6.3925247
    C3 = -9.6778430E-03
    C4 = 6.2215701E-07
    C5 = 2.0747825E-09
    C6 = -9.4840240E-13
    C7 = 4.1635019
    C8 = -5.8002206E+03
    C9 = 1.3914993
    C10 = -4.8640239E-02
    C11 = 4.1764768E-05
    C12 = -1.4452093E-08
    C13 = 6.5459673
    if (T > 173.15 && T < 273.15)
        p_ws = exp(C1/T + C2 + C3*T + C4*T^2 + C5*T^3 + C6*T^4 + C7*log(T))
    elseif (T >= 273.15 && T < 473.15)
        p_ws = exp(C8/T + C9 + C10*T + C11*T^2 + C12*T^3 + C13*log(T))
    else
        error("Temperature should be between -100 and 200 deg C")
    end #if
    return p_ws/1000
end

"""
    humidity_ratio(p, t, str, val)

Calculates humidity ratio provided `p` total pressure in ``kPa``, `t` Temperature in kelvin, `str` string could be (``rh`` for relative humidty, ``wbt`` for wet bulb temperature, ``dpt`` for dewpoint temperature, ``w`` for humidity ratio), `val` value of variable specified in the string
"""
function humidity_ratio(p, t, str, val)
    t0 = 273.15
    if str == "rh" || str == "RH"
        rh = val
        p_ws = psat_water(t)
        w = 0.621945*rh*(p_ws)/(p-rh*p_ws)
    elseif str == "WBT" || str == "wbt"
        wbt = val
        p_ws_wbt = psat_water(wbt)
        ws_wbt = 0.621945*(p_ws_wbt)/(p-p_ws_wbt)
        if t >= t0
            w = ((2501-2.326*(wbt-t0))*ws_wbt - 1.006*((t-t0) - (wbt-t0)))/(2501 + 1.86*(t-t0) - 4.186*(wbt-t0))
        else
            w = ((2830-0.24*(wbt-t0))*ws_wbt - 1.006*((t-t0) - (wbt-t0)))/(2830 + 1.86*(t-t0) - 2.1*(wbt-t0))
        end #if
    elseif str == "dpt" || str == "DPT"
        dpt = val
        p_ws = psat_water(dpt)
        w = 0.621945*(p_ws)/(p-p_ws)
    elseif str == "w" || str == "W"
        w = val
    else
        error("use 'rh' for relative humidity,'wbt' for wet bulb temp,'dpt' for dew point temp")
    end #if
    return w
end #function

"""
    density_moistair(p,t,str,val)

Calculates density of moist air in ``kg/m^3`` provided `p` total pressure in ``kPa``, `t` Temperature in kelvin, `str` string could be (``rh`` for relative humidty, ``wbt`` for wet bulb temperature, ``dpt`` for dewpoint temperature, ``w`` for humidity ratio), `val` value of variable specified in the string
"""
function density_moistair(p,t,str,val)
    if str == "w" || str == "w"
        w = val
    else
        w = humidity_ratio(p, t, str, val)
    end
    v = 0.287042*(t)*(1 + 1.607858*w)/p
    return 1/v
end
"""
    dewpoint(p,t,str,val)

Calculates dewpoint temperature of moist air in ``kelvin`` provided `p` total pressure in ``kPa``, `t` Temperature in ``kelvin``, `str` string could be (``rh`` for relative humidty, ``wbt`` for wet bulb temperature, ``w`` for humidity ratio), `val` value of variable specified in the string
"""
function dewpoint(p,t,str,val)
    if str == "w" || str == "W"
        w = val
    else
        w = humidity_ratio(p, t, str, val)
    end
    p_w = p/(0.621945/w +1)
    t_dpt = Roots.find_zero(x -> (psat_water(x)-p_w),(t-100,t))
end
"""
    enthalpy_moistair(p,t,str,val)

Calculates enthalpy of moist air in ``kJ/kg`` provided `p` total pressure in ``kPa``, `t` Temperature in ``kelvin``, `str` string could be (``rh`` for relative humidty, ``wbt`` for wet bulb temperature, ``dpt`` for dewpoint temperature, ``w`` for humidity ratio), `val` value of variable specified in the string
"""
function enthalpy_moistair(p, t, str, val)
    if str == "W" || str == "w"
        w = val
        h = 1.006*(t-273.15)+ w*(2501+1.86*(t-273.15))
    else
        w = humidity_ratio(p, t, str, val)
        h = 1.006*(t-273.15)+ w*(2501+1.86*(t-273.15))
    end
    return h
end
"""
    wetbulb(p,t,str,val)

Calculates wetbulb temperature of moist air in ``Kelvin`` provided `p` total pressure in ``kPa``, `t` Temperature in ``kelvin``, `str` string could be (``rh`` for relative humidty, ``dpt`` for dewpoint temperature, ``w`` for humidity ratio), `val` value of variable specified in the string
"""
function wetbulb(p,t,str,val)
    w = humidity_ratio(p,t,str,val)
    if str == "dpt" || str == "DPT"
        dpt = val
    else
        dpt = dewpoint(p,t,str,val)
    end
    t_wbt = Roots.find_zero(temp -> (humidity_ratio(p,t,"wbt",temp)-w), (dpt,t))
end

export psat_water, humidity_ratio, density_moistair, dewpoint, enthalpy_moistair, wetbulb

end # module

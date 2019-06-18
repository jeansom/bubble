lines = {
    "O VI": {
        'continuum': ((1110, 1140), (1180, 1200)), # continuum bounds, [ang]
        'line': (1150, 1175), # line bounds, [ang]
        'p0': (0.4,1160,1) # initial values for gauss. fit (norm [erg/cm^2/s], mu [ang], sigma [ang])
    },
    "HI Ly$\\alpha$ (atmosphere)": {
        'continuum': ((1175, 1200), (1230, 1300)), 
        'line': (1200, 1230), 
        'p0': (53.22743067, 1216.1198963, 2.88220604)
    },
    "HI Ly$\\alpha$": {
        'continuum': ((1300, 1345), (1410, 1500)), 
        'line': (1360, 1375),
        'p0': (3, 1360, 1)
    },
    "N V": {
        'continuum': ((1300, 1345), (1410, 1500)), 
        'line': (1383, 1410), 
        'p0': (3, 1390, 1)
    },
    "C IV": {
        'continuum': ((1700, 1730), (1750, 1800)),
        'line': (1733, 1748),
        'p0': (3, 1740, 1)
    },
    "He II": {
        'continuum': ((1800, 1830), (1850, 1900)),
        'line': (1836, 1847),
        'p0': (3, 1840, 1)
    },
    "Si IV": {
        'continuum': ((1500, 1550), (1600, 1700)),
        'line': (1570, 1590),
        'p0': (3, 1577, 1)
    }
}
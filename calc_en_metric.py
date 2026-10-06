import numpy as np

def calc_en_metric(k, X, sX, Y, sY, sU):
    N = np.abs(X - Y) / (k * np.sqrt(sX**2 + sY**2 + sU**2))
    RMSE = np.sqrt(np.mean((X.flatten() - Y.flatten())**2))
    
    if len(X) > 1:
        R = np.corrcoef(X.flatten(), Y.flatten())
        R2 = R[1, 0]**2  # R[1, 0] gives the correlation coefficient
    else:
        R2 = np.nan
    
    return N, RMSE, R2

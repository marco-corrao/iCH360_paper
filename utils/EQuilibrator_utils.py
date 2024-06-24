import pandas as pd
import numpy as np
from component_contribution.linalg import LINALG
def compute_preprocessing_vectors(G:np.ndarray,
                                            S:np.ndarray,
                                            b:np.array):
    #1 Preprocessing vectors to estimate the mean
    #Linear Algebra magic
    GtS = G.T @ S
    inv_S, rank_S, P_R_S, P_N_S = LINALG.invert_project(S)
    inv_StG, rank_StG, P_R_StG, P_N_StG = LINALG.invert_project(GtS)
    inv_SWS, _, _, _ = LINALG.invert_project(S @ S.T)
    inv_GSWGS, _, _, _ = LINALG.invert_project(GtS @ GtS.T)

    # Compute pre_processing vectors (Equation 7)
    v_r = (P_R_S @ inv_S.T + P_N_S @ G @ inv_StG.T) @ b
    v_g = inv_StG.T @ b
    
    #===================================================
    #2. the C matrix for covariance estimation
    MSE_inf = 1e10
    # Equation 3
    Nr = S.shape[1]
    e_rc = (S.T @ P_R_S @ inv_S.T - np.eye(Nr)) @ b
    MSE_rc = (e_rc.T @ e_rc) / (Nr - rank_S)
    e_gc = (GtS.T @ inv_StG.T - np.eye(Nr)) @ b
    MSE_gc = (e_gc.T @ e_gc) / (Nr - rank_StG)
    C_rc = P_R_S @ inv_SWS @ P_R_S
    C_gc = P_N_S @ G @ inv_GSWGS @ G.T @ P_N_S
    C_inf = P_N_S @ G @ P_N_StG @ G.T @ P_N_S

    # Equation 12
    C_1 = MSE_rc * C_rc + MSE_gc * C_gc + MSE_inf * C_inf
    C_2 = MSE_gc * P_N_S @ G @ inv_GSWGS + MSE_inf * G @ P_N_StG
    C_3 = MSE_gc * inv_GSWGS + MSE_inf * P_N_StG
    C = np.block([[C_1, C_2], [C_2.T, C_3]])

    #===============================
    #The matrix Lq to compute the square root of the covariance
    #DrG0 covariance square root (Eq. 13/14)
    a_rc=np.sqrt(MSE_rc)
    a_gc=np.sqrt(MSE_gc)
    a_inf=np.sqrt(MSE_inf)
    L=np.block([[a_rc*inv_S@P_R_S,             np.zeros((inv_S.shape[0],inv_StG.shape[1]))],
                [a_gc*inv_StG @ G.T @ P_N_S,   a_gc*inv_StG],
                [a_inf*P_N_StG@G.T,            a_inf*P_N_StG]]).T


    #Compute Lq
    Lq=LINALG.qr_rank_deficient(L.T).T
    return v_r,v_g,C,Lq



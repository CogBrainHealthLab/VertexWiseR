"""Python function to run SLM with user-defined template data directory
   (as fetch_template_surface set to default in the original function)
   Modified from BrainStat toolbox https://brainstat.readthedocs.io/en/master/_modules/brainstat/stats/SLM.html#SLM
   © Copyright 2021, MICA Lab, CNG Lab Revision 1f3068fb
   VertexWiseR-related modifications include:
   -added 'brainstat.stats' to lines 37-39 to avoid relative import
   -added the "data_dir" argument to the fetch_template_surface instead of going with default
   -qc() was removed as it was causing issues with reticulate import
   -the fetch_parcellation and fetch_yeo_networks_metadata functions and Yeo parcellation data fetching were removed as not needed here
   -A new "inverse" option was included to allow surface to be treated as a predictor and the contrast variable as its dependent variable, so there are more checks and manually computed parameters whether the effect is fixed or mixed.
"""

# type: ignore
""" Standard Linear regression models. """
import warnings
from cmath import sqrt
from pathlib import Path
from pprint import pformat
from typing import Optional, Sequence, Tuple, Union

import numpy as np
import numpy.linalg as la
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy.stats import skew, kurtosis
from brainspace.mesh.mesh_elements import get_cells, get_points
from brainspace.vtk_interface.wrappers.data_object import BSPolyData
from nibabel.nifti1 import Nifti1Image

from brainstat._typing import ArrayLike
from brainstat.datasets import fetch_template_surface
from brainstat.datasets.base import fetch_template_surface
from brainstat.mesh.utils import _mask_edges, mesh_edges
from brainstat.stats.terms import FixedEffect, MixedEffect
from brainstat.stats.utils import apply_mask, undo_mask
from brainstat.stats._linear_model import _get_design_matrix, _check_constant_term, _model_univariate_fixed_effects, _model_univariate_mixed_effects, _compute_resls 
from brainstat.stats._t_test import _t_test

class SLM:
    """Core Class for running BrainStat linear models"""

    # Import class methods
    from brainstat.stats._linear_model import _linear_model
    from brainstat.stats._multiple_comparisons import _fdr, _random_field_theory
    from brainstat.stats._t_test import _t_test

    def __init__(
        self,
        model: Union[FixedEffect, MixedEffect],
        contrast: ArrayLike,
        surf: Optional[Union[str, dict, BSPolyData]] = None,
        mask: Optional[ArrayLike] = None,
        *,
        correction: Optional[Union[str, Sequence[str]]] = None,
        thetalim: float = 0.01,
        drlim: float = 0.1,
        two_tailed: bool = True,
        cluster_threshold: float = 0.001,
        data_dir: Optional[Union[str, Path]] = None,
        inverse: bool = False
    ) -> None:
        """Constructor for the SLM class.

        Parameters
        ----------
        model : brainstat.stats.terms.FixedEffect, brainstat.stats.terms.MixedEffect
            The linear model to be fitted of dimensions (observations, predictors).
            Note that, for volumetric input, BrainStat follows Fortran (MATLAB)
            convention for ordering voxels, i.e. the first dimension changes
            first.
        contrast : array-like
            Vector of contrasts in the observations.
        surf : str, dict, BSPolyData, Nifti1Image, optional
            A surface provided as either a dictionary with keys 'tri' for its
            faces (n-by-3 array) and 'coord' for its coordinates (3-by-n array),
            or as a BrainSpace BSPolyData object, a string containing a template
            name accepted by fetch_template_surface, or a Nifti1Image wherein 0
            denotes excluded voxels and any other value denotes included voxels,
            by default None.
        mask : array-like, optional
            A mask containing True for vertices to include in the analysis, by
            default None.
        correction : str, Sequence, optional
            String or sequence of strings. If it contains "rft" a random field
            theory multiple comparisons correction will be run. If it contains
            "fdr" a false discovery rate multiple comparisons correction will be
            run. Both may be provided. By default None.
        thetalim : float, optional
            Lower limit on variance coefficients in standard deviations, by
            default 0.01.
        drlim : float, optional
            Step of ratio of variance coefficients in standard deviations, by
            default 0.1.
        two_tailed : bool, optional
            Determines whether to return two-tailed or one-tailed p-values. Note
            that multivariate analyses can only be two-tailed, by default True.
        cluster_threshold : float, optional
            P-value threshold or statistic threshold for defining clusters in
            random field theory, by default 0.001.
        data_dir : str, pathlib.Path, optional
            Path to the location to store BrainStat data files, defaults to
            $HOME_DIR/brainstat_data.
        inverse : bool, optional
            Determines whether to test the surface data as contrast instead of 
            as dependent variable (vertex-wise value swapped with contrast in each 
            vertex-wise linear model, covariates remaining in place), by default 
            False.
        """
        # Input arguments.
        self.inverse = inverse
        self.model = model
        self.contrast = np.array(contrast)
        """ 
        VertexWiseR change: self.data_dir is defined earlier, so it can be used as 
        an argument in fetch_template_surface()
        """
        self.data_dir = data_dir
        if isinstance(surf, str):
            self.surf_name = surf
            self.surf = fetch_template_surface(self.surf_name, data_dir=self.data_dir)
        else:
            self.surf_name = None
            self.surf = surf

        self.mask = mask
        self.correction = [correction] if isinstance(correction, str) else correction
        self.niter = 1
        self.thetalim = thetalim
        self.drlim = drlim
        self.two_tailed = two_tailed
        self.cluster_threshold = cluster_threshold


        # Error check
        if self.surf is None:
            if self.correction is not None and "rft" in self.correction:
                raise ValueError("Random Field Theory corrections require a surface.")

        # We have to initialize fit parameters for our unit tests here.
        self._reset_fit_parameters()


    def fit(self, Y: np.ndarray) -> None:
        """Fits the SLM model

        Parameters
        ----------
        Y : numpy.array
            Input data (observation, vertex, variate)

        Raises
        ------
        ValueError
            An error will be thrown when multivariate data is provided and a
            one-tailed test is requested.
        """
        if Y.ndim < 2 or Y.ndim > 3:
            raise ValueError("Input data must be two or three dimensional.")
        elif Y.ndim == 3:
            if (not self.two_tailed) and Y.shape[2] > 1:
                raise NotImplementedError(
                    "One-tailed tests are not implemented for multivariate data."
                )
            if Y.shape[2] > 3 and "rft" in self.correction:
                raise NotImplementedError(
                    "Random Field Theory corrections are not implemented for more than three variates."
                )
            student_t_test = Y.shape[2] == 1
        else:
            student_t_test = True

        self._reset_fit_parameters()
        if self.mask is not None:
            Y_masked = apply_mask(Y, self.mask, axis=1)
        else:
            Y_masked = Y.copy()
            
        if not self.inverse:    
            self._linear_model(Y_masked)
            self._t_test()
        #inverse orientation
        else:
            #mimics linear_model but modifies Y and X to reverse contrast
            #and DV
            n_samples = Y_masked.shape[0]
            
            self.X, self.V = _get_design_matrix(self, n_samples)
            _check_constant_term(self.X)
            X_base, V_base = self.X.copy(), (None if self.V is None else self.V.copy())
            self.df = n_samples - la.matrix_rank(self.X)
            #read predictors and remove the contrast variable from among them
            contrast_idx = [j for j in range(self.X.shape[1]) if np.allclose(self.X[:, j], self.contrast.reshape(-1))]
            covariates = np.delete(self.X, contrast_idx, axis=1)
            n_vertices = Y_masked.shape[1]
            #make the contrast the DV
            DV = self.X[:, contrast_idx].reshape(-1, 1) #as 2D array
            #prepare covariates (always the same besides vertex IV)
            X_cov = np.concatenate([covariates, np.empty((covariates.shape[0], 1), dtype=covariates.dtype)], axis=1)
            #For every vertex, test covariates+vertex against DV
            t_stats = np.zeros(n_vertices)
            coefs_all = np.zeros((covariates.shape[1]+1, n_vertices))
            residuals_all=np.empty((Y_masked.shape[0], n_vertices))
            for v in range(n_vertices):
                #add vertex column as predictor
                X_cov[:, -1] = Y_masked[:, v]
                self.X = X_cov
                self.Y = DV
                self.V = V_base
                #run model manually for single-column Y matrix
                if isinstance(self.model, MixedEffect):
                    residuals,  self.V, self.coef, self.SSE, self.r, self.dr = _model_univariate_mixed_effects(self,DV)
                else:
                    residuals, self.V, self.coef, self.SSE = _model_univariate_fixed_effects(self,DV)
                
                coefs_all[:, v] = self.coef[:, 0]
                residuals_all[:, v] = residuals[:, 0]
                self.contrast = self.X[:,-1] #vertex variable
                self._t_test() #update self.t 
                t_stats[v] = self.t[0, 0]
            
            #restore base definitions to
            self.X = X_base
            self.V = V_base
            self.t = t_stats.reshape(1, -1) #needs to be 2D to match
            self.coef = coefs_all
            #computing resels, normally in _linear_model
            if self.surf is not None:  
                self.resl, mesh_connections = _compute_resls(self, residuals_all)
                key = list(mesh_connections.keys())[0]
                setattr(self, key, mesh_connections[key])
        
        if self.mask is not None:
            self._unmask()
        if self.correction is not None:
            self.multiple_comparison_corrections(student_t_test=True)

    def multiple_comparison_corrections(self, student_t_test: bool) -> None:
        """Performs multiple comparisons corrections. If a (one-sided) student-t
        test was run, then make it two-tailed if requested."""
        P1, Q1 = self._run_multiple_comparisons()

        if self.two_tailed and student_t_test:
            self.t = -self.t
            P2, Q2 = self._run_multiple_comparisons()
            self.t = -self.t
            self.P = _merge_rft(P1, P2)
            self.Q = _merge_fdr(Q1, Q2)
        else:
            self.P = P1
            # Make variable type consistent with two one tailed outputs.
            for field in ["peak", "clus", "clusid"]:
                if self.P is not None and field in self.P:
                    if field == "clusid":
                        self.P[field] = [self.P[field]]
                    else:
                        for field2 in self.P[field].keys():
                            self.P[field][field2] = [self.P[field][field2]]
            self.Q = Q1

        self._surfstat_to_brainstat_rft()


    def _run_multiple_comparisons(self) -> Tuple[Optional[dict], Optional[np.ndarray]]:
        """Runs the multiple comparisons tests and returns their outputs.

        Returns
        -------
        dict, None
            Results of random_field_theory. None if not requested.
        np.array, None
            Results of fdr. None if not requested.
        """
        P = None
        Q = None
        if "rft" in self.correction:
            P = {}
            P["pval"], P["peak"], P["clus"], P["clusid"] = self._random_field_theory()
        if "fdr" in self.correction:
            Q = self._fdr()
        return P, Q

    def _reset_fit_parameters(self) -> None:
        """Sets empty parameters before fitting. Prevents issues arising from
        using the same object to fit twice.
        """
        self.X = None
        self.t = None
        self.df = None
        self.SSE = None
        self.coef = None
        self.V = None
        self.k = None
        self.r = None
        self.dr = None
        self.resl = None
        self.c = None
        self.ef = None
        self.sd = None
        self.dfs = None
        self.P = None
        self.Q = None
        self.du = None

    def _surfstat_to_brainstat_rft(self) -> None:
        """Convert SurfStat RFT output to BrainStat RFT output.
        In short, we convert the peak/clus arrays to pandas dataframes and add Yeo networks.
        """
        if self.P is not None:
            if "peak" in self.P:
                for i in range(len(self.P["peak"]["t"])):
                    if self.P["peak"]["vertid"][i] is not None:
                        self.P["peak"]["vertid"][i] = self.P["peak"]["vertid"][
                            i
                        ].astype(int)
                    if self.P["peak"]["clusid"][i] is not None:
                        self.P["peak"]["clusid"][i] = self.P["peak"]["clusid"][
                            i
                        ].astype(int)

            if "clus" in self.P:
                for i in range(len(self.P["clus"]["clusid"])):
                    if self.P["clus"]["clusid"][i] is not None:
                        self.P["clus"]["clusid"][i] = self.P["clus"]["clusid"][
                            i
                        ].astype(int)

            none_squeeze = lambda x: np.squeeze(x) if x is not None else None
            for field in ["peak", "clus"]:
                P_tmp = []
                for i in range(len(self.P[field]["P"])):
                    tail_dict = {
                        key: none_squeeze(value[i])
                        for key, value in self.P[field].items()
                    }
                    if tail_dict["P"] is not None:
                        if tail_dict["P"].size == 1:
                            P_tmp.append(pd.DataFrame.from_dict([tail_dict]))
                        else:
                            P_tmp.append(pd.DataFrame.from_dict(tail_dict))
                            P_tmp[i].sort_values(by="P", ascending=True)
                    else:
                        P_tmp.append(pd.DataFrame(columns=tail_dict.keys()))
                self.P[field] = P_tmp

    def _unmask(self) -> None:
        """Changes all masked parameters to their input dimensions."""
        simple_unmask_parameters = ["t", "coef", "SSE", "r", "ef", "sd", "dfs"]
        for key in simple_unmask_parameters:
            attr = getattr(self, key)
            if attr is not None:
                setattr(self, key, undo_mask(attr, self.mask, axis=1))

        # slm.resl unmask
        if self.resl is not None:
            edges = mesh_edges(self.surf)
            _, idx = _mask_edges(edges, self.mask)
            self.resl = undo_mask(self.resl, idx, axis=0)

    def __str__(self) -> str:
        """Returns a string representation of the model."""
        return pformat(
            {key: value for key, value in self.__dict__.items() if not callable(value)}
        )

    """ Property specifications. """

    @property
    def surf(self) -> Union[BSPolyData, dict, Nifti1Image, None]:
        return self._surf

    @surf.setter
    def surf(self, value):
        self._surf = value
        if self.surf is not None:
            if isinstance(self.surf, BSPolyData):
                self.tri = np.array(get_cells(self.surf)) + 1
                self.coord = np.array(get_points(self.surf)).T
            elif isinstance(self.surf, Nifti1Image):
                self.lat = self.surf.get_fdata() != 0
            else:
                if "tri" in value:
                    self.tri = value["tri"]
                    self.coord = value["coord"]
                elif "lat" in value:
                    self.lat = value["lat"]
                    self.coord = value["coord"]

    @surf.deleter
    def surf(self):
        del self._surf

    @property
    def tri(self) -> np.ndarray:
        if hasattr(self, "_tri"):
            return self._tri
        else:
            return None

    @tri.setter
    def tri(self, value):
        if value is not None and np.any(value < 0):
            raise ValueError("Triangle indices must be non-negative.")
        self._tri = value

    @tri.deleter
    def tri(self):
        del self._tri

    @property
    def lat(self) -> np.ndarray:
        if hasattr(self, "_lat"):
            return self._lat
        else:
            return None

    @lat.setter
    def lat(self, value):
        self._lat = value

    @lat.deleter
    def lat(self):
        del self._lat



def _merge_rft(P1: dict, P2: dict) -> dict:
    """Merge two one-tailed outputs of the random_field_theory function.

    Parameters
    ----------
    P1 : dict
        Output dict of random_field_theory
    P2 : dict
        Output dict of random_field_theory

    Returns
    -------
    dict
        Two-tailed version of the inputs.
    """
    if P1 is None and P2 is None:
        return None

    P = {}
    for key1 in P1 or P2:
        P[key1] = {}
        if key1 == "clusid":
            P[key1] = [P1[key1], P2[key1]]
            continue
        for key2 in P1[key1]:
            if key1 == "pval":
                P[key1][key2] = _onetailed_to_twotailed(P1[key1][key2], P2[key1][key2])
            else:
                P[key1][key2] = [P1[key1][key2], P2[key1][key2]]
    return P


def _merge_fdr(Q1: Optional[ArrayLike], Q2: Optional[ArrayLike]) -> np.ndarray:
    """Merge two one-tailed outputs of the fdr function.

    Parameters
    ----------
    Q1 : array-like, None
        Q-values
    Q2 : array-like, None
        Q-values

    Returns
    -------
    np.ndarray
        Two-tailed FDR p-values
    """
    if Q1 is None and Q2 is None:
        return None
    return _onetailed_to_twotailed(Q1, Q2)


def _onetailed_to_twotailed(p1: ArrayLike, p2: ArrayLike) -> np.ndarray:
    """Converts two one-tailed tests to a two-tailed test"""
    if p1 is None and p2 is None:
        return None
    elif p1 is None:
        p1 = p2
    elif p2 is None:
        p2 = p1

    return np.minimum(np.minimum(p1, p2) * 2, 1)


def f_test(slm1: SLM, slm2: SLM) -> SLM:
    """F-statistics for comparing two uni- or multi-variate fixed effects models.

    Parameters
    ----------
    slm1 : brainstat.stats.SLM.SLM
        Standard linear model returned by the t_test function; see Notes for
        details.
    slm2 : brainstat.stats.SLM.SLM
        Standard linear model returned by the t_test function; see Notes for
        details.

    Returns
    -------
    brainstat.stats.SLM.SLM
        Standard linear model with f-test results included.

    """

    if slm1.r is not None or slm2.r is not None:
        warnings.warn("Mixed effects models not programmed yet.")

    slm = SLM(FixedEffect(1), FixedEffect(1))
    if slm1.df > slm2.df:
        X1 = slm1.X
        X2 = slm2.X
        df1 = slm1.df
        df2 = slm2.df
        SSE1 = slm1.SSE
        SSE2 = slm2.SSE
        for key in slm2.__dict__:
            setattr(slm, key, getattr(slm2, key))
    else:
        X1 = slm2.X
        X2 = slm1.X
        df1 = slm2.df
        df2 = slm1.df
        SSE1 = slm2.SSE
        SSE2 = slm1.SSE
        for key in slm1.__dict__:
            setattr(slm, key, getattr(slm1, key))

    r = X1 - np.dot(np.dot(X2, np.linalg.pinv(X2)), X1)
    d = np.sum(r.flatten() ** 2) / np.sum(X1.flatten() ** 2)

    if d > np.spacing(1):
        print("Models are not nested.")
        return

    slm.df = np.array([[df1 - df2, df2]])
    h = SSE1 - SSE2

    # if slm['coef'] is 3D and third dimension is 1, then squeeze it to 2D
    if np.ndim(slm.coef) == 3 and np.shape(slm.coef)[2] == 1:
        x1, x2, x3 = np.shape(slm.coef)
        slm.coef = slm.coef.reshape(x1, x2)

    if np.ndim(slm.coef) == 2:
        slm.k = np.array(1)
        slm.t = np.dot(h / (SSE2 + (SSE2 <= 0)) * (SSE2 > 0), df2 / (df1 - df2))
    elif np.ndim(slm.coef) > 2:
        k2, v = np.shape(SSE2)
        k = np.around((np.sqrt(1 + 8 * k2) - 1) / 2)
        slm.k = np.array(k)
        if k > 3:
            print("Roy's max root for k>3 not programmed yet.")
            return

        l = min(k, df1 - df2)
        slm.t = np.zeros((int(l), int(v)))

        if k == 2:
            det = SSE2[0, :] * SSE2[2, :] - SSE2[1, :] ** 2
            a11 = SSE2[2, :] * h[0, :] - SSE2[1, :] * h[1, :]
            a21 = SSE2[0, :] * h[1, :] - SSE2[1, :] * h[0, :]
            a12 = SSE2[2, :] * h[1, :] - SSE2[1, :] * h[2, :]
            a22 = SSE2[0, :] * h[2, :] - SSE2[1, :] * h[1, :]
            a0 = a11 * a22 - a12 * a21
            a1 = (a11 + a22) / 2
            s1 = np.array([sqrt(x) for x in (a1**2 - a0)]).real
            d = (df2 / (df1 - df2)) / (det + (det <= 0)) * (det > 0)
            slm.t[0, :] = (a1 + s1) * d
            if l == 2:
                slm.t[1, :] = (a1 - s1) * d
        if k == 3:
            det = (
                SSE2[0, :] * (SSE2[2, :] * SSE2[5, :] - SSE2[4, :] ** 2)
                - SSE2[5, :] * SSE2[1, :] ** 2
                + SSE2[3, :] * (SSE2[1, :] * SSE2[4, :] * 2 - SSE2[2, :] * SSE2[3, :])
            )
            m1 = SSE2[2, :] * SSE2[5, :] - SSE2[4, :] ** 2
            m3 = SSE2[0, :] * SSE2[5, :] - SSE2[3, :] ** 2
            m6 = SSE2[0, :] * SSE2[2, :] - SSE2[1, :] ** 2
            m2 = SSE2[3, :] * SSE2[4, :] - SSE2[1, :] * SSE2[5, :]
            m4 = SSE2[1, :] * SSE2[4, :] - SSE2[2, :] * SSE2[3, :]
            m5 = SSE2[1, :] * SSE2[3, :] - SSE2[0, :] * SSE2[4, :]
            a11 = m1 * h[0, :] + m2 * h[1, :] + m4 * h[3, :]
            a12 = m1 * h[1, :] + m2 * h[2, :] + m4 * h[4, :]
            a13 = m1 * h[3, :] + m2 * h[4, :] + m4 * h[5, :]
            a21 = m2 * h[0, :] + m3 * h[1, :] + m5 * h[3, :]
            a22 = m2 * h[1, :] + m3 * h[2, :] + m5 * h[4, :]
            a23 = m2 * h[3, :] + m3 * h[4, :] + m5 * h[5, :]
            a31 = m4 * h[0, :] + m5 * h[1, :] + m6 * h[3, :]
            a32 = m4 * h[1, :] + m5 * h[2, :] + m6 * h[4, :]
            a33 = m4 * h[3, :] + m5 * h[4, :] + m6 * h[5, :]
            a0 = (
                -a11 * (a22 * a33 - a23 * a32)
                + a12 * (a21 * a33 - a23 * a31)
                - a13 * (a21 * a32 - a22 * a31)
            )
            a1 = a22 * a33 - a23 * a32 + a11 * a33 - a13 * a31 + a11 * a22 - a12 * a21
            a2 = -(a11 + a22 + a33)
            q = a1 / 3 - a2**2 / 9
            r = (a1 * a2 - 3 * a0) / 6 - a2**3 / 27
            s1 = (r + [sqrt(x) for x in (q**3 + r**2)]) ** (1 / 3)
            z = np.zeros((3, v))
            z[0, :] = 2 * s1.real - a2 / 3
            z[1, :] = -s1.real - a2 / 3 + np.sqrt(3) * s1.imag
            z[2, :] = -s1.real - a2 / 3 - np.sqrt(3) * s1.imag

            if not np.count_nonzero(z) == 0:
                z.sort(axis=0)
                z = z[::-1]
            d = df2 / (df1 - df2) / (det + (det <= 0)) * (det > 0)

            for j in range(0, l):
                slm.t[j, :] = z[j, :] * d
    return slm

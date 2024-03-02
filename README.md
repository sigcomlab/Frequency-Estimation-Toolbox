# Frequency Estimation Toolbox
**Introduction**

The problem of estimating the amplitude, phase and frequency of multiple tones in additive white Gaussian noise (AWGN) has received significant attention for a number of years because of its relevance in various fields, including radar systems and wireless communications. In the last years, some solutions to this problem have been proposed. Those available in this toolbox, which are described in detail in [1,2,3,4] and are based on approximate maximum likelihood (ML) estimation, exhibit the following advantages

1. turn a complicated multidimensional problem (whose dimensionality is usually unknown a priori) into a sequence of lower dimensional subproblems. Consequently, their overall complexity is proportional to that required to solve each of such subproblems and is usually much lower than that of: a) parametric estimation methods, like the MUSIC [5], the ESPRIT [6]; b) non parametric spectral estimators, like the Capon method [7], the amplitude and phase estimation of a sinusoid (APES) [8], the iterative adaptive approach for amplitude and phase estimation (IAA-APES) [9].
   
2. perform better than independently estimating the tones associated with the largest peaks of the original periodogram. In fact, they allow to identify peaks that are initially masked by the leakage due to nearby stronger tones.
   
3. are able to estimate an unknown L in a simple fashion. In fact, this result can be achieved setting the initial value of this parameter to zero and applying a suitable test to establish whether, at each repetition of the first step, the largest peak detected in the periodogram of the last residual is significant [10] or whether, at each repetition of the second step, the energy of the new residual is large enough [11]. If one of these conditions is satisfied, the estimate of L is incremented by one and, then, the next step is carried out; otherwise, the estimation process is terminated. It is worth stressing that various estimation methods (e.g., the MUSIC and the ESPRIT estimators) require prior knowledge of L and that, in these cases, the use of some methods, like the generalized Akaike information criterion [12] or the minimum description length [13] is commonly proposed for the estimation of this parameter; however, the computational effort they require is not negligible.

**Content**

In this toolbox, the following frequency estimation algorithms based on approximate ML are collected:

1.	Single Frequency Estimation and Cancellation (SFEC) algorithm [1] –
  
2.	Complex Single Frequency Estimation and Cancellation (SFEC) algorithm [2] –
  
3.	BLA BLA BLA [3] –
   
4.	BLA BLA BLA [4] –
   
Each algorithm, implemented by a dedicated function, can be tested through its test main script (also included in the toolbox). All the test main scripts contain two examples employing the frequency estimation algorithm to search for the strongest tones in: a) synthetically generated signal, and; b) data from real world measurements.

**References**

[1] P. Di Viesti, A. Davoli, G. Guerzoni and G. M. Vitetta, "Novel Deterministic Detection and Estimation Algorithms for Colocated Multiple-Input Multiple-Output Radars," in IEEE Access, vol. 10, pp. 2216-2255, 2022, doi: 10.1109/ACCESS.2021.3139200.

[2] Pasquale Di Viesti , ALESSANDRO DAVOLI , Giorgio Guerzoni , et al. “Novel Methods for Approximate Maximum Likelihood Estimation of Multiple Superimposed Undamped Tones and Their Application to Radar systems”. TechRxiv. August 03, 2021. DOI: 10.36227/techrxiv.15054321.v2

[3] CITA ARTICOLO MICHELE

[4] CITA ARTICOLO MICHELE

[5] R. Schmidt, “Multiple emitter location and signal parameter estimation,” IEEE Trans. Antennas Propag., vol. 34, no. 3, pp. 276–280, Mar. 1986.

[6] R. Roy and T. Kailath, “ESPRIT-estimation of signal parameters via rotational invariance techniques,” IEEE Trans. Acoust. Speech Signal Process., vol. 37, no. 7, pp. 984–995, Jul. 1989.

[7] J. Capon, “High-Resolution Frequency-Wavenumber Spectrum Analysis,” Proc. IEEE, vol. 57, no. 8, pp. 1408–1418, Aug. 1969. 

[8] J. Li and P. Stoica, “An Adaptive Filtering Approach to Spectral Estimation and SAR Imaging,” IEEE Trans. Signal Process., vol. 44, no. 6, pp. 1469–1484, Jun. 1996. 

[9] T. Yardibi, J. Li, P. Stoica, M. Xue, and A. B. Baggeroer, “Source Localization and Sensing: A Nonparametric Iterative Adaptive Approach Based on Weighted Least Squares,” IEEE Trans. Aerosp. Electron. Syst., vol. 46, no. 1, pp. 425–443, Jan. 2010.

[10] P. Whittle, “The simultaneous estimation of a time series harmonic components and covariance structure,” p. 15, 1952.

[11] P. Gough, “A Fast Spectral Estimation Algorithm Based on the FFT,” IEEE Trans. Signal Process., vol. 42, no. 6, pp. 1317–1322, Jun. 1994. 

[12] J. Li and P. Stoica, “Efficient Mixed-Spectrum Estimation with Applications to Target Feature Extraction,” IEEE Trans. Signal Process., vol. 44, no. 2, pp. 281–295, Feb. 1996. 

[13] J. Rissanen, “A Universal Prior for Integers and Estimation by Minimum Description Length,” Ann. Stat., vol. 11, no. 2, pp. 416–431, Jun. 1983.

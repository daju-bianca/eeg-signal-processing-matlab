# eeg-signal-processing-matlab
Automated MATLAB &amp; EEGLAB pipeline for non-invasive diabetes risk estimation using EEG signal processing (FIR, ASR, ICA, FFT), Sample Entropy, PCA, and K-Means clustering.
# EEG Signal Processing and Unsupervised Clustering for Non-Invasive Diabetes Risk Assessment

This repository contains an automated batch processing pipeline implemented in **MATLAB** and **EEGLAB** for the preprocessing, spectral analysis, nonlinear complexity evaluation, and unsupervised clustering of electroencephalographic (EEG) signals. The objective of this study is to identify atypical neurophysiological responses to standardized nociceptive stimuli associated with metabolic and vascular alterations in diabetes mellitus.

---

## 1. Dataset and Experimental Protocol

EEG signals were recorded using a **Bitbrain EEG** system across **12 channels** (`Fp1`, `Fp2`, `F3`, `F4`, `C3`, `C4`, `P3`, `P4`, `O1`, `O2`, `Fz`, `Pz`) positioned according to the International 10-20 system, at a sampling rate of **256 Hz**. Raw recordings (`.csv`) were converted into native EEGLAB format (`.set`) across two standardized experimental protocols:

| Experimental Protocol | Stimulus Type | Subjects | Mean Duration / Subject | Target Neurophysiological Response |
| :--- | :--- | :---: | :---: | :--- |
| **Cold Pressor Test (CPT)** | Thermal stress (ice water, 0–4°C) | 41 | ~169 s (2 min 49 s) | Sympathetic nervous system activation, Theta and Beta band dynamics |
| **Mechanical Stimulation** | Sharp nociceptive pinprick stimulus | 36 | ~11.7 s | Acute startle motor reflex, Beta and Gamma band reactivity |

---

## 2. Methodology and Processing Pipeline

### 2.1. EEG Preprocessing (MATLAB & EEGLAB)
To remove biological and environmental artifacts while preserving signal phase, the following automated preprocessing chain was applied:
1. **Digital FIR Filtering:** High-pass filter at **1 Hz** (baseline drift removal) and low-pass filter at **40 Hz**.
2. **Notch Filter (50 Hz):** Power line interference suppression.
3. **Artifact Subspace Reconstruction (ASR):** Automated removal of transient high-amplitude spikes.
4. **Re-referencing (AVG):** Common average reference computation across all 12 electrodes.
5. **Independent Component Analysis (ICA):** Computed at the `STUDY` level for blind source separation. Scalp topographic maps were used to isolate and reject ocular artifacts (frontal EOG activity on `Fp1`/`Fp2`, e.g., `IC1`, `IC3`) and muscular or poor electrode contact noise (peripheral EMG activity, e.g., `IC11`), retaining valid cortical sources (e.g., `IC2`).

### 2.2. Frequency-Domain Spectral Analysis (FFT)
Single-sided amplitude spectra were computed via Fast Fourier Transform (FFT) to extract absolute and relative power across the 5 clinical EEG frequency bands:
* **Delta (0.5–4 Hz):** Associated in clinical literature with increased amplitude in diabetic patients.
* **Theta (4–8 Hz):** Marker of cognitive processing and thermal discomfort perception.
* **Alpha (8–13 Hz):** Relaxed wakefulness (characteristic spectral peak confirmed at ~8.5 Hz; typically attenuated in diabetic cohorts).
* **Beta (13–30 Hz):** Active concentration, alertness, and cortical stress (12–30 Hz).
* **Gamma (30–40 Hz):** Sensory integration and rapid cognitive processing.

### 2.3. Nonlinear Complexity Analysis (Sample Entropy)
To quantify time-series irregularity and cortical adaptability across each of the 12 channels, **Sample Entropy** was computed:

$$\text{SampEn}(m, r, N) = -\ln\left(\frac{A}{B}\right)$$

Implemented parameters:
* Window length: $N = 3000$ samples
* Embedding dimension: $m = 2$
* Similarity tolerance threshold: $r = 0.2 \times \text{SD}$ (standard deviation of the signal)

High entropy values indicate complex, adaptable cortical dynamics, whereas low values reflect increased regularity, synchronization, or pathological rigidity.

### 2.4. Dimensionality Reduction and Unsupervised Clustering (PCA + K-Means)
* **Principal Component Analysis (PCA):** Applied to the feature matrix (`Subjects × 12 SampEn channels`) to project high-dimensional entropy profiles onto the first two principal components (`PC1` and `PC2`).
* **K-Means Clustering:** Unsupervised partitioning based on Euclidean distance minimization.
* **Cluster Validation:**
  * **Elbow Method (WSS):** Analysis of the *Within-Cluster Sum of Squares* confirmed the optimal inflection point at **$k = 2$ clusters** across both datasets.
  * **Silhouette Score:** Used to quantify intra-cluster cohesion and inter-cluster separation.

---

## 3. Results and Technical Interpretation

| Experimental Protocol | Explained Variance (PCA) | Silhouette Score | Cluster 1 (Typical Response) | Cluster 2 (Atypical Response) | Technical & Neurophysiological Findings |
| :--- | :---: | :---: | :--- | :--- | :--- |
| **Thermal Stress (CPT)** | PC1: 72.6%<br>PC2: 7.6% | **0.755** | **35 / 41 subjects**<br>Mean SampEn: `0.5 – 0.7` | **6 / 41 subjects**<br>Mean SampEn: `> 0.8` | Cluster 1 exhibits controlled neural synchronization and normal adaptation to thermal stress. Within the 6 atypical cases in Cluster 2, signal inspection revealed 3 hardware-induced outliers and 3 cases of genuine cortical maladaptation. |
| **Mechanical Stimulus** | PC1: 50.2%<br>PC2: 14.9% | **0.510** | **25 / 36 subjects**<br>Mean SampEn: `0.6 – 0.8` | **11 / 36 subjects**<br>Mean SampEn: `0.4 – 0.5` | Cluster 1 demonstrates rapid recovery and high cortical flexibility following the acute startle reflex. Cluster 2 exhibits a sharp post-stimulus entropy collapse, indicating neurophysiological rigidity and impaired rapid self-regulation. |

---

## 4. Repository Structure

* **`scripts/`** — MATLAB scripts (`.m`) for batch preprocessing, FIR/Notch filtering, ICA decomposition, FFT spectral extraction, Sample Entropy computation, and K-Means clustering.
* **`results/`** — Exported spreadsheets (`.xlsx` / `.csv`) containing band power metrics, channel-wise Sample Entropy values, and cluster assignments.
* **`figures/`** — Generated plots (filter frequency responses, ICA topographic maps, single-sided FFT spectra, Elbow curve, and PCA + Silhouette clustering diagrams).
* **`docs/`** — Project documentation and technical presentation slides.

---

## 5. Requirements

* **MATLAB** (R2022b or newer)
* **EEGLAB Toolbox** (for `.set` data structures, ASR, and STUDY-level ICA)
* **Signal Processing Toolbox** (digital filter design and FFT analysis)
* **Statistics and Machine Learning Toolbox** (PCA, K-Means, and Silhouette evaluation)

---

## 6. Future Work

* Clinical correlation of the unsupervised cluster assignments with blood glucose and HbA1c measurements.
* Extension to *Multiscale Sample Entropy (MSE)* to evaluate signal complexity across multiple temporal scales.
* Implementation of supervised classifiers (SVM, Neural Networks) upon clinical labeling of the cohort.

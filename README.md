# FFT_NTT_C

## GOAL
    Using C++ to implement Variant FFT/IFFT and NTT/INTT 

<br>

## Description
    FFT / NTT algorithm is used to accelerate the polynomial multiplication operation.

<br>

## File category
- software
    - naive_polymulti.cpp
        - naive polynomial multiplication
        - time complexity : $O(n^2)$
    - FFT_org.cpp
        - fast fourier transform
        - only using DIT-FFT (Cooley-Tukey) to implement
        - with preprocessing (reverse bit)
    - FFT_GSCT.cpp
        - fast fourier transform
        - using DIF-FFT (Centleman-Sande) and DIT-FFT (Cooley-Tukey) to implement FFT and IFFT respectively
        - without preprocessing (reverse bit)
    - FFT.cpp
        - still working ...
    - NTT_org.cpp
        - number theoretic transform
        - only using DIT-FFT (Cooley-Tukey) to implement
        - with preprocessing (reverse bit)
    - NTT_GSCT.cpp
        - number theoretic transform
        - using DIF-FFT (Centleman-Sande) and DIT-FFT (Cooley-Tukey) to implement NTT and INTT respectively
        - without preprocessing (reverse bit)
    - NTT.cpp
        - still working ...



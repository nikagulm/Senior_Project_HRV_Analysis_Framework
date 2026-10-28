#ifndef FFT_H
#define FFT_H
#include <vector>
#include <complex>

class FFT {
public:
    //radix-2 Cooley-Tukey FFT (length must be a power of 2)
    //inverse = true gives the inverse FFT (divided by N)
    static void fft(std::vector<std::complex<double>>& data, bool inverse);

    //DFT of any length (power of 2: fft(), otherwise Bluestein's algorithm)
    //inverse = true gives the inverse DFT (divided by N)
    static void dft(std::vector<std::complex<double>>& data, bool inverse);
};

#endif
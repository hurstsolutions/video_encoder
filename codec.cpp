#include "codec.hpp"

const int LUMINANCE_QUANT_TABLE[BLOCK_MATRIX_SIZE] = {
    16, 11, 10, 16, 24, 40, 51, 61,
    12, 12, 14, 19, 26, 58, 60, 55,
    14, 13, 16, 24, 40, 57, 69, 56,
    14, 17, 22, 29, 51, 87, 80, 62,
    18, 22, 37, 56, 68, 109, 103, 77,
    24, 35, 55, 64, 81, 104, 113, 92,
    49, 64, 78, 87, 103, 121, 120, 101,
    72, 92, 95, 98, 112, 100, 103, 99
};

const int ZIG_ZAG_INDEX_ORDER[BLOCK_MATRIX_SIZE]{
    0, 1, 5, 6, 14, 15, 27, 28,
    2, 4, 7, 13, 16, 26, 29, 42,
    3, 8, 12, 17, 25, 30, 41, 43,
    9, 11, 18, 24, 31, 40, 44, 53,
    10, 19, 23, 32, 39, 45, 52, 54,
    20, 22, 33, 38, 46, 51, 55, 60,
    21, 34, 37, 47, 50, 56, 59, 61,
    35, 36, 48, 49, 57, 58, 62, 63
};

static std::array<double, DCT_N> make_normal_factors(){
    std::array<double, DCT_N> a{};
    for(int i = 0; i < DCT_N; ++i){
        a[i] = (i==0) ? sqrt(1.0/DCT_N) : sqrt(2.0/DCT_N);
    }
    return a;
}

static std::array<std::array<double, DCT_N>, DCT_N> make_cos_table(){
    std::array<std::array<double, DCT_N>, DCT_N> t{};
    for( int k = 0; k < DCT_N; ++k){
        for (int n = 0; n < DCT_N; ++n){
            t[k][n] = cos((2*n + 1) * k * M_PI / (2.0 * DCT_N));
        }
    }
    return t;
}

const std::array<double, 8> NORMALIZATION_FACTORS = make_normal_factors();
const std::array<std::array<double, 8>, 8> COSINE_TABLE = make_cos_table();
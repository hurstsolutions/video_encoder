#pragma once
#include <vector>

constexpr int BLOCK_MATRIX_SIZE = 64;

struct Image {
    int width;
    int height;
    std::vector<unsigned char> data;
};

struct EncodedImage {
    int width;
    int height;
    std::vector<short> data;
};

struct YCbCrImage{
    int width;
    int height;
    std::vector<unsigned char> y_data;
    std::vector<unsigned char> cb_data;
    std::vector<unsigned char> cr_data;
};

struct Block{
    std::vector<double> data;
    Block() : data(BLOCK_MATRIX_SIZE, 0.0){}
};

extern const int LUMINANCE_QUANT_TABLE[BLOCK_MATRIX_SIZE];
extern const int ZIG_ZAG_INDEX_ORDER[BLOCK_MATRIX_SIZE];
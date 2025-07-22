#pragma once
#define _USE_MATH_DEFINES

#include <vector>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>

constexpr int BLOCK_MATRIX_SIZE = 64;
constexpr int DCT_N = 8;
constexpr int BLOCK_SIZE = 8;

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

struct MotionVector {
    int x;
    int y;
};

struct PFrameBlockData{
    MotionVector motion_vector;
    Block comp_residual_coeffs;
};

extern const int LUMINANCE_QUANT_TABLE[BLOCK_MATRIX_SIZE];
extern const int ZIG_ZAG_INDEX_ORDER[BLOCK_MATRIX_SIZE];

extern const std::array<double, 8> NORMALIZATION_FACTORS;
extern const std::array<std::array<double, 8>, 8> COSINE_TABLE;

extern Image read_ppm(const std::string& filename);
extern YCbCrImage rgb_to_ycbcr(const Image& rgb_img);
extern Block get_block_from_frame(const YCbCrImage& frame, int top_left_pixel_x, int top_left_pixel_y);
extern Block get_predicted_block(const YCbCrImage& reference_frame, int current_pixel_x, int current_pixel_y, const MotionVector& motion_vector);
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>

#include "codec.hpp"



YCbCrImage rgb_to_ycbcr(const Image& rgb_img){
    YCbCrImage ycbcr_img;
    ycbcr_img.width = rgb_img.width;
    ycbcr_img.height = rgb_img.height;

    const size_t num_pixels = rgb_img.width * rgb_img.height;
    ycbcr_img.y_data.resize(num_pixels);
    ycbcr_img.cb_data.resize(num_pixels);
    ycbcr_img.cr_data.resize(num_pixels);

    for (size_t i=0; i < num_pixels; ++i){
       // std::cout << "i = " << i << std::endl;
        double r = rgb_img.data[i * 3 + 0];
        double g = rgb_img.data[i * 3 + 1];
        double b = rgb_img.data[i * 3 + 2];

        //Convert to YCbCr
        double y = 0.299 * r + 0.587 * g + 0.114 * b;
        double cb = -0.168736 * r - 0.331264 * g + 0.5 * b + 128;
        double cr = 0.5 * r - 0.418688 * g - 0.081312 * b + 128;

        //std::cout << "<R G B> <" << r << " " << g << " " << b << ">" << std::endl;
        //std::cout << "<Y Cb Cr> <" << y << " " << cb << " " << cr << ">" << std::endl; 

        ycbcr_img.y_data[i] = static_cast<unsigned char>(std::max(0.0, std::min(255.0, y)));
        ycbcr_img.cb_data[i] = static_cast<unsigned char>(std::max(0.0, std::min(255.0, cb)));
        ycbcr_img.cr_data[i] = static_cast<unsigned char>(std::max(0.0, std::min(255.0, cr)));
    } 
    
    return ycbcr_img;
}

Image read_ppm(const std::string& filename){
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()){
        std::cerr << "Error could not open file " << filename << std::endl;
        exit(1);
    }

    std::string magic_number;
    file >> magic_number;
    if (magic_number != "P6"){
        std::cerr << "Error: not a P6 PPM file." << std::endl;
        exit(1);
    }

    Image img;
    int max_value;
    file >> img.width >> img.height >> max_value;
    std::cout << "Width - " << img.width << std::endl;
    std::cout << "Height - " << img.height << std::endl;
    file.ignore(1, '\n');

    // total size of image, because there are three bytes per pixel(RGB)
    const size_t data_size = img.width * img.height * 3;

    std::cout << "Data Size - " << data_size << std::endl;
    img.data.resize(data_size);

    file.read(reinterpret_cast<char*>(img.data.data()), data_size);

    if(!file){
        std::cerr << "Error - failed to read pixel data from " << filename << std::endl;
        exit(1);
    }
    file.close();

    return img;
}

void perform_dct(Block& block){
    std::vector<double> temp(BLOCK_MATRIX_SIZE);
    //store the row results in this temp vector and then store the column results from the temp v
    //vector into the original block. Transforms the original block.

    // 1D DCT on rows
    for (int i = 0; i < DCT_N; ++i){ // Process each row
        for (int j = 0; j < DCT_N; ++j){ //j is the frequency coefficient
            double sum = 0.0;
            for (int k = 0; k < DCT_N; ++k){ // Process each pixel
                sum += block.data[i * DCT_N + k] * cos((2*k+1) * j * M_PI / (2.0 * DCT_N));
            }
            double c = (j==0) ? sqrt(1.0/DCT_N) : sqrt(2.0/DCT_N);
            temp[i*DCT_N+j] = c * sum;
        }
    }

    // 1D DCT on columns
    for (int j = 0; j < DCT_N; ++j){
        for (int i = 0; i< DCT_N; ++i){
            double sum = 0.0;
            for (int k = 0; k < DCT_N; ++k){
                sum += temp[k*DCT_N+j] * cos((2*k+1) * i * M_PI / (2.0*DCT_N));
            }
            double c = (i==0) ? sqrt(1.0/DCT_N) : sqrt(2.0 / DCT_N);
            block.data[i * DCT_N + j] = c * sum;
        }
    }
}

void quantize_block(Block& block){
    double new_value;
    for (int i = 0; i < BLOCK_MATRIX_SIZE; ++i){
        new_value = round(block.data[i] / LUMINANCE_QUANT_TABLE[i]);
        block.data[i] = new_value;
    }
}

std::vector<double> zig_zag_scan(Block& block){
    std::vector<double> scanned_vector(BLOCK_MATRIX_SIZE);
    for (int i = 0; i < BLOCK_MATRIX_SIZE; ++i){
        scanned_vector[i] = block.data[ZIG_ZAG_INDEX_ORDER[i]];
    }
    return scanned_vector;
}

void intra_frame_code(const std::string& input_filename){
    Image rgb_image = read_ppm(input_filename);
    std::cout << "Encoder starting..." << std::endl;
    std::cout << "Successfully read image " << input_filename << std::endl;
    std::cout << "Dimensions: " << rgb_image.width << "x" << rgb_image.height << std::endl;

    YCbCrImage ycbcr_image = rgb_to_ycbcr(rgb_image);
    std::cout << "Successfully converted RGB to YCbCr" << std::endl;

    std::cout << "Size of Luma (y) data: " << ycbcr_image.y_data.size() << " bytes." << std::endl;
    
    
    const int block_size = 8;
    int width_in_blocks = ycbcr_image.width / block_size;
    int height_in_blocks = ycbcr_image.height / block_size;

    std::vector<Block> y_dct_blocks;

    std::cout << "Processing Y Plane..." << std::endl;

    // Visit each block
    for(int by = 0; by < height_in_blocks; ++by){
        for (int bx = 0; bx < width_in_blocks; ++bx){
            Block current_block;

            // Visit each pixel in each block
            for(int y = 0; y < block_size; ++y){
                for (int x = 0; x < block_size; ++x){
                    int pixel_x = bx * block_size + x;
                    int pixel_y = by * block_size + y;
                    int pixel_index = pixel_y * ycbcr_image.width + pixel_x;

                    current_block.data[y*block_size+x] = static_cast<double>(ycbcr_image.y_data[pixel_index]) - 128.0;
                }
            }

            perform_dct(current_block);
            quantize_block(current_block);
            y_dct_blocks.push_back(current_block);
        }
    }

    std::cout << "Finished DCT and Quantization for " << y_dct_blocks.size() << " Y Blocks." << std::endl;

    if(!y_dct_blocks.empty()){
        std::cout << "First 8x8 Y DCT blocks - Quantized:" << std::endl;
        for(int y = 0; y < block_size; ++y){
            for(int x = 0; x < block_size; ++x){
                std::cout << std::setw(8) << static_cast<int>(y_dct_blocks[0].data[y * block_size + x]);
            }
            std::cout << std::endl;
        }
    } else{
        std::cout << "No blocks found." << std::endl;
    }

    //Collect all the Zig Zag Scan data into one vector
    std::vector<short> all_zig_zag_scans;
    std::vector<double> scan;
    const int y_dct_blocks_size = y_dct_blocks.size();
    for(int i = 0; i < y_dct_blocks_size; ++i){
        scan = zig_zag_scan(y_dct_blocks[i]);
        for (double coeff : scan){
            all_zig_zag_scans.push_back(static_cast<short>(coeff));
        }
    }

    std::ofstream output_file("output.bin", std::ios::binary);
    if (!output_file.is_open()){
        std::cerr << "Error could not open output file " << std::endl;
        exit(1);
    }
    output_file.write(reinterpret_cast<char*>(&ycbcr_image.width), sizeof(ycbcr_image.width));
    output_file.write(reinterpret_cast<char*>(&ycbcr_image.height), sizeof(ycbcr_image.height));
    output_file.write(reinterpret_cast<char*>(all_zig_zag_scans.data()), all_zig_zag_scans.size()*sizeof(short));

    output_file.close();

    std::cout << "Successfully wrote compressed data." << std::endl;
    std::cout << "File size: " << sizeof(ycbcr_image.width) + sizeof(ycbcr_image.height) + all_zig_zag_scans.size() * sizeof(short) << " bytes. " << std::endl;
    
}

int SAD_score(const YCbCrImage& reference_image, const YCbCrImage& current_image, int ref_x, int ref_y, int cur_x, int cur_y, size_t block_size){
    int sum = 0;
    int width = reference_image.width;
    for(int x = 0; x < block_size; ++x){
        for(int y = 0; y < block_size; ++y){
            int reference_pixel_index = (y + ref_y) * width + (x + ref_x);
            int current_pixel_index = (y + cur_y) * width + (x + cur_x);
            sum += std::abs(current_image.y_data[current_pixel_index] - reference_image.y_data[reference_pixel_index]); 
        }
    }
    return sum;
}

MotionVector find_best_match(int block_x, int block_y, int search_area, size_t block_size, const YCbCrImage& current_frame, const YCbCrImage& reference_frame){
    MotionVector current_best_match{0, 0};
    int current_sad_score = INT_MAX;
    // Convert the indicies of the blocks we got in to pixel cooridinates
    const int pixel_x = block_x * block_size;
    const int pixel_y = block_y * block_size;

    for(int search_x = -search_area; search_x <= search_area; ++search_x){
        //Out of bounds check
        const int reference_block_x = pixel_x + search_x;
        if(reference_block_x < 0 || reference_block_x + block_size > current_frame.width || reference_block_x + block_size > reference_frame.width){
            continue;
        }
        for(int search_y = search_area; search_y <= search_area; ++search_y){
            //Out of bounds check
            const int reference_block_y = pixel_y + search_y;
            if(reference_block_y < 0 || reference_block_y + block_size > current_frame.height || reference_block_y + block_size > reference_frame.height){
                continue;
            }
            int sad_score = SAD_score(reference_frame, current_frame, reference_block_x, reference_block_y, pixel_x, pixel_y, block_size);
            if (sad_score < current_sad_score){
                current_sad_score = sad_score;
                current_best_match.x = reference_block_x;
                current_best_match.y = reference_block_y;
            }
        }
    }
    return current_best_match;
}

int main(){
    // Intra Frame Coder
    // const std::string input_file_name = "input2.ppm";
    // intra_frame_code(input_file_name);

    //Inter Frame Coder
    const std::string input_file_1 = "input_frame_001.ppm";
    const std::string input_file_2 = "input_frame_002.ppm";

    Image frame_1_image = read_ppm(input_file_1);
    Image frame_2_image = read_ppm(input_file_2);

    YCbCrImage frame_1 = rgb_to_ycbcr(frame_1_image);
    YCbCrImage frame_2 = rgb_to_ycbcr(frame_2_image);

    const int block_size = 8;
    const int search_area = 16;
    int width_in_blocks = frame_2.width / block_size;
    int height_in_blocks = frame_2.height / block_size;

    std::vector<MotionVector> motion_vectors;
    for (int block_y = 0; block_y < height_in_blocks; ++block_y){
        for(int block_x = 0; block_x < width_in_blocks; ++block_x){
            MotionVector best_match = find_best_match(block_x, block_y, search_area, block_size, frame_2, frame_1);
            motion_vectors.push_back(best_match);
        }
    }

    // Print some motion vectors to check work;
    std::cout << "Printing first 24 motion vectors." << std::endl;

    for(int i = 0; i < 24; ++i){
        std::cout << std::setw(8) << "(" << motion_vectors[i].x << ", " << motion_vectors[i].y <<")";
        if(!(i % 8) && i){
            std::cout << std::endl;
        }
    }




    return 0;
}
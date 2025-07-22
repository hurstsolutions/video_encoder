#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>

#include "codec.hpp"



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
    
    
    int width_in_blocks = ycbcr_image.width / BLOCK_SIZE;
    int height_in_blocks = ycbcr_image.height / BLOCK_SIZE;

    std::vector<Block> y_dct_blocks;

    std::cout << "Processing Y Plane..." << std::endl;

    // Visit each block
    for(int by = 0; by < height_in_blocks; ++by){
        for (int bx = 0; bx < width_in_blocks; ++bx){
            Block current_block;

            // Visit each pixel in each block
            for(int y = 0; y < BLOCK_SIZE; ++y){
                for (int x = 0; x < BLOCK_SIZE; ++x){
                    int pixel_x = bx * BLOCK_SIZE + x;
                    int pixel_y = by * BLOCK_SIZE + y;
                    int pixel_index = pixel_y * ycbcr_image.width + pixel_x;

                    current_block.data[y*BLOCK_SIZE+x] = static_cast<double>(ycbcr_image.y_data[pixel_index]) - 128.0;
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
        for(int y = 0; y < BLOCK_SIZE; ++y){
            for(int x = 0; x < BLOCK_SIZE; ++x){
                std::cout << std::setw(8) << static_cast<int>(y_dct_blocks[0].data[y * BLOCK_SIZE + x]);
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

int SAD_score(const YCbCrImage& reference_image, const YCbCrImage& current_image, int ref_x, int ref_y, int cur_x, int cur_y){
    int sum = 0;
    int width = reference_image.width;
    for(int x = 0; x < BLOCK_SIZE; ++x){
        for(int y = 0; y < BLOCK_SIZE; ++y){
            int reference_pixel_index = (y + ref_y) * width + (x + ref_x);
            int current_pixel_index = (y + cur_y) * width + (x + cur_x);
            sum += std::abs(current_image.y_data[current_pixel_index] - reference_image.y_data[reference_pixel_index]); 
        }
    }
    return sum;
}

MotionVector find_best_match(int block_x, int block_y, int search_area, const YCbCrImage& current_frame, const YCbCrImage& reference_frame){
    MotionVector current_best_match{0, 0};
    int current_sad_score = INT_MAX;
    // Convert the indicies of the blocks we got in to pixel cooridinates
    const int pixel_x = block_x * BLOCK_SIZE;
    const int pixel_y = block_y * BLOCK_SIZE;

    for(int search_x = -search_area; search_x <= search_area; ++search_x){
        //Out of bounds check
        const int reference_block_x = pixel_x + search_x;
        if(reference_block_x < 0 || reference_block_x + BLOCK_SIZE > current_frame.width || reference_block_x + BLOCK_SIZE > reference_frame.width){
            continue;
        }
        for(int search_y = search_area; search_y <= search_area; ++search_y){
            //Out of bounds check
            const int reference_block_y = pixel_y + search_y;
            if(reference_block_y < 0 || reference_block_y + BLOCK_SIZE > current_frame.height || reference_block_y + BLOCK_SIZE > reference_frame.height){
                continue;
            }
            int sad_score = SAD_score(reference_frame, current_frame, reference_block_x, reference_block_y, pixel_x, pixel_y);
            if (sad_score < current_sad_score){
                current_sad_score = sad_score;
                current_best_match.x = reference_block_x;
                current_best_match.y = reference_block_y;
            }
        }
    }
    return current_best_match;
}




Block create_residual_block(const Block& current_block, const Block& predicted_block){
    Block residual_block;
    for(int i = 0; i < BLOCK_MATRIX_SIZE; ++i){
        residual_block.data[i] = current_block.data[i] - predicted_block.data[i];
    }
    return residual_block;
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

    const int search_area = 16;
    int width_in_blocks = frame_2.width / BLOCK_SIZE;
    int height_in_blocks = frame_2.height / BLOCK_SIZE;

    // std::vector<MotionVector> motion_vectors;
    std::vector<PFrameBlockData> p_frame_blocks;
    for (int block_y = 0; block_y < height_in_blocks; ++block_y){
        for(int block_x = 0; block_x < width_in_blocks; ++block_x){
            int pixel_x_of_actual_block = block_x * BLOCK_SIZE;
            int pixel_y_of_actual_block = block_y * BLOCK_SIZE; 
            Block actual_block = get_block_from_frame(frame_2, pixel_x_of_actual_block, pixel_y_of_actual_block);
            MotionVector best_match = find_best_match(block_x, block_y, search_area, frame_2, frame_1);
            // motion_vectors.push_back(best_match);

            Block predicted_block = get_predicted_block(frame_1, pixel_x_of_actual_block, pixel_y_of_actual_block, best_match); 
            Block residual_block = create_residual_block(actual_block, predicted_block);
            perform_dct(residual_block);
            quantize_block(residual_block);
            p_frame_blocks.push_back({best_match, residual_block});
        }
    }

    // Print some motion vectors to check work;
    // std::cout << "Printing first 24 motion vectors." << std::endl;

    // for(int i = 0; i < 24; ++i){
    //     std::cout << std::setw(8) << "(" << motion_vectors[i].x << ", " << motion_vectors[i].y <<")";
    //     if(!(i % 8) && i){
    //         std::cout << std::endl;
    //     }
    // }
    
    //Verification of Motion Compensation and residual building.
    std::cout << "First Block's motion vector - (" << p_frame_blocks[0].motion_vector.x << ", " << p_frame_blocks[0].motion_vector.y << ")" << std::endl;
    std::cout << "First 8 Quantized Coefficients: " << std::endl;
    for (int i = 0; i < 8; ++i){
        std::cout << std::setw(8) << p_frame_blocks[0].comp_residual_coeffs.data[i];
    }
    std::cout << std::endl;

    std::ofstream output_file ("output_frame.bin", std::ios::binary);
    if(!output_file.is_open()){
        std::cerr << "Error - could not open output file." << std::endl;
        exit(1);
    }
    output_file.write(reinterpret_cast<char*>(&frame_2.width), sizeof(frame_2.width));
    output_file.write(reinterpret_cast<char*>(&frame_2.width), sizeof(frame_2.width));
    for(PFrameBlockData& p_frame : p_frame_blocks){
        output_file.write(reinterpret_cast<char*>(&p_frame.motion_vector.x), sizeof(p_frame.motion_vector.x));
        output_file.write(reinterpret_cast<char*>(&p_frame.motion_vector.y), sizeof(p_frame.motion_vector.y));
        std::vector<short> short_coeffs(64);
        for(int i = 0; i < p_frame.comp_residual_coeffs.data.size(); ++i){
            short_coeffs[i] = static_cast<short>(p_frame.comp_residual_coeffs.data[i]);
        }
        output_file.write(reinterpret_cast<char*>(&short_coeffs), short_coeffs.size() * sizeof(short));
    }
    output_file.close();

    std::cout << "Successfully wrote p-frame data." << std::endl;
    


    return 0;
}
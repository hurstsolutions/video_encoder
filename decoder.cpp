#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>

#include "codec.hpp"

EncodedImage read_binary(const std::string& filename){
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()){
        std::cerr << "Error - could not open file. " << filename << std::endl;
        exit(1);
    }

    EncodedImage img;
    file.read(reinterpret_cast<char*>(&img.width), sizeof(img.width));
    file.read(reinterpret_cast<char*>(&img.height), sizeof(img.height));
    std::cout << "Width - " << img.width << std::endl;
    std::cout << "Height - " << img.height << std::endl;

    file.seekg(0, std::ios::end); // Go to the end of the file
    size_t file_size = file.tellg(); // Return at what byte the reader's at
    size_t coeff_data_size = file_size - sizeof(img.width) - sizeof(img.height); // calculate the size of the rest of the coefficient data
    file.seekg(sizeof(img.width) + sizeof(img.height), std::ios::beg); // go to the start of the coefficient data
    img.data.resize(coeff_data_size / sizeof(short)); // Prep the vector by initializing the size
    file.read(reinterpret_cast<char*>(img.data.data()), coeff_data_size);

    if(!file){
        std::cerr << "Error - failed to read coefficient data." << std::endl;
        exit(1);
    }

    file.close();

    std::cout << "Succesffuly read " << img.data.size() << " coefficients. " << std::endl;

    return img;
}

std::vector<Block> inverse_zig_zag(const std::vector<short>& zz_data){
    //Vector to return at the end of the function
    std::vector<Block> all_blocks;
    //Structures to help processing in the loops
    Block new_block;
    size_t num_blocks = zz_data.size() / BLOCK_MATRIX_SIZE;
    all_blocks.reserve(num_blocks);

    for(size_t b = 0; b < num_blocks; ++b){
        Block new_block;
        size_t block_start = b * BLOCK_MATRIX_SIZE;
        for (int j = 0; j < BLOCK_MATRIX_SIZE; ++j){
            int coeff = zz_data[block_start + j];
            int destination = ZIG_ZAG_INDEX_ORDER[j];
            new_block.data[destination] = static_cast<double>(coeff);
        }

        all_blocks.push_back(std::move(new_block));
    }
    return all_blocks;
}

void dequantize_block(Block& block){
    for(int i = 0; i < BLOCK_MATRIX_SIZE; ++i){
        block.data[i] = block.data[i] * LUMINANCE_QUANT_TABLE[i];
    }
}

void inverse_dct(Block& block){
    double temp[64];
    //Column Pass
    for(int col = 0; col < DCT_N ; ++col){
        for(int x = 0; x < DCT_N; ++x){
            double sum = 0.0;
            for(int u = 0; u < DCT_N; ++u){
                sum += NORMALIZATION_FACTORS[u] * block.data[u*DCT_N+col] * COSINE_TABLE[u][x];
            }
            temp[x*DCT_N+col] = sum;
        }
    }
    //Row Pass
    for (int row = 0; row < DCT_N; ++row){
        for(int y = 0; y < DCT_N; ++y){
            double sum = 0.0;
            for(int v = 0; v < DCT_N; ++v){
                sum += NORMALIZATION_FACTORS[v] * temp[row*DCT_N+v] * COSINE_TABLE[v][y];
            }
            block.data[row*DCT_N+y] = sum;
        }
    }
}

void write_ppm(const std::string& filename, int width, int height, const std::vector<unsigned char>& data){
    std::ofstream output_file(filename, std::ios::binary);
    if(!output_file.is_open()){
        std::cerr << "Error could not open output file. " << std::endl;
        exit(1);
    }
    output_file << "P6\n" << width << " " << height << "\n255\n";
    output_file.write(reinterpret_cast<const char*>(data.data()), data.size());
    std::cout << "Successfully wrote to " << filename << std::endl;
    output_file.close();
}

void intra_frame_decoder(const std::string& output_filename){

    EncodedImage encoded_image = read_binary(output_filename);
    std::cout << "Decoder Starting...." << std::endl;
    std::cout << "Successfully read encoded image " << output_filename << std::endl;
    int encoded_width = encoded_image.width;
    int encoded_height = encoded_image.height;
    std::cout << "Dimensions: " << encoded_width << "x" << encoded_height << std::endl;

    // Perform the zig zag scan in reverse to rebuild each block
    std::vector<Block> inversed_blocks = inverse_zig_zag(encoded_image.data);
    //Dequantize the blocks
    for (Block& block : inversed_blocks){
        dequantize_block(block);
    }
    // Print the first block to verify
    if(!inversed_blocks.empty()){
        std::cout << "First 8x8 Y DCT blocks - DeQuantized:" << std::endl;
        for(int y = 0; y < BLOCK_SIZE; ++y){
            for(int x = 0; x < BLOCK_SIZE; ++x){
                std::cout << std::setw(8) << static_cast<int>(inversed_blocks[0].data[y * BLOCK_SIZE + x]);
            }
            std::cout << std::endl;
        }
    } else{
        std::cout << "No blocks found." << std::endl;
    }

    std::cout << "De-Quantized " << inversed_blocks.size() << " blocks" << std::endl;

    // Inverse DCT and Shift/Clamp values
    Image final_image;
    final_image.width = encoded_width;
    final_image.height = encoded_height;
    final_image.data.resize(encoded_width * encoded_height * 3);
    int width_in_blocks = encoded_width/ BLOCK_SIZE;
    for(size_t i = 0; i < inversed_blocks.size(); ++i){ // Go through each block
        int block_y = i / width_in_blocks;
        int block_x = i % width_in_blocks;
        Block& current_block = inversed_blocks[i];

        inverse_dct(current_block);

        for(int y=0; y < BLOCK_SIZE; ++y){ // Go through each pixel in each block
            for (int x=0; x < BLOCK_SIZE; ++x){
                int pixel_x = block_x * BLOCK_SIZE + x;
                int pixel_y = block_y * BLOCK_SIZE + y;

                double val = current_block.data[y * BLOCK_SIZE + x] + 128.0;
                unsigned char clamped_val = static_cast<unsigned char>(std::clamp(val, 0.0, 255.0));

                //Flip 180 and then fill the final image.
                int rot_x = (encoded_width - 1) - pixel_x;
                int rot_y = (encoded_height -1) - pixel_y;
                size_t final_index = (rot_y * encoded_width + rot_x) * 3;
                final_image.data[final_index + 0] = clamped_val;
                final_image.data[final_index + 1] = clamped_val;
                final_image.data[final_index + 2] = clamped_val;
            }
        }
    }

    write_ppm(".\\Images\\output.ppm", encoded_width, encoded_height, final_image.data);
}

void add_blocks(Block& destination, Block& block_1, Block& block_2){
    
    if(block_1.data.size() != block_2.data.size()){
        std::cerr << "Blocks must be the same size to add together." << std::endl;
        exit(1);
    }
    for(int i = 0; i < block_1.data.size(); ++i){
        destination.data[i] = block_1.data[i] + block_2.data[i];
    }

}

int main(){
    const std::string input_file_1 = ".\\Images\\input_frame_001.ppm";
    Image frame_1_image = read_ppm(input_file_1);
    YCbCrImage reference_frame = rgb_to_ycbcr(frame_1_image);
    std::cout << "Loaded Reference Frame." << std::endl;

    const std::string binary_file_name = ".\\Images\\output_frame.bin";
    std::ifstream binary_file(binary_file_name, std::ios::binary);
    if(!binary_file.is_open()){
        std::cerr << "Error - could not open file. " << std::endl;
        exit(1);
    }

    int p_frame_width;
    int p_frame_height;
    binary_file.read(reinterpret_cast<char*>(&p_frame_width), sizeof(p_frame_width));
    binary_file.read(reinterpret_cast<char*>(&p_frame_height), sizeof(p_frame_height));
    std::cout << "Width - " << p_frame_width << std::endl;
    std::cout << "Height - " << p_frame_height << std::endl;
    
    const int size_of_header = sizeof(p_frame_height) + sizeof(p_frame_width);
    const int increment = (sizeof(int) * 2) + (sizeof(short) * BLOCK_MATRIX_SIZE);
    
    binary_file.seekg(0, std::ios::end); //Go to the end of the file
    size_t file_size = binary_file.tellg(); // Return where the reader's at i.e. the total size of the file
    size_t p_frame_block_data_size = file_size - size_of_header; // Calculate what we haven't read yet.
    binary_file.seekg(size_of_header, std::ios::beg); // jump to start of p_frame block data
    std::vector<PFrameBlockData> p_frame_vector;
    PFrameBlockData p_frame;
    std::vector<short> short_coeffs(BLOCK_MATRIX_SIZE);
    while(binary_file.read(reinterpret_cast<char*>(&p_frame.motion_vector.x), sizeof(int)) &&
            binary_file.read(reinterpret_cast<char*>(&p_frame.motion_vector.y), sizeof(int)) &&
            binary_file.read(reinterpret_cast<char*>(short_coeffs.data()), (sizeof(short)*BLOCK_MATRIX_SIZE)))
    {
        for(int i = 0; i<BLOCK_MATRIX_SIZE; ++i){
            p_frame.comp_residual_coeffs.data[i] = static_cast<double>(short_coeffs[i]);
        }
        p_frame_vector.push_back(p_frame);
    }


    binary_file.close();

    std::cout << "Successfully read - " << p_frame_vector.size() << " p_frames. " << std::endl;

    Image final_image;
    final_image.width = p_frame_width;
    final_image.height = p_frame_height;
    final_image.data.resize(final_image.width * final_image.height * 3);
    const int width_in_blocks = p_frame_width/BLOCK_SIZE;
    for(int p_index = 0; p_index < p_frame_vector.size(); ++p_index){
        PFrameBlockData& p_frame = p_frame_vector[p_index];
        //Decode the residual coefficients
        dequantize_block(p_frame.comp_residual_coeffs);
        inverse_dct(p_frame.comp_residual_coeffs);

        Block& decoded_residual_block = p_frame.comp_residual_coeffs;
        MotionVector& mv = p_frame.motion_vector;
        //Calculate top left pixel coordinates for each block 
        const int top_left_pixel_x = (p_index % width_in_blocks) * BLOCK_SIZE; //Find which Block and then * 8 for the pixel cooridinate
        const int top_left_pixel_y = (p_index / width_in_blocks) * BLOCK_SIZE;
        Block predicted_block = get_predicted_block(reference_frame, top_left_pixel_x, top_left_pixel_y, mv);
        
        //Add the Decoded Residual Block and the Predicted Block together. Then assemble the image.
        Block final_block;
        add_blocks(final_block, decoded_residual_block, predicted_block);
        for(int y=0; y < BLOCK_SIZE; ++y){ // Go through each pixel in each block
            for (int x=0; x < BLOCK_SIZE; ++x){
                int pixel_x = top_left_pixel_x + x;
                int pixel_y = top_left_pixel_y + y;

                double val = final_block.data[y * BLOCK_SIZE + x] + 128.0;
                unsigned char clamped_val = static_cast<unsigned char>(std::clamp(val, 0.0, 255.0));


                size_t final_index = (pixel_y * final_image.width + pixel_x) * 3;
                final_image.data[final_index + 0] = clamped_val;
                final_image.data[final_index + 1] = clamped_val;
                final_image.data[final_index + 2] = clamped_val;
            }
        }
    }

    write_ppm(".\\Images\\output_frame2.ppm", final_image.width, final_image.height, final_image.data);

    

    return 0;
}
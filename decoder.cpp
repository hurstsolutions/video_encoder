#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

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
            int coeff = zz_data[block_start + ZIG_ZAG_INDEX_ORDER[j]];
            int destination = ZIG_ZAG_INDEX_ORDER[j];
            new_block.data[j] = static_cast<double>(coeff);
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


int main(){
    const std::string output_filename = "output.bin";
    EncodedImage encoded_image = read_binary(output_filename);
    std::cout << "Decoder Starting...." << std::endl;
    std::cout << "Successfully read encoded image " << output_filename << std::endl;
    int encoded_width = encoded_image.width;
    int encoded_height = encoded_image.height;
    std::cout << "Dimensions: " << encoded_width << "x" << encoded_height << std::endl;

    // Perform the zig zag scan in reverse to rebuild each block
    const int block_size = 8;
    std::vector<Block> inversed_blocks = inverse_zig_zag(encoded_image.data);
    //Dequantize the blocks
    for (Block& block : inversed_blocks){
        dequantize_block(block);
    }
    // Print the first block to verify
    if(!inversed_blocks.empty()){
        std::cout << "First 8x8 Y DCT blocks - DeQuantized:" << std::endl;
        for(int y = 0; y < block_size; ++y){
            for(int x = 0; x < block_size; ++x){
                std::cout << std::setw(8) << static_cast<int>(inversed_blocks[0].data[y * block_size + x]);
            }
            std::cout << std::endl;
        }
    } else{
        std::cout << "No blocks found." << std::endl;
    }

    std::cout << "De-Quantized " << inversed_blocks.size() << " blocks" << std::endl;


    return 1;
}
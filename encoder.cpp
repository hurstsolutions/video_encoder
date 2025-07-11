#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <cmath>

struct Image {
    int width;
    int height;
    std::vector<unsigned char> data;
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
    Block() : data(64, 0.0){}
};

YCbCrImage rgb_to_ycbcr(const Image& rgb_img){
    YCbCrImage ycbcr_img;
    ycbcr_img.width = rgb_img.width;
    ycbcr_img.height = rgb_img.height;

    const size_t num_pixels = rgb_img.width * rgb_img.height;
    ycbcr_img.y_data.resize(num_pixels);
    ycbcr_img.cb_data.resize(num_pixels);
    ycbcr_img.cr_data.reserve(num_pixels);

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

    return img;
}

void perform_dct(Block& block){
    const int N = 8;
    std::vector<double> temp(64);
    //store the row results in this temp vector and then store the column results from the temp v
    //vector into the original block. Transforms the original block.

    // 1D DCT on rows
    for (int i = 0; i < N; ++i){ // Process each row
        for (int j = 0; j < N; ++j){ //j is the frequency coefficient
            double sum = 0.0;
            for (int k = 0; k < N; ++k){ // Process each pixel
                sum += block.data[i * N + k] * cos((2*k+1) * j * M_PI / (2.0 * N));
            }
            double c = (j==0) ? sqrt(1.0/N) : sqrt(2.0/N);
            temp[i*N+j] = c * sum;
        }
    }

    // 1D DCT on columns
    for (int j = 0; j < N; ++j){
        for (int i = 0; i< N; ++i){
            double sum = 0.0;
            for (int k = 0; k < N; ++k){
                sum += temp[k*N+j] * cos((2*k+1) * i * M_PI / (2.0*N));
            }
            double c = (i==0) ? sqrt(1.0/N) : sqrt(2.0 / N);
            block.data[i * N + j] = c * sum;
        }
    }
}


int main(){
    const std::string input_filename = "input2.ppm";
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
        for (int bx = 0; bx < height_in_blocks; ++bx){
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

            y_dct_blocks.push_back(current_block);
        }
    }

    std::cout << "Finished DCT for " << y_dct_blocks.size() << " Y Blocks." << std::endl;

    if(!y_dct_blocks.empty()){
        std::cout << "First 8x8 Y DCT blocks (top-left is DC coeff):" << std::endl;
        for(int y = 0; y < block_size; ++y){
            for(int x = 0; x < block_size; ++x){
                std::cout << std::fixed << std::setprecision(1) << std::setw(8) << y_dct_blocks[0].data[y * block_size + x];
            }
            std::cout << std::endl;
        }
    } else{
        std::cout << "No blocks found." << std::endl;
    }
    
    
    return 0;
}
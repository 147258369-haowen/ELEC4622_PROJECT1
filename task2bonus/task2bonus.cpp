// Task_1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include "io_bmp.h"
#include "image_comps.h"
extern float h1D[3];
extern float h2D[3];
int main(int argc, char* argv[]) {
    STATE state = CHECK_INPUT;
    int filterChoose = 0;
    bmp_in in;
    bmp_out out;
    my_image_comp* input_comps = nullptr;
    my_image_comp* output_comps = nullptr;
    io_byte* line = nullptr;
    float** matrix = nullptr;
    my_image_comp* temp_comps = nullptr;
    my_image_comp* dog_x = nullptr;
    my_image_comp* dog_y = nullptr;
    ImageParam imageParam;
    float sigma = 0;
    while (true) {
        switch (state) {
        case CHECK_INPUT:
            CheckInput(argc, argv, &sigma, &filterChoose, &imageParam);
            state = filterChoose ? LOAD_PICTURE : MOVING_AVERAGE_CHECK;
            break;
        case LOAD_PICTURE:
            //printf("%f\r\n", FilterInit(matrix, HIGHT, WIDTH));//归一化
            imageParam.GaussianDimension = GaussianWindowDimensionChoose(sigma);
            LoadImage(&in, &input_comps, &output_comps, &line, &imageParam, &filterChoose, argv);
            state = filterChoose ? DOG : MOVING_AVERAGE;
            break;
        case GAUSSIAN_FILTER:
            printf("Gaussan dimension: %d\r\n", imageParam.GaussianDimension);
            matrix = allocateMatrix(imageParam.GaussianDimension);
            LoadGaussianValue(matrix, sigma, imageParam.GaussianDimension);// load the gaussian value to the matrix
            for (int n = 0; n < imageParam.num_comp; n++)
                input_comps[n].perform_boundary_extension();
            for (int n = 0; n < imageParam.num_comp; n++)
                output_comps[n].init(imageParam.height, imageParam.width, 0); // Don't need a border for output
            temp_comps = new  my_image_comp[imageParam.num_comp];
            for (int n = 0; n < imageParam.num_comp; n++)
                temp_comps[n].init(imageParam.height, imageParam.width, (imageParam.GaussianDimension - 1) / 2); // Don't need a border for output
            for (int n = 0; n < imageParam.num_comp; n++)
                horizontal(input_comps + n, temp_comps + n, matrix, imageParam.GaussianDimension);
            for (int n = 0; n < imageParam.num_comp; n++)
                vertical(temp_comps + n, output_comps + n, matrix, imageParam.GaussianDimension);
            /*for (int n = 0; n < imageParam.num_comp; n++)
                apply_filter_modified(input_comps + n, output_comps + n, matrix, imageParam.GaussianDimension);*/
            state = (imageParam.gradientFlag) ? IMAGE_GRADIENT : OUTPUT_PICTURE;
            delete[] temp_comps;
            break;
        case MOVING_AVERAGE_CHECK:
            VarianceLoopCheck(sigma, &imageParam);
            state = LOAD_PICTURE;
            break;
        case MOVING_AVERAGE:
            matrix = allocateMatrix(imageParam.MV_Dimension);
            MovingAverageSetValue(matrix, imageParam.MV_Dimension);
            for (int n = 0; n < imageParam.num_comp; n++)
                input_comps[n].perform_boundary_extension();
            for (int n = 0; n < imageParam.num_comp; n++)
                output_comps[n].init(imageParam.height, imageParam.width, 0); // Don't need a border for output
            temp_comps = new  my_image_comp[imageParam.num_comp];
            for (int n = 0; n < imageParam.num_comp; n++)
                temp_comps[n].init(imageParam.height, imageParam.width, (imageParam.MV_Dimension - 1) / 2); // Don't need a border for output
            for (int n = 0; n < imageParam.num_comp; n++)
                horizontal(input_comps + n, temp_comps + n, matrix, imageParam.MV_Dimension);
            for (int n = 0; n < imageParam.num_comp; n++)
                vertical(temp_comps + n, output_comps + n, matrix, imageParam.MV_Dimension);
            state = (imageParam.gradientFlag) ? IMAGE_GRADIENT : OUTPUT_PICTURE;
            delete[] temp_comps;
            break;
        case IMAGE_GRADIENT:
            printf("Gradient\r\n");
            temp_comps = new  my_image_comp[imageParam.num_comp];
            for (int n = 0; n < imageParam.num_comp; n++)
                temp_comps[n].init(imageParam.height, imageParam.width, 1); // Don't need a border for output
            for (int n = 0; n < imageParam.num_comp; n++)
                temp_comps[n].GrradientHorizontalFilter(output_comps + n, 3, imageParam.alpha,h1D);
            for (int n = 0; n < imageParam.num_comp; n++)
                output_comps[n].GrradientverticalFilter(temp_comps + n, 3, imageParam.alpha,h2D);
            state = OUTPUT_PICTURE;
            break;
        case DOG:
            printf("DOG\n");
            for (int n = 0; n < imageParam.num_comp; n++)
                output_comps[n].init(imageParam.height, imageParam.width, 0); // Don't need a border for output
            imageParam.dog_x = new float[imageParam.GaussianDimension];
            imageParam.dog_y = new float[imageParam.GaussianDimension];
            LoadDerivativeGaussianValue(imageParam.dog_x, sigma, imageParam.GaussianDimension,1); //1 for x
            LoadDerivativeGaussianValue(imageParam.dog_y, sigma, imageParam.GaussianDimension, 0); //0 for y
            printf("\r\n");
            dog_x = new  my_image_comp[imageParam.num_comp];
            for (int n = 0; n < imageParam.num_comp; n++)
                dog_x[n].init(imageParam.height, imageParam.width, (imageParam.GaussianDimension - 1) / 2);
            dog_y = new  my_image_comp[imageParam.num_comp];
            for (int n = 0; n < imageParam.num_comp; n++)
                dog_y[n].init(imageParam.height, imageParam.width, (imageParam.GaussianDimension - 1) / 2);

            for (int n = 0; n < imageParam.num_comp; n++)
                dog_x[n].GrradientHorizontalFilter(input_comps + n, imageParam.GaussianDimension, imageParam.alpha, imageParam.dog_x);
            for (int n = 0; n < imageParam.num_comp; n++)
                dog_y[n].GrradientverticalFilter(input_comps + n, imageParam.GaussianDimension, imageParam.alpha, imageParam.dog_y);

            for (int i = 0; i < 3; i++) {
                for (int y = 0; y < imageParam.height; y++) {
                    for (int x = 0; x < imageParam.width; x++) {
                        float* ip = dog_x[i].buf + y * dog_x[0].stride + x;
                        float* ip_2 = dog_y[i].buf + y * dog_y[0].stride + x;
                        float* op = output_comps[i].buf + y * output_comps[0].stride + x;
                        float magnitude = imageParam.alpha *sqrt((*ip) * (*ip) + (*ip_2) * (*ip_2));
                        CLAMP_TO_BYTE(magnitude);
                        *op = magnitude;
                       
                    }
                }
                
            }

            state = OUTPUT_PICTURE;
            break;
        case OUTPUT_PICTURE:
            OutputImage(&out, input_comps, &output_comps, &line, &imageParam, argv);
            bmp_in__close(&in);
            delete[] line;
            delete[] input_comps;
            delete[] output_comps;

            return 1;
        default:
            break;
        }
    }
}
// Task_1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include "io_bmp.h"
#include "image_comps.h"
extern float h1[5][5];
extern float h2[9][9];
extern float h3[9][9];
int main(int argc, char* argv[]) {
    STATE state = CHECK_INPUT;
    int filterChoose = 0;
    bmp_in in;
    bmp_out out;
    my_image_comp* input_comps = nullptr;
    my_image_comp* output_comps = nullptr;
    io_byte* line = nullptr;
    ImageParam imageParam;
    float sigma = 0;
    float** matrix = allocateMatrix(WIDTH);
    while (true) {
        switch (state) {
        case CHECK_INPUT:
            CheckInput(argc, argv,&sigma,&filterChoose);
            state = filterChoose ? LOAD_PICTURE : MOVING_AVERAGE_CHECK;
            break;
        case LOAD_PICTURE:
            //printf("%f\r\n", FilterInit(matrix, HIGHT, WIDTH));//归一化
            LoadImage(&in, &input_comps, &output_comps,&line, &imageParam, &filterChoose,argv);
            state = filterChoose ? GAUSSIAN_FILTER : MOVING_AVERAGE;
            break;
        case GAUSSIAN_FILTER:
            LoadGaussianValue(matrix, sigma);// load the gaussian value to the matrix
            FilterInit(matrix, HIGHT, WIDTH);//normalize
            for (int n = 0; n < imageParam.num_comp; n++)
                input_comps[n].perform_boundary_extension();
            for (int n = 0; n < imageParam.num_comp; n++)
                output_comps[n].init(imageParam.height, imageParam.width, 0); // Don't need a border for output
            for (int n = 0; n < imageParam.num_comp; n++)
                apply_filter_modified(input_comps + n, output_comps + n, matrix, WIDTH);
            state = OUTPUT_PICTURE;
            break;
        case MOVING_AVERAGE_CHECK:
            VarianceLoopCheck(sigma, &imageParam);
            state = LOAD_PICTURE;
            break;
        case MOVING_AVERAGE:       
            FreeMatrix(matrix,WIDTH);
            matrix = allocateMatrix(imageParam.MV_Dimension);
            MovingAverageSetValue(matrix, imageParam.MV_Dimension);
            for (int n = 0; n < imageParam.num_comp; n++)
                input_comps[n].perform_boundary_extension();
            for (int n = 0; n < imageParam.num_comp; n++)
                output_comps[n].init(imageParam.height, imageParam.width, 0); // Don't need a border for output
            for (int n = 0; n < imageParam.num_comp; n++)
                apply_filter_modified(input_comps + n, output_comps + n, matrix, imageParam.MV_Dimension);
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
// Task_1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include "io_bmp.h"
#include "image_comps.h"
extern float laplacianKernel[3][3];
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
    my_image_comp* temp_comps_2 = nullptr;
    my_image_comp* temp_comps_out = nullptr;
    ImageParam imageParam;
    float sigma = 0;
    //float** lapkernel = allocateMatrix(3);
    float* lapkernel = new float[3*3];
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
            state = filterChoose ? GAUSSIAN_FILTER : MOVING_AVERAGE;
            break;
        case GAUSSIAN_FILTER:
            printf("Gaussan dimension: %d\r\n", imageParam.GaussianDimension);
            matrix = allocateMatrix(imageParam.GaussianDimension);
            LoadGaussianValue(matrix, sigma, imageParam.GaussianDimension);// load the gaussian value to the matrix
            for (int n = 0; n < imageParam.num_comp; n++)
                input_comps[n].perform_boundary_extension();
            for (int n = 0; n < imageParam.num_comp; n++)
                output_comps[n].init(imageParam.height, imageParam.width, 1); // Don't need a border for output
            temp_comps = new  my_image_comp[imageParam.num_comp];
            for (int n = 0; n < imageParam.num_comp; n++)
                temp_comps[n].init(imageParam.height, imageParam.width, (imageParam.GaussianDimension - 1) / 2); // Don't need a border for output
            for (int n = 0; n < imageParam.num_comp; n++)
                horizontal(input_comps + n, temp_comps + n, matrix, imageParam.GaussianDimension, 0);
            for (int n = 0; n < imageParam.num_comp; n++)
                vertical(temp_comps + n, output_comps + n, matrix, imageParam.GaussianDimension, 0);
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
                output_comps[n].init(imageParam.height, imageParam.width, 1); // Don't need a border for output
            temp_comps = new  my_image_comp[imageParam.num_comp];
            for (int n = 0; n < imageParam.num_comp; n++)
                temp_comps[n].init(imageParam.height, imageParam.width, (imageParam.MV_Dimension - 1) / 2); // Don't need a border for output
            for (int n = 0; n < imageParam.num_comp; n++)
                horizontal(input_comps + n, temp_comps + n, matrix, imageParam.MV_Dimension, 1);
            for (int n = 0; n < imageParam.num_comp; n++)
                vertical(temp_comps + n, output_comps + n, matrix, imageParam.MV_Dimension, 1);
            state = (imageParam.gradientFlag) ? IMAGE_GRADIENT : OUTPUT_PICTURE;
            delete[] temp_comps;
            break;
        case IMAGE_GRADIENT:
            printf("Gradient\r\n");
#if LAPLACIAN == 1
            //for (int i = 0; i < 3; i++) {
            //    for (int j = 0; j < 3; j++) {
            //        lapkernel[i*3 + j] = laplacianKernel[i][j];
            //    }
            //}
            //temp_comps = new  my_image_comp[imageParam.num_comp];
            //for (int n = 0; n < imageParam.num_comp; n++)
            //    temp_comps[n].init(imageParam.height, imageParam.width, 1); // Don't need a border for output
            //for (int n = 0; n < imageParam.num_comp; n++)
            //    apply_filter_modified(output_comps+n, temp_comps +n, lapkernel,3);
            ///*for (int n = 0; n < imageParam.num_comp; n++)
            //    apply_filter_modified(temp_comps + n, output_comps + n, lapkernel, 3);*/
            // for (int i = 0; i < 3; i++) {
            //    for (int j = 0; j < (temp_comps[i].stride * (temp_comps[i].height + 2 * temp_comps[i].border)); j++) {
            //        output_comps[i].handle[j] = imageParam.alpha*(temp_comps[i].handle[j] + 128);
            //        CLAMP_TO_BYTE(output_comps[i].handle[j]);
            //    }
            // }

            /************************************************************************************/
            temp_comps = new  my_image_comp[imageParam.num_comp];
            for (int n = 0; n < imageParam.num_comp; n++)
                temp_comps[n].init(imageParam.height, imageParam.width, 1); // Don't need a border for output
            temp_comps_2 = new  my_image_comp[imageParam.num_comp];
            for (int n = 0; n < imageParam.num_comp; n++)
                temp_comps_2[n].init(imageParam.height, imageParam.width, 1); // Don't need a border for output
            /************************************************************************************/
            for (int n = 0; n < imageParam.num_comp; n++)//f11
                temp_comps[n].SecondGrradientHorizontalFilter(output_comps + n, 3, &imageParam, imageParam.alpha);
            for (int n = 0; n < imageParam.num_comp; n++)//f22
                temp_comps_2[n].SecondGrradientverticalFilter(output_comps + n, 3, &imageParam, imageParam.alpha);
            //alpha*(f11 + f22) + 128
            for (int i = 0; i < 3; i++) {
                for (int r = 0; r < output_comps[i].height; r++) {  //进行卷积操作
                    for (int c = 0; c < output_comps[i].width; c++)
                    {
                        float* f11 = temp_comps[i].buf + r * temp_comps[i].stride + c;
                        float* f22 = temp_comps_2[i].buf + r * temp_comps_2[i].stride + c;
                        float* op = output_comps[i].buf + r * output_comps[i].stride + c;
                        float sum = (imageParam.alpha*(*f11 + *f22) + 128);
                        CLAMP_TO_BYTE(sum);
                        *op = sum;
                    }
                }
                
            }


            /************************************************************************************/
#else 
            temp_comps = new  my_image_comp[imageParam.num_comp];
            for (int n = 0; n < imageParam.num_comp; n++)
                temp_comps[n].init(imageParam.height, imageParam.width, 1); // Don't need a border for output
            for (int n = 0; n < imageParam.num_comp; n++)
                temp_comps[n].GrradientHorizontalFilter(output_comps + n, 3, imageParam.alpha);
            for (int n = 0; n < imageParam.num_comp; n++)
                output_comps[n].GrradientverticalFilter(temp_comps + n, 3, imageParam.alpha);
#endif
            state = OUTPUT_PICTURE;
            break;
        case OUTPUT_PICTURE:
#if LAPLACIAN == 1
            //OutputImage(&out, input_comps, &temp_comps, &line, &imageParam, argv);
            OutputImage(&out, input_comps, &output_comps, &line, &imageParam, argv);
            bmp_in__close(&in);
            delete[] line;
            delete[] input_comps;
            delete[] temp_comps;

#else
            OutputImage(&out, input_comps, &output_comps, &line, &imageParam, argv);
            bmp_in__close(&in);
            delete[] line;
            delete[] input_comps;
            delete[] output_comps;
#endif


            return 1;
        default:
            break;
        }
    }
}
//for (int i = 0; i < 3; i++) {
//    for (int j = 0; j < 3; j++) {
//        lapkernel[i][j] = laplacianKernel[i][j];
//    }
//}
//temp_comps = new  my_image_comp[imageParam.num_comp];
//for (int n = 0; n < imageParam.num_comp; n++)
//    temp_comps[n].init(imageParam.height, imageParam.width, 1); // Don't need a border for output
//for (int n = 0; n < imageParam.num_comp; n++)
//    apply_filter_modified(output_comps+n, temp_comps +n, lapkernel,3);
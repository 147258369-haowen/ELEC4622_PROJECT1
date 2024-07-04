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
    float* lapkernel = new float[3 * 3];
    float* LOGBuffer = nullptr;
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
            state = (imageParam.gradientFlag) ? THREE_CHANNEL : OUTPUT_PICTURE;
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
            state = (imageParam.gradientFlag) ? THREE_CHANNEL : OUTPUT_PICTURE;
            delete[] temp_comps;
            break;
        case LOG:
            LOGBuffer = input_comps->LOGKernelCreat(imageParam.GaussianDimension);
            input_comps->LoGKernelLoadValue(LOGBuffer, imageParam.GaussianDimension, sigma);
            for (int n = 0; n < imageParam.num_comp; n++)
                output_comps[n].init(imageParam.height, imageParam.width, 1); // Don't need a border for output
            for (int n = 0; n < imageParam.num_comp; n++)
                apply_filter_modified(input_comps + n, output_comps + n, LOGBuffer, imageParam.GaussianDimension);
            state = OUTPUT_PICTURE;
            //for (int i = 0; i < imageParam.GaussianDimension; i++) {
            //    for (int j = 0; j < imageParam.GaussianDimension; j++) {
            //        printf("%f, ",LOGBuffer[i* imageParam.GaussianDimension + j]);
            //    }
            //    printf("\n");
            //}
            break;
        case THREE_CHANNEL://0 -> Blue 1->Green 2-> Red
            //for (int n = 0; n < imageParam.num_comp; n++)
            //    output_comps[n].init(imageParam.height, imageParam.width, 1);
            ///**********BLUE YB[n] = X[n]*********************************/    



            ///**********GREEN*********************************/


            temp_comps = new  my_image_comp;
            temp_comps->init(imageParam.height, imageParam.width, 1); // Don't need a border for output
            temp_comps_2 = new  my_image_comp;
            temp_comps_2->init(imageParam.height, imageParam.width, 1); // Don't need a border for output
            temp_comps->SecondGrradientHorizontalFilter(output_comps + 1, 3, &imageParam, imageParam.beta);
            temp_comps_2->SecondGrradientverticalFilter(output_comps+1, 3, &imageParam, imageParam.beta);
            for (int r = 0; r < output_comps[1].height; r++) {
                for (int c = 0; c < output_comps[1].width; c++)
                {
                    float* a = temp_comps->buf + r * temp_comps->stride + c;
                    float* b = temp_comps_2->buf + r * temp_comps_2->stride + c;
                    float* output = output_comps[1].buf + r * output_comps[1].stride + c;                
                    float sum = (imageParam.beta * (*a + *b) + 128);
                    CLAMP_TO_BYTE(sum);
                    *output = sum;

                }
            }
            delete temp_comps;
            delete temp_comps_2;
            /**********RED*********************************/
            temp_comps = new  my_image_comp;
            temp_comps->init(imageParam.height, imageParam.width, 1); // Don't need a border for output
            temp_comps->GradientFilter(output_comps + 2,3, imageParam.alpha);
            for (int r = 0; r < output_comps[2].height; r++) {  //进行卷积操作
                for (int c = 0; c < output_comps[2].width; c++)
                {
                    float* a = temp_comps->buf + r * temp_comps->stride + c;
                    float* output = output_comps[2].buf + r * output_comps[2].stride + c;
                    *output = *a;

                }
            }
            delete temp_comps;
            state = OUTPUT_PICTURE;
            break;
        case IMAGE_GRADIENT:
            printf("Gradient\r\n");
#if LAPLACIAN == 1
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

            temp_comps = new  my_image_comp[imageParam.num_comp];
            for (int n = 0; n < imageParam.num_comp; n++)
                temp_comps[n].init(imageParam.height, imageParam.width, 1); // Don't need a border for output
            /************************************************************************************/
            for (int n = 0; n < imageParam.num_comp; n++)
                temp_comps[n].SecondGrradientHorizontalFilter(output_comps + n, 3, &imageParam, imageParam.alpha);
            for (int n = 0; n < imageParam.num_comp; n++)
                output_comps[n].SecondGrradientverticalFilter(temp_comps + n, 3, &imageParam, imageParam.alpha);
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
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
    my_image_comp* temp_comps_out = nullptr;
    ImageParam imageParam;
    float sigma = 0;
    float** lapkernel = allocateMatrix(3);
    
    while (true) {
        switch (state) {
        case CHECK_INPUT:
            CheckInput(argc, argv,&sigma,&filterChoose,&imageParam);
            state = filterChoose ? LOAD_PICTURE : MOVING_AVERAGE_CHECK;
            break;
        case LOAD_PICTURE:
            //printf("%f\r\n", FilterInit(matrix, HIGHT, WIDTH));//归一化
            imageParam.GaussianDimension = GaussianWindowDimensionChoose(sigma);
            LoadImage(&in, &input_comps, &output_comps,&line, &imageParam, &filterChoose,argv);
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
                horizontal(input_comps + n, temp_comps + n, matrix, imageParam.GaussianDimension,0);
            for (int n = 0; n < imageParam.num_comp; n++)
                vertical(temp_comps + n, output_comps + n, matrix, imageParam.GaussianDimension,0);
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
                temp_comps[n].init(imageParam.height, imageParam.width, (imageParam.MV_Dimension - 1)/2); // Don't need a border for output
            for (int n = 0; n < imageParam.num_comp; n++)
                horizontal(input_comps + n, temp_comps + n, matrix, imageParam.MV_Dimension,1);
            for (int n = 0; n < imageParam.num_comp; n++)
                vertical(temp_comps + n, output_comps + n, matrix, imageParam.MV_Dimension,1);
            state =(imageParam.gradientFlag)? IMAGE_GRADIENT:OUTPUT_PICTURE;
            delete[] temp_comps;
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
                temp_comps[n].SecondGrradientHorizontalFilter(output_comps + n, 3,&imageParam, imageParam.alpha);
            for (int n = 0; n < imageParam.num_comp; n++)
                output_comps[n].SecondGrradientverticalFilter(temp_comps + n, 3,&imageParam, imageParam.alpha);
            /************************************************************************************/
#else 
            temp_comps = new  my_image_comp[imageParam.num_comp];
            for (int n = 0; n < imageParam.num_comp; n++)
                temp_comps[n].init(imageParam.height, imageParam.width, 1); // Don't need a border for output
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < (temp_comps[i].stride * (temp_comps[i].height + 2 * temp_comps[i].border)); j++) {
                    temp_comps[i].handle[j] = output_comps[i].handle[j];
                }
            }
             for (int n = 0; n < imageParam.num_comp; n++)
                output_comps[n].GradientFilter(temp_comps + n,3, imageParam.alpha);
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
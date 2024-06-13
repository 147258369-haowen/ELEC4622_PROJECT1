/*****************************************************************************/
// File: filtering_main.cpp
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#include "io_bmp.h"
#include "image_comps.h"
#include <iostream>
/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/*****************************************************************************/
/*                  my_image_comp::perform_boundary_extension                */
/*****************************************************************************/

void my_image_comp::perform_boundary_extension()
{
    int r, c;
    // First extend upwards
    float* first_line = buf;
    for (r = 1; r <= border; r++)
        for (c = 0; c < width; c++)
            first_line[-r * stride + c] = first_line[c];
    // Now extend downwards
    float* last_line = buf + (height - 1) * stride;
    for (r = 1; r <= border; r++)
        for (c = 0; c < width; c++)
            last_line[r * stride + c] = last_line[c];
    // Now extend all rows to the left and to the right
    float* left_edge = buf - border * stride;
    float* right_edge = left_edge + width - 1;
    for (r = height + 2 * border; r > 0; r--, left_edge += stride, right_edge += stride)
        for (c = 1; c <= border; c++)
        {
            left_edge[-c] = left_edge[0];
            right_edge[c] = right_edge[0];
        }
}


/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */

/*****************************************************************************/
/*                                apply_filter                               */
/*****************************************************************************/

void apply_filter(my_image_comp* in, my_image_comp* out)
{
#define FILTER_EXTENT 4//向左或者向右延申的长度
#define FILTER_DIM (2*FILTER_EXTENT+1)//卷积核宽
#define FILTER_TAPS (FILTER_DIM*FILTER_DIM)//卷积核点的数量

    // Create the filter kernel as a local array on the stack, which can accept
    // row and column indices in the range -FILTER_EXTENT to +FILTER_EXTENT.
    float filter_buf[FILTER_TAPS];
    float* mirror_psf = filter_buf + (FILTER_DIM * FILTER_EXTENT) + FILTER_EXTENT;
    // `mirror_psf' points to the central tap in the filter
    int r, c;
    for (r = -FILTER_EXTENT; r <= FILTER_EXTENT; r++)
        for (c = -FILTER_EXTENT; c <= FILTER_EXTENT; c++)
            mirror_psf[r * FILTER_DIM + c] = 1.0F / FILTER_TAPS;

    // Check for consistent dimensions
    assert(in->border >= FILTER_EXTENT);
    assert((out->height <= in->height) && (out->width <= in->width));

    // Perform the convolution
    for (r = 0; r < out->height; r++)
        for (c = 0; c < out->width; c++)
        {
            float* ip = in->buf + r * in->stride + c;
            float* op = out->buf + r * out->stride + c;
            float sum = 0.0F;
            for (int y = -FILTER_EXTENT; y <= FILTER_EXTENT; y++)
                for (int x = -FILTER_EXTENT; x <= FILTER_EXTENT; x++)
                    sum += ip[y * in->stride + x] * mirror_psf[y * FILTER_DIM + x];
            *op = sum;
        }
}
float FilterInit(float** input, int height, int width) {
    float sum = 0;
    float temp = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            sum += input[i][j];
        }
    }
    if (sum != 1) {
        temp = (1 / sum);
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                input[i][j] = input[i][j] * temp;
            }
        }
        return temp;

    }
    return 1;
}
void apply_filter_modified(my_image_comp* in, my_image_comp* out, float** inputfilter, int width) {
    int filter_extent = (width - 1) / 2;
    int filter_dim = width;
    int filter_taps = width * width;
    float* filter_buf = new float[filter_taps];
    float* mirror_psf = filter_buf + (filter_dim * filter_extent) + filter_extent;//中间点
    for (int i = -filter_extent; i <= filter_extent; i++) {//加载卷积核
        for (int j = -filter_extent; j <= filter_extent; j++) {
            mirror_psf[i * filter_dim + j] = inputfilter[i + filter_extent][j + filter_extent];
        }
    }
    for (int r = 0; r < out->height; r++)//进行卷积操作
        for (int c = 0; c < out->width; c++)
        {
            float* ip = in->buf + r * in->stride + c;
            float* op = out->buf + r * out->stride + c;
            float sum = 0.0F;
            for (int y = -filter_extent; y <= filter_extent; y++)//列
                for (int x = -filter_extent; x <= filter_extent; x++)//行
                    sum += ip[y * in->stride + x] * mirror_psf[y * filter_dim + x];
            *op = sum;
        }
    delete[] filter_buf;
}
void unsharp_mask_filter(float** inputfilter, int width, float alpha) {
    float* impulseSignal = new float[width * width];

    for (int i = 0; i < width; i++) {
        for (int j = 0; j < width; j++) {
            impulseSignal[i * width + j] = 0;
            if (i == j && i == ((width - 1) / 2) && j == ((width - 1) / 2)) {
                impulseSignal[i * width + j] = 1;
            }
        }
    }

    for (int i = 0; i < width; i++) {
        for (int j = 0; j < width; j++) {
            float temp;
            temp = impulseSignal[i * width + j] + alpha * (impulseSignal[i * width + j] - inputfilter[i][j]);
            inputfilter[i][j] = temp;
        }
    }

    delete[] impulseSignal;
}

float h1[5][5] = {
    {0, 1 / 3.0, 1 / 2.0, 1 / 3.0, 0},
    {1 / 3.0, 1 / 2.0, 1, 0.5, 1 / 3.0},
    {1 / 2.0, 1, 1, 1, 1 / 2.0},
    {1 / 3.0, 1 / 2.0, 1, 1 / 2.0, 1 / 3.0},
    {0, 1 / 3.0, 1 / 2.0, 1 / 3.0, 0}
};
float h2[9][9]{
    {1 / 4.0,1 / 2.0,1 / 2.0,1 / 4.0,0,0,0,0,0},
    {1 / 2.0,1,1,1 / 2.0,0,0,0,0,0},
    {1 / 2.0,1,1,1 / 2.0,0,0,0,0,0},
    {1 / 4.0,1 / 2.0,1 / 2.0,1 / 4.0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0},
};
float h3[9][9]{
    {0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,1 / 4.0,1 / 2.0,1 / 2.0,1 / 4.0},
    {0,0,0,0,0,1 / 2.0,1,1,1 / 2.0},
    {0,0,0,0,0,1 / 2.0,1,1,1 / 2.0},
    {0,0,0,0,0,1 / 4.0,1 / 2.0,1 / 2.0,1 / 4.0},
};
#define WIDTH 5
#define HIGHT 5
/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

void CheckInput(int argc, char* argv[],float* sigma ,int* filterChooseFlag) {
    if (argc < 4)
    {
        fprintf(stderr, "Usage: %s <in bmp file> <out bmp file> <sigma> (optional)<-w>\n", argv[0]);
        exit(-1);
    }
    else if (atof(argv[3]) <= 0.0) {
        fprintf(stderr, "The sigma value is not correct\n");
        exit(-1);
    }
    if (argc == 4) {
        *sigma = atof(argv[3]);
        fprintf(stdout, "Gaussian\n");
        *filterChooseFlag = GAUSSIAN;
    }
    else if (argc == 5 && strcmp(argv[4],"-w\n")) {
        *sigma = atof(argv[3]);
        fprintf(stdout, "Moving average\n");
        *filterChooseFlag = MOVINGAVERAGE;
    }
    
}
void LoadImage(bmp_in* in, my_image_comp** input_comps, my_image_comp** output_comps,io_byte** line, 
    ImageParam* imageParam, int* filterChoose,char** argv) {
    int err_code = 0;
    try {
        if ((err_code = bmp_in__open(in, argv[1])) != 0)
            throw err_code;
        int width = in->cols, height = in->rows;
        int n, num_comps = in->num_components;
        imageParam->height = height;
        imageParam->width = width;
        imageParam->num_comp = num_comps;
        printf("height: %d\n", height);
        printf("width: %d\n", width);
        printf("num_comps: %d\n", num_comps);
        *input_comps = new my_image_comp[num_comps];
        *output_comps = new my_image_comp[num_comps];
        if (*filterChoose) {//gaussian
            for (n = 0; n < num_comps; n++)
                (*input_comps)[n].init(height, width, 4); // Leave a border of 4
        }
        else {
            for (n = 0; n < num_comps; n++)
                (*input_comps)[n].init(height, width, (imageParam->MV_Dimension - 1)/2); // Leave a border of 4
        }
        printf("stride: %d\n", (*input_comps)->stride);
        int r; // Declare row index
        *line = new io_byte[width * num_comps];
        for (r = height - 1; r >= 0; r--) {
            // "r" holds the true row index we are reading, since the image is
            // stored upside down in the BMP file.
            if ((err_code = bmp_in__get_line(in, *line)) != 0)
                throw err_code;
            for (n = 0; n < num_comps; n++) {
                io_byte* src = *line + n; // Points to first sample of component n
                float* dst = (*input_comps)[n].buf + r * (*input_comps)[n].stride;
                for (int c = 0; c < width; c++, src += num_comps)
                    dst[c] = (float)*src; // The cast to type "float" is not
                // strictly required here, since bytes can always be
                // converted to floats without any loss of information.
            }
        }
        //bmp_in__close(in);
    }
    catch (int exc) {
        if (exc == IO_ERR_NO_FILE)
            fprintf(stderr, "Cannot open supplied input or output file.\n");
        else if (exc == IO_ERR_FILE_HEADER)
            fprintf(stderr, "Error encountered while parsing BMP file header.\n");
        else if (exc == IO_ERR_UNSUPPORTED)
            fprintf(stderr, "Input uses an unsupported BMP file format.\n  Current "
                "simple example supports only 8-bit and 24-bit data.\n");
        else if (exc == IO_ERR_FILE_TRUNC)
            fprintf(stderr, "Input or output file truncated unexpectedly.\n");
        else if (exc == IO_ERR_FILE_NOT_OPEN)
            fprintf(stderr, "Trying to access a file which is not open!(?)\n");
    }
}
int OutputImage(bmp_out* out, my_image_comp* input_comps, my_image_comp** output_comps, io_byte** line, ImageParam* imageParam, char** argv) {
    int err_code = 0;
    try {
        printf("enter output\n");
        int width = imageParam->width, height = imageParam->height;
        int n, num_comps = imageParam->num_comp;
        //*output_comps = new my_image_comp[num_comps];
        //for (n = 0; n < num_comps; n++)
        //    (*output_comps)[n].init(height, width, 0); // Don't need a border for output

        printf("height: %d\n", (*output_comps)[0].height);
        printf("width: %d\n", (*output_comps)[0].width);
        printf("stride: %d\n", (*output_comps)[0].stride);
        if ((err_code = bmp_out__open(out, argv[2], width, height, num_comps)) != 0)
            throw err_code;
        int r;
        for (r = height - 1; r >= 0; r--) {
            // "r" holds the true row index we are writing, since the image is
            // written upside down in BMP files.
            for (n = 0; n < num_comps; n++) {
                io_byte* dst = *line + n; // Points to first sample of component n
                float* src = (*output_comps)[n].buf + r * (*output_comps)[n].stride;
                for (int c = 0; c < width; c++, dst += num_comps) {
                    *dst = (io_byte)src[c];
                }
            }
            bmp_out__put_line(out, *line);
        }
        bmp_out__close(out);

        printf("Output success\n");
    }
    catch (int exc) {
        if (exc == IO_ERR_NO_FILE)
            fprintf(stderr, "Cannot open supplied input or output file.\n");
        else if (exc == IO_ERR_FILE_HEADER)
            fprintf(stderr, "Error encountered while parsing BMP file header.\n");
        else if (exc == IO_ERR_UNSUPPORTED)
            fprintf(stderr, "Input uses an unsupported BMP file format.\n  Current "
                "simple example supports only 8-bit and 24-bit data.\n");
        else if (exc == IO_ERR_FILE_TRUNC)
            fprintf(stderr, "Input or output file truncated unexpectedly.\n");
        else if (exc == IO_ERR_FILE_NOT_OPEN)
            fprintf(stderr, "Trying to access a file which is not open!(?)\n");
        return -1;
    }
    return 0;
}
float GaussianFillKernel(int x, int y, float sigma) {
    float coffienct = 1.0f / (2.0f*PI*sigma*sigma);
    float e_part = exp(-(x * x + y * y) / 2.0f * sigma * sigma);
    return coffienct * e_part;
}

int LoadGaussianValue(float** matrix, float sigma) {
    int x_offset =((WIDTH - 1) / 2);
    int y_offset = ((HIGHT - 1) / 2);
    int x = (WIDTH - 1);
    int y = 0;
    for (int i = -x_offset; i <= x_offset; i++) {
        for (int j = -y_offset; j <= y_offset; j++) {
            matrix[x][y] = GaussianFillKernel(i, j,sigma);
            printf("%f,", matrix[x][y]);
            x --;
        }
        printf("\n");
        x = (WIDTH - 1);
        y++;
    }

    return 1;
}

void FreeMatrix(float** matrix, int width) {
    for (int i = 0; i < width; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
}
float** allocateMatrix(int dimension) {
    float** matrix = new float* [dimension];
    for (int i = 0; i < dimension; i++) {
        matrix[i] = new float[dimension];
        for (int j = 0; j < dimension; j++) {
            matrix[i][j] = 0.0f;
        }
    }
    return matrix;
}

static float MovingAverageMapDimension2variance(int dimension,int mode) {
    int x_offset = (dimension - 1) / 2;
    int y_offset = (dimension - 1) / 2;
    float temp = 0.0f;
    float weight = 1.0f / (dimension * dimension);
    for (int i = -x_offset; i <= x_offset; i++) {
        for (int j = -y_offset; j <= y_offset; j++) {
            temp += weight*((float)i * (float)i + (float)j * (float)j);
        }
    }
    float result =  temp;
    if (mode == 1) {
        printf("dimension: %d ,variance: %f,root variance: %f\n", dimension, result, sqrt(result));
    }
    return result;
}

static bool Compare(float variance,float sigma) {
    return (variance >= (sigma*sigma)) ? true : false;
}
static float Abs(float a, float b) {
    float temp = a - b;
    if (temp <= 0) {
        return temp = -temp;  
    }
    else {
        return temp;
    }
}
//Compare the variance with the square of sigma and find the dimension of the sliding 
//average kernel corresponding to the variance closest to the square of sigma
int VarianceLoopCheck(float sigma, ImageParam* imageParam) {
    int dimension = 3;
    float variance = 0.0f;
    bool check = false;
    float sigmaSquare = (sigma * sigma);
    while (1) {
        variance = MovingAverageMapDimension2variance(dimension,1);
        check = Compare(variance, sigma);
        if (check == false) dimension += 2;
        else {
            int previousDimension = (dimension - 2);
            int variancePreviousOne = MovingAverageMapDimension2variance(previousDimension,0);
            int differenceCurrent = Abs(variance,sigmaSquare);
            int differencePrevious = Abs(variancePreviousOne,sigmaSquare);
            float closestVariance = (differenceCurrent >= differencePrevious) ? variancePreviousOne : variance;
            int output = (differenceCurrent >= differencePrevious) ? previousDimension : dimension;
            printf("sigma: %f,the closest dimension is %d,variance:%f\r\n", sigma,output, closestVariance);
            imageParam->MV_Dimension = output;
            return output;
        }
    }
}
void MovingAverageSetValue(float** matrix,int dimension) {
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            matrix[i][j] = 1.0f / ((float)dimension * (float)dimension);
        }
    }
}
//{ 0, 1, 2, 3, 4 }
//{ 0, 1, 2, 3, 4 }
//{ 0, 1, 2, 3, 4 }
//{ 0, 1, 2, 3, 4 }
//{ 0, 1, 2, 3, 4 }
//[5,0][4,0],3 0
//[5,1][4]
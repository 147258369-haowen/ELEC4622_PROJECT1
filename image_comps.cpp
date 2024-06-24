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
#include <emmintrin.h> // Include SSE2 processor intrinsic functions
#include <pmmintrin.h>
/* ========================================================================= */
/*                 Implementation of `my_image_comp' functions               */
/* ========================================================================= */

/*****************************************************************************/
/*                  my_image_comp::perform_boundary_extension                */
/*****************************************************************************/
float h1D[3] = { -0.5, 0.0, 0.5 };
float h2D[3] = { -0.5, 0.0, 0.5 };
float laplacianKernel[3][3] = {
        { 0.0f,-1.0f,0.0f},
        { -1.0f, 4.0f, -1.0f },
        { 0.0f,-1.0f,0.0f}
};
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

void listshift(float* buffer, float** ptr, int width) {
    static int flag = 0;
    if (flag == (width - 1)) {
        flag = 0;
    }
    else {
        flag++;
    }
    *ptr = buffer + flag;

}
void horizontal(my_image_comp* in, my_image_comp* out, float** inputfilter, int width,int G_MF_flag) {

    int filter_extent = (width - 1) / 2;
    int filter_dim = width;
    int filter_taps = width * width;
    float* filter_buf = new float[filter_taps];
    float* mirror_psf = filter_buf + (filter_dim * filter_extent) + filter_extent;//中间点
    float temp = 0.0f;
    for (int i = -filter_extent; i <= filter_extent; i++) {//加载卷积核
        for (int j = -filter_extent; j <= filter_extent; j++) {
            mirror_psf[i * filter_dim + j] = inputfilter[i + filter_extent][j + filter_extent];
        }
    }
    
    for (int j = -filter_extent; j <= filter_extent; j++) {
        temp += mirror_psf[j];
    }
    if (temp > 1.03f || temp < 0.98) {
        for (int j = -filter_extent; j <= filter_extent; j++) {
            mirror_psf[j] = (mirror_psf[j] * (1.0f / temp));
        }
    }
    int inital_flag = 0;
    float sum = 0.0F;
    float* buffer = new float[width];
    float* bufptr = buffer;
    for (int r = 0; r < out->height; r++)//进行卷积操作

        for (int c = 0; c < out->width; c++)
        {
            float* ip = in->buf + r * in->stride + c;
            float* op = out->buf + r * out->stride + c;
            //mirror_psf = mirror_psf + r * out->stride + c;
            if (!G_MF_flag) {//gaussian
                sum = 0.0F;
            }        
           // for (int y = -filter_extent; y <= filter_extent; y++)//列
            if (inital_flag == 0) {
                for (int x = -filter_extent; x <= filter_extent; x++)//行
                {
                    sum += ip[x] * mirror_psf[x];
                    *bufptr = ip[x] * mirror_psf[x];
                    listshift(buffer,&bufptr,width);
                }

            }
            else {
                sum += (ip[filter_extent] * mirror_psf[filter_extent]);
                sum -= *bufptr;
                *bufptr = (ip[filter_extent] * mirror_psf[filter_extent]);
                listshift(buffer, &bufptr, width);
            }
            if (sum > 255) sum = 255;
            else if (sum < 0) sum = 0;
            inital_flag = G_MF_flag ? 1 : 0;
            *op = sum;
            //if (c == (out->width - 1)) {
            //    sum = 0;
            //    inital_flag = 0;
            //}
        }
    delete[] filter_buf;
}
void vertical(my_image_comp* in, my_image_comp* out, float** inputfilter, int width, int G_MF_flag) {

    int filter_extent = (width - 1) / 2;
    int filter_dim = width;
    int filter_taps = (2* filter_extent + 1);
    float* filter_buf = new float[filter_taps];
    float* mirror_psf = filter_buf  + filter_extent;//中间点
    float temp = 0.0f;
    for (int j = -filter_extent; j <= filter_extent; j++) {
        mirror_psf[j] = inputfilter[(width-1)/2][j + filter_extent];
        temp += mirror_psf[j];
    }
    if (temp > 1.03f || temp <0.98) {
        for (int j = -filter_extent; j <= filter_extent; j++) {
            mirror_psf[j] = (mirror_psf[j] * (1.0f / temp));
        }   
    }
    float* buffer = new float[width];
    float* bufptr = buffer;
    int inital_flag = 0;
    float sum = 0.0F;
    for (int r = 0; r < out->width; r++)//进行卷积操作
        for (int c = 0; c < out->height; c++)
        {
            //float* ip = in->buf + r * in->stride + c;
            //float* op = out->buf + r * out->stride + c;
            float* ip = in->buf + c * in->stride + r;
            float* op = out->buf + c * out->stride + r;

            if (!G_MF_flag) {
                sum = 0.0F;
            }
            if (!inital_flag) {
                for (int y = -filter_extent; y <= filter_extent; y++)//
                {
                    
                    sum += ip[y * in->stride] * mirror_psf[y];
                    *bufptr = ip[y * in->stride] * mirror_psf[y];
                    listshift(buffer, &bufptr, width);
                }
            }
            else {
                sum += ip[filter_extent * in->stride] * mirror_psf[filter_extent];
                sum -= *bufptr;
                *bufptr = (ip[filter_extent * in->stride] * mirror_psf[filter_extent]);
                listshift(buffer, &bufptr, width);
            }
            if (sum > 255) sum = 255;
            else if (sum < 0) sum = 0;
            inital_flag = G_MF_flag ? 1:0;
            *op = sum;
  /*          if (c == (out->height - 1)) {
                sum = 0;
                inital_flag = 0;
            }*/
        }
    delete[] filter_buf;
}
float FilterNormalized(float** input, int dimension) {
    float sum = 0;
    float temp = 0;
#ifdef DEBUG
    printf("After mormalize:\r\n");
#endif
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            sum += input[i][j];
        }
    }
    for (int i = 0; i < dimension; ++i) {
        for (int j = 0; j < dimension; ++j) {
            input[i][j] /= sum;
#ifdef DEBUG
            printf("%f,", input[i][j]);
#endif 
        }
#ifdef DEBUG
        printf("\n");
#endif 
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
void my_image_comp:: apply_filter_modified_simo(my_image_comp* in, my_image_comp* out, float** inputfilter, int width) {
    int filter_extent = (width - 1) / 2;
    int filter_dim = width;
    int filter_taps = width * width;

    __m128* filter_buf = (__m128*) _mm_malloc(filter_taps * sizeof(__m128), 16);
    __m128* mirror_psf = filter_buf + (filter_dim * filter_extent) + filter_extent; // 中间点

    for (int i = -filter_extent; i <= filter_extent; i++) { // 加载卷积核
        for (int j = -filter_extent; j <= filter_extent; j++) {
            float filter_value = inputfilter[i + filter_extent][j + filter_extent];
            mirror_psf[i * filter_dim + j] = _mm_set_ps1(filter_value); // 将单个float值转换为__m128
        }
    }

    // 卷积操作
    for (int r = 0; r < out->height; r++) {
        for (int c = 0; c < out->width; c++) {
            float* ip = in->buf + r * in->stride + c;
            float* op = out->buf + r * out->stride + c;
            __m128 sum = _mm_setzero_ps(); // 初始化为0的__m128类型
            for (int y = -filter_extent; y <= filter_extent; y++) { // 列
                for (int x = -filter_extent; x <= filter_extent; x++) { // 行
                    // 处理边界条件
                    int yy = r + y;
                    int xx = c + x;
                    if (yy < 0) yy = 0;
                    if (yy >= in->height) yy = in->height - 1;
                    if (xx < 0) xx = 0;
                    if (xx >= in->width) xx = in->width - 1;

                    float pixel_value = in->buf[yy * in->stride + xx];
                    __m128 ip_val = _mm_set_ps1(pixel_value); // 将ip的值加载为__m128
                    sum = _mm_add_ps(sum, _mm_mul_ps(ip_val, mirror_psf[y * filter_dim + x]));
                }
            }
            // 将__m128的值合并成一个float
            float result[4];
            _mm_storeu_ps(result, sum);
            *op = result[0] + result[1] + result[2] + result[3];
        }
    }

    // 使用完毕后记得释放内存
    _mm_free(filter_buf);
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

void CheckInput(int argc, char* argv[],float* sigma ,int* filterChooseFlag, ImageParam* param) {
    param->gradientFlag = false;
    if (argc < 4)
    {
        fprintf(stderr, "Usage: %s <in bmp file> <out bmp file> <sigma> (optional)<-w>\n", argv[0]);
        exit(-1);
    }
    else if (atof(argv[3]) <= 0.0) {
        fprintf(stderr, "The sigma value is not correct\n");
        exit(-1);
    }
    *sigma = atof(argv[3]);
    if (argc == 4) {
       
        fprintf(stdout, "Gaussian\n");
        *filterChooseFlag = GAUSSIAN;
    }
    else if (argc == 5 && !strcmp(argv[4],"-w")) {
        printf("1\r\n");
        fprintf(stdout, "Moving average\n");
        *filterChooseFlag = MOVINGAVERAGE;
    }
    else if (argc == 5 && (atof(argv[4]) > 0.0f)) {
        printf("2\r\n");
        param->alpha = atof(argv[4]);
        if (param->alpha <= 0.0f) {
            fprintf(stderr, "The alpha value is not correct\n");
            exit(-1);
        }
        *filterChooseFlag = GAUSSIAN;
        param->gradientFlag = true;
    }
    if (argc == 6) {
        printf("3\r\n");
        param->alpha = atof(argv[5]);
        if (param->alpha <= 0.0f) {
            fprintf(stderr, "The alpha value is not correct\n");
            exit(-1);
        }
        param->gradientFlag = true;
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
                (*input_comps)[n].init(height, width, (imageParam->GaussianDimension-1)/2); // Leave a border of 4
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
    float coefficient = 1.0f / (2.0f * PI * sigma * sigma);
    float e_part = exp(-((float)x * (float)x + (float)y * (float)y) / (2.0f * sigma * sigma));
    return coefficient * e_part;
}
int GaussianWindowDimensionChoose(float sigma) {
    float windowSize = 2 * (3 * sigma) + 1;
    int dimension = (int)windowSize;
    return dimension;
}
int LoadGaussianValue(float** matrix, float sigma,int dimension) {
    int x_offset =((dimension - 1) / 2);
    int y_offset = ((dimension - 1) / 2);
    int x = (dimension - 1);
    int y = 0;
    for (int i = -x_offset; i <= x_offset; i++) {
        for (int j = -y_offset; j <= y_offset; j++) {
            matrix[x][y] = GaussianFillKernel(i, j,sigma);
#ifdef DEBUG
            printf("%f,", matrix[x][y]);
#endif
            x --;
        }
#ifdef DEBUG
        printf("\n");
#endif
        x = (dimension - 1);
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
void my_image_comp::vector_filter(my_image_comp* in,int dimension)
{
//#define FILTER_EXTENT 4
//#define FILTER_TAPS (2*FILTER_EXTENT+1)
    int FILTER_EXTENT_ = (dimension-1)/2;//半径
    int FILTER_TAPS_ = (2 * FILTER_EXTENT_ + 1);
    // Create the vertical filter PSF as a local array on the stack.
    __m128* filter_buf = (__m128*)_mm_malloc(FILTER_TAPS_ * sizeof(__m128), 16);
    
    __m128* mirror_psf = filter_buf + FILTER_EXTENT_;//中心点
    // `mirror_psf' points to the central tap in the filter
    for (int t = -FILTER_EXTENT_; t <= FILTER_EXTENT_; t++)
        mirror_psf[t] = _mm_set1_ps(1.0F / FILTER_TAPS_);

    // Check for consistent dimensions
    assert(in->border >= FILTER_EXTENT_);
    assert((this->height <= in->height) && (this->width <= in->width));
    assert(((stride & 3) == 0) && ((in->stride & 3) == 0));
    int vec_stride_in = in->stride / 4;
    int vec_stride_out = this->stride / 4;
    int vec_width_out = (this->width + 3) / 4; // Big enough to cover the width

    // Do the filtering
    __m128* line_out = (__m128*) buf;
    __m128* line_in = (__m128*)(in->buf);
    for (int r = 0; r < height; r++,
        line_out += vec_stride_out, line_in += vec_stride_in)
        for (int c = 0; c < vec_width_out; c++)
        {
            __m128* ip = (line_in + c) - vec_stride_in * FILTER_EXTENT_;
            __m128 sum = _mm_setzero_ps();
            for (int y = -FILTER_EXTENT_; y <= FILTER_EXTENT_; y++, ip += vec_stride_in)
                sum = _mm_add_ps(sum, _mm_mul_ps(mirror_psf[y], *ip));
            line_out[c] = sum;
        }
}
void my_image_comp::vector_horizontal_filter(my_image_comp* in, int dimension)
{
    int FILTER_EXTENT_ = (dimension - 1) / 2; // 半径
    int FILTER_TAPS_ = (2 * FILTER_EXTENT_ + 1);

    // 动态分配垂直滤波PSF作为堆上的本地数组。
    __m128* filter_buf = (__m128*)_mm_malloc(FILTER_TAPS_ * sizeof(__m128), 16);
    __m128* mirror_psf = filter_buf + FILTER_EXTENT_; // 中心点

    // `mirror_psf' 指向滤波器的中心抽头
    for (int t = -FILTER_EXTENT_; t <= FILTER_EXTENT_; t++)
        mirror_psf[t] = _mm_set1_ps(1.0F / FILTER_TAPS_);

    // 检查维度是否一致
    assert(in->border >= FILTER_EXTENT_);
    assert((this->height <= in->height) && (this->width <= in->width));
    assert(((stride & 3) == 0) && ((in->stride & 3) == 0));
    int vec_stride_in = in->stride / 4;
    int vec_stride_out = this->stride / 4;
    int vec_width_out = (this->width + 3) / 4; // 足够覆盖宽度

    // 进行滤波操作
    __m128* line_out = (__m128*)buf;
    __m128* line_in = (__m128*)(in->buf);
    for (int r = 0; r < height; r++,
        line_out += vec_stride_out, line_in += vec_stride_in)
    {
        for (int c = 0; c < vec_width_out; c++)
        {
            __m128 sum = _mm_setzero_ps();
            for (int x = -FILTER_EXTENT_; x <= FILTER_EXTENT_; x++)
            {
                int index = c + x;
                // 确保索引在有效范围内
                if (index >= 0 && index < vec_width_out)
                {
                    __m128* ip = line_in + index;
                    sum = _mm_add_ps(sum, _mm_mul_ps(mirror_psf[x], *ip));
                }
            }
            line_out[c] = sum;
        }
    }

    // 释放动态分配的内存
    _mm_free(filter_buf);
}
void my_image_comp::SecondGrradientHorizontalFilter(my_image_comp* in, int dimension, ImageParam* imagP, int alpha) {
    int radius = (dimension - 1) / 2;
    float* centralPoint = &h1D[1];
    my_image_comp* temp = new my_image_comp;
    temp->init(this->height, this->width, 1);
    float* outBuf = temp->buf;
    float* inBuf = in->buf;
    for (int i = 1; i <= 2; i++) {
        for (int r = 0; r < this->height; r++) {  //进行卷积操作
            for (int c = 0; c < this->width; c++)
            {
                float* ip = inBuf + r * in->stride + c;
                float* op = outBuf + r * this->stride + c;
                float sum = 0.00F;
                for (int x = -radius; x <= radius; x++)//行
                {
                    float temp_ = ip[x] * centralPoint[x];
                    if (i == 1) {
                        sum += (temp_);//(float)alpha *
                    }
                    else {
                        sum += ((float)alpha * temp_ + 128);
                    }
                    
                }
                if (sum <= 0.0f) sum = 0.0f;
                else if (sum >= 255.0f) sum = 255.0f;
                *op = sum;
            }
        }
        outBuf = this->buf;
        inBuf = temp->buf;
    }
    delete temp;

}
void my_image_comp::SecondGrradientverticalFilter(my_image_comp* in, int dimension, ImageParam* imagP, int alpha) {
    float* centralPoint = &h2D[1];
    int radius = (dimension - 1) / 2;
    my_image_comp* temp_comp = new my_image_comp;
    temp_comp->init(this->height, this->width, 1);
    float* outBuf = temp_comp->buf;
    float* inBuf = in->buf;
    for (int i = 1; i <= 2; i++) {
        for (int r = 0; r < this->width; r++) {//进行卷积操作
            for (int c = 0; c < this->height; c++)
            {
                float* ip = inBuf + c * in->stride + r;
                float* op = outBuf + c * this->stride + r;

                float sum = 0.00F;

                for (int y = -radius; y <= radius; y++)//
                {
                    float temp = ip[y * in->stride] * centralPoint[y];
                    if (i == 1) {
                        sum += ( temp);//(float)alpha *
                    }
                    else {
                        sum += ((float)alpha * temp + 128);
                    }
                }
                if (sum <= 0.0f) sum = 0.0f;
                else if (sum >= 255.0f) sum = 255.0f;

                *op = sum;
            }
        }
        outBuf = this->buf;
        inBuf = temp_comp->buf;
    }
    delete temp_comp;
}
void my_image_comp::GrradientHorizontalFilter(my_image_comp* in, int dimension,int alpha) {
    int radius = (dimension - 1) / 2;
    float* centralPoint = &h1D[1];
    for (int r = 0; r < this->height; r++)//进行卷积操作
        for (int c = 0; c < this->width; c++)
        {
            float* ip = in->buf + r * in->stride + c;
            float* op = this->buf + r * this->stride + c;
            //mirror_psf = mirror_psf + r * out->stride + c;
            float sum = 0.00F;
            // for (int y = -filter_extent; y <= filter_extent; y++)//列
            for (int x = -radius; x <= radius; x++)//行
            {
                float temp = ip[x] * centralPoint[x];
                sum += ((float)alpha * temp);
            }
            if (sum <= 0.0f) sum = 0.1f;
            else if (sum >= 255.0f) sum = 254.0f;
            *op = sum;
        }

}
void my_image_comp::GrradientverticalFilter(my_image_comp* in, int width,int alpha) {
    float* centralPoint = &h2D[1];
    int radius = (width - 1) / 2;
    for (int r = 0; r < this->height; r++)//进行卷积操作
        for (int c = 0; c < this->width; c++)
        {
            float* ip = in->buf + r * in->stride + c;
            float* op = this->buf + r * this->stride + c;

            float sum = 0.00F;

            for (int y = -radius; y <= radius; y++)//
            {
                float temp = ip[y * in->stride] * centralPoint[y];
                sum += ((float)alpha* temp);
            }
            if (sum <= 0.0f) sum = 0.1f;
            else if (sum >= 255.0f) sum = 254.0f;
            
            *op = sum;
        }
}
//{ 0, 1, 2, 3, 4 }
//{ 0, 1, 2, 3, 4 }
//{ 0, 1, 2, 3, 4 }
//{ 0, 1, 2, 3, 4 }
//{ 0, 1, 2, 3, 4 }
//[5,0][4,0],3 0
//[5,1][4]
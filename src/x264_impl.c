#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <signal.h>

#include <x264.h>
#include "const.h"

FILE * fd;

int width, height;
x264_param_t param;
x264_picture_t pic;
x264_picture_t pic_out;
x264_t* encoder;
int i_frame;
int i_frame_size;
x264_nal_t *nal;
int i_nal;

void fill_x264_buf(uint32_t *pixels)
{
    // int luma_size = width * height;
    // int chroma_size = luma_size / 4;
    // TODO : better setup size

    // convert to YUV
    uint8_t *Y = pic.img.plane[0];
    uint8_t *U = pic.img.plane[1];
    uint8_t *V = pic.img.plane[2];

    #define R(i) ((pixels[i]&0x00FF0000)>>16)
    #define G(i) ((pixels[i]&0x0000FF00)>>8)
    #define B(i) ((pixels[i]&0x000000FF))

    for(int i = 0; i < SCREEN_WIDTH * SCREEN_HEIGHT; i++)
        Y[i] = ((66 * R(i) + 129 * G(i) + 25 * B(i) + 128) >> 8) + 16;
        
    for(int i = 0; i < SCREEN_WIDTH; i+=2)
        for(int j = 0; j < SCREEN_HEIGHT; j+=2){
            int w = SCREEN_WIDTH;
            // uint8_t r = (R(j*w+i) + R(j*w+i+1) + R((j+1)*w+i) + R((j+1)*w+i+1)) / 4;
            // uint8_t g = (G(j*w+i) + G(j*w+i+1) + G((j+1)*w+i) + G((j+1)*w+i+1)) / 4;
            // uint8_t b = (B(j*w+i) + B(j*w+i+1) + B((j+1)*w+i) + B((j+1)*w+i+1)) / 4;
            uint8_t r = R(j*w+i);
            uint8_t g = G(j*w+i);
            uint8_t b = B(j*w+i);
            
            U[i/2 + j*w/4] = ((-38*r -74*g +112*b + 128)>>8)+128;
            V[i/2 + j*w/4] = ((112*r-94*g-18*b+128)>>8)+128;
        }
        
    #undef R
    #undef G
    #undef B
}

void intHandler(int dummy) {
    close_x264();
    exit(0);
}

void init_x264(char * filename)
{
    fd = fopen(filename, "wb");
    width = SCREEN_WIDTH;
    height = SCREEN_HEIGHT;
    i_frame = 0;

    // param
    x264_param_default_preset(&param, "medium", NULL);
    param.i_threads = 1;
    param.i_width = width;
    param.i_height = height;
    param.i_fps_num = 25; //fps
    param.i_fps_den = 1;
    // Intra refres:
    param.i_keyint_max = 25; //fps
    param.b_intra_refresh = 1;
    //Rate control:
    param.rc.i_rc_method = X264_RC_CRF;
    param.rc.f_rf_constant = 25;
    param.rc.f_rf_constant_max = 35;
    //For streaming:
    param.b_repeat_headers = 1;
    param.b_annexb = 1;
    x264_param_apply_profile(&param, "high");

    x264_picture_alloc( &pic, param.i_csp, param.i_width, param.i_height );
    encoder = x264_encoder_open(&param);

    signal(SIGINT, intHandler);
}

int update_x264(uint32_t *pixels)
{
    fill_x264_buf(pixels);

    pic.i_pts = i_frame++;
    i_frame_size = x264_encoder_encode( encoder, &nal, &i_nal, &pic, &pic_out );
    if( i_frame_size < 0 )
        return 1;
    else if( i_frame_size )
    {
        if( !fwrite( nal->p_payload, i_frame_size, 1, fd ) )
            return 1;
    }
    return 0;
}

void close_x264()
{
    while( x264_encoder_delayed_frames( encoder ) )
    {
        i_frame_size = x264_encoder_encode( encoder, &nal, &i_nal, NULL, &pic_out );
        if( i_frame_size < 0 )
            break;
        else if( i_frame_size )
        {
            if( !fwrite( nal->p_payload, i_frame_size, 1, fd ) )
                break;
        }
    }

    x264_encoder_close( encoder );
    x264_picture_clean( &pic );
    fclose(fd);
}

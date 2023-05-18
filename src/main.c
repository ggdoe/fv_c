#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_events.h>

#define SCREEN_WIDTH 640
#define SCREEN_HEIGHT 480
// #define SCREEN_WIDTH 2560
// #define SCREEN_HEIGHT 1600

int main(int argc, char ** argv)
{   
    bool quit = false;
    SDL_Event event;
    SDL_Init(SDL_INIT_VIDEO);

    SDL_Window * window = SDL_CreateWindow("SDL2 Pixel Drawing",
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);
    SDL_Renderer * renderer = SDL_CreateRenderer(window, -1, 0);
    SDL_Texture * texture = SDL_CreateTexture(renderer,
            SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STATIC, SCREEN_WIDTH, SCREEN_HEIGHT);
    
    Uint32 *pixels = malloc(SCREEN_WIDTH * SCREEN_HEIGHT * sizeof(Uint32));

    ///

    bool leftMouseButtonDown = false;
    memset(pixels, 255, SCREEN_HEIGHT * SCREEN_WIDTH * sizeof(Uint32));
    while (!quit)
    {
        SDL_UpdateTexture(texture, NULL, pixels, SCREEN_WIDTH * sizeof(Uint32));

        SDL_RenderClear(renderer);
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);

        // SDL_Delay(20);
        while(SDL_PollEvent(&event)){
            switch (event.type)
            {
                case SDL_QUIT:
                    quit = true;
                    break;
                case SDL_MOUSEBUTTONUP:
                    if (event.button.button == SDL_BUTTON_LEFT)
                        leftMouseButtonDown = false;
                    break;
                case SDL_MOUSEBUTTONDOWN:
                    if (event.button.button == SDL_BUTTON_LEFT)
                        leftMouseButtonDown = true;
                case SDL_MOUSEMOTION:
                    if (leftMouseButtonDown)
                    {
                        int mouseX = event.motion.x;
                        int mouseY = event.motion.y;

                        int width, height;
                        SDL_GetWindowSize(window, &width, &height);
                        mouseX /= ((float)width  / SCREEN_WIDTH);
                        mouseY /= ((float)height / SCREEN_HEIGHT);

                        for(int j=-10; j < 10; j++)
                        for(int i=-10; i < 10; i++){
                            if(i*i + j*j < 50)
                                continue;
                            if(i*i + j*j > 100)
                                continue;
                            if((mouseY + 10) * SCREEN_WIDTH + mouseX + 10 > SCREEN_WIDTH * SCREEN_HEIGHT)
                                break;
                            if((mouseY - 10) * SCREEN_WIDTH + mouseX - 10 < 0)
                                break;
                            pixels[(mouseY+j) * SCREEN_WIDTH + mouseX + i] = rand();
                        }
                        // printf("%d %d\n", mouseX, mouseY);
                    }
                    break;
            }
        }
    }

    free(pixels);
    SDL_DestroyTexture(texture);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
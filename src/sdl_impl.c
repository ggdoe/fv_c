#include <SDL2/SDL.h>
// #include <SDL2/SDL_events.h>

#include "const.h"

SDL_Window * window;
SDL_Renderer * renderer;
SDL_Texture * texture;
SDL_Event event;
extern bool quit;

void init_sdl()
{
    SDL_Init(SDL_INIT_VIDEO);

    window = SDL_CreateWindow("fv_c",
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, 0);
    renderer = SDL_CreateRenderer(window, -1, 0);
    texture = SDL_CreateTexture(renderer,
            SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STATIC, SCREEN_WIDTH, SCREEN_HEIGHT);
}

void update_sdl(uint32_t *pixels)
{
    SDL_UpdateTexture(texture, NULL, pixels, SCREEN_WIDTH * sizeof(uint32_t));

    SDL_RenderClear(renderer);
    SDL_RenderCopy(renderer, texture, NULL, NULL);
    SDL_RenderPresent(renderer);

    while(SDL_PollEvent(&event)){
        switch (event.type)
        {
            case SDL_QUIT:
                quit = true;
                break;
        }
    }
    SDL_Delay(20);
}

void close_sdl()
{
    SDL_DestroyTexture(texture);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}


#include <SDL2/SDL.h>
// #include <SDL2/SDL_events.h>

#include "const.h"

SDL_Window * window;
SDL_Renderer * renderer;
SDL_Texture * texture;
SDL_Event event;
SDL_Rect render_rect;
extern bool quit;

inline static 
void update_render_rect()
{
    int w, h;
    SDL_GetWindowSize(window, &w, &h);
    int min = (w < h) ? w : h;
    render_rect.x = (w - min)/2;
    render_rect.y = (h - min)/2;
    render_rect.w = min;
    render_rect.h = min;
    
    SDL_RenderClear(renderer);
}

void init_sdl()
{
    SDL_Init(SDL_INIT_VIDEO);

    window = SDL_CreateWindow("fv_c",
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_RESIZABLE | SDL_WINDOW_SHOWN);
    renderer = SDL_CreateRenderer(window, -1, 0);
    texture = SDL_CreateTexture(renderer,
            SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STATIC, SCREEN_WIDTH, SCREEN_HEIGHT);
    update_render_rect();
}

void update_sdl(uint32_t *pixels)
{
    SDL_UpdateTexture(texture, NULL, pixels, SCREEN_WIDTH * sizeof(uint32_t));

    // SDL_RenderClear(renderer);
    
    SDL_RenderCopy(renderer, texture, NULL, &render_rect);
    SDL_RenderPresent(renderer);

    while(SDL_PollEvent(&event)){
        switch (event.type)
        {
            case SDL_QUIT:
                quit = true;
                break;
            case SDL_WINDOWEVENT:
                if(event.window.event == SDL_WINDOWEVENT_SIZE_CHANGED)
                    update_render_rect();
                break;
        }
    }
    // SDL_Delay(20);
}

void close_sdl()
{
    SDL_DestroyTexture(texture);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}


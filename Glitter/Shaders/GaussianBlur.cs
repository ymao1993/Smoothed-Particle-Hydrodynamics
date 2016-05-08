#version 450 core

//image size
layout(location = 0) uniform int width;
layout(location = 1) uniform int height;

//input&output image texture
layout(binding = 0, rgba32f) readonly uniform image2D inputImg;
layout(binding = 1, rgba32f) writeonly uniform image2D outputImg;

//local work group size
layout(local_size_x = 32,
       local_size_y = 32) in;

void main()
{
    //acquire ID
    vec2 coord = gl_GlobalInvocationID.xy;

    //only compute pixels falling inside the image
    if(coord.x < width && coord.y < height)
    {
      vec4 value = vec4(0,0,0,0); 
      value += imageLoad(inputImg, ivec2(coord.x - 1, coord.y - 1)) * 0.0625;
      value += imageLoad(inputImg, ivec2(coord.x    , coord.y - 1)) * 0.125;
      value += imageLoad(inputImg, ivec2(coord.x + 1, coord.y - 1)) * 0.0625;
      value += imageLoad(inputImg, ivec2(coord.x -1 , coord.y))     * 0.125;
      value += imageLoad(inputImg, ivec2(coord.x    , coord.y))     * 0.25;
      value += imageLoad(inputImg, ivec2(coord.x +1 , coord.y))     * 0.125;
      value += imageLoad(inputImg, ivec2(coord.x -1 , coord.y + 1)) * 0.0625;
      value += imageLoad(inputImg, ivec2(coord.x    , coord.y + 1)) * 0.125;
      value += imageLoad(inputImg, ivec2(coord.x +1 , coord.y + 1)) * 0.0625;
      imageStore(outputImg, ivec2(coord), value);
    }



}

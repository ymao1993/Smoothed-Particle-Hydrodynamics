#version 430 core

layout(location = 2) uniform mat4 w2v;
layout(location = 4) uniform sampler2D normalMap;

out vec4 color;

void main()
{  
   //compute point coord
   vec2 pointCoord = 2 * (gl_PointCoord - vec2(0.5,0.5));
   
   //compute normal (in camera space)
   float mag = dot(pointCoord.xy, pointCoord.xy);
   if(mag > 1) discard;
   
   //fetch normal
   vec3 normal = texelFetch(normalMap, ivec2(gl_FragCoord.xy), 0).xyz;
   color = vec4(normal,1);

}


#version 430 core

layout(location = 2) uniform mat4 w2v;

out vec4 color;

void main()
{  
   //compute point coord
   vec2 pointCoord = 2 * (gl_PointCoord - vec2(0.5,0.5)) * vec2(1,-1);
   
   //compute normal (in camera space)
   float mag = dot(pointCoord.xy, pointCoord.xy);
   if(mag > 1) discard;
   vec3 normal = vec3(pointCoord, sqrt(1-mag));
   
   //visualize normal
   color = vec4(normal,1);
}


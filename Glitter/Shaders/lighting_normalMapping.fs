#version 430 core

layout(location = 2) uniform mat4 w2v;

//lighting
layout(location = 4) uniform vec3 ambient;
layout(location = 5) uniform vec3 diffuse;
layout(location = 6) uniform vec3 specular;
layout(location = 7) uniform vec3 ldir;  //light direction in view space

uniform sampler2D normalMap;

in vec3 vpos; //view space position

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

   //phong lighting
   color = vec4(ambient,1);
   color += max(0,dot(normal, ldir)) * vec4(diffuse,1);
   //vec3 eye = -normalize(vpos);
   //color = pow(max(0,dot(eye, reflect(ldir, normal))),10) * vec4(specular,1);
   //color = vec4(specular,1);
   color.a = 1;

}


//#version 120
//
// simple.frag
//
//uniform mat4 mat;
uniform sampler2D texture0;
uniform vec3 lightdir;
varying vec3 normal;
varying vec4 color;
varying vec2 texcoord;

void main (void)
{
  gl_FragColor = max(texture2D(texture0, texcoord),color.xyzw);
  float a = 0.5; //coefficient of ambient light
  float d = a; // Reflection coefficient
  gl_FragColor *= d;

}


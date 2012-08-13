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
  //gl_FragColor = max(texture2D(texture0, texcoord),color.xyzw);
  gl_FragColor = texture2D(texture0, texcoord);
  //gl_FragColor = vec4(0.3, 0.8, 0.2, 1.0);
  //vec2 dfdx = abs(dFdx(texcoord));
  //vec2 dfdy = abs(dFdy(texcoord));
  //float d = max(max(dfdx.x,dfdx.y),max(dfdy.x,dfdy.y));
  //gl_FragColor *= d;
}


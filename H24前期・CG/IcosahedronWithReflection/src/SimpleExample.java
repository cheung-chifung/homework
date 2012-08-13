import javax.media.opengl.*;
import com.jogamp.opengl.util.PMVMatrix; 

public class SimpleExample extends SimpleExampleBase{
  Object3D obj;
  PMVMatrix mats;
  Shader shader;
  int uniformMat;
  int uniformLight;
  float t=0;

  public SimpleExample(){
    obj = new Icosahedron();
    mats = new PMVMatrix();
    shader = new Shader("resource/simple.vert", "resource/simple.frag");
  }

  public void init(GLAutoDrawable drawable){
    drawable.setGL(new DebugGL2(drawable.getGL().getGL2()));
    final GL2 gl = drawable.getGL().getGL2();
    //drawable.getGL().getGL2();
    gl.glViewport(0, 0, SCREENW, SCREENH);

    // Clear color buffer with black
    // gl.glClearColor(1.0f, 0.5f, 1.0f, 1.0f);
    gl.glClearColor(.0f, .0f, .0f, 1.0f);
    gl.glClearDepth(1.0f);
    gl.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT);
    gl.glEnable(GL2.GL_DEPTH_TEST);
    gl.glPixelStorei(GL.GL_UNPACK_ALIGNMENT,1);

    gl.glCreateShader(GL2GL3.GL_VERTEX_SHADER);
    shader.init(gl);
    int programName =shader.getID();
    gl.glBindAttribLocation(programName,Object3D.VERTEXPOSITION, "inposition");
    gl.glBindAttribLocation(programName,Object3D.VERTEXCOLOR, "incolor");
    gl.glBindAttribLocation(programName,Object3D.VERTEXNORMAL, "innormal");
    gl.glBindAttribLocation(programName,Object3D.VERTEXTEXCOORD0,"intexcoord0");
    shader.link(gl);
    uniformMat = gl.glGetUniformLocation(programName, "mat");
    uniformLight = gl.glGetUniformLocation(programName, "lightdir");
    gl.glUseProgram(programName);
    gl.glUniform3f(uniformLight, 0f, 10f, -10f);
    obj.init(gl, mats, programName);
    gl.glUseProgram(0);

  }

  public void display(GLAutoDrawable drawable){
    final GL2 gl = drawable.getGL().getGL2();
    gl.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT);
    mats.glMatrixMode(GL2.GL_MODELVIEW);
    mats.glLoadIdentity();
    mats.glTranslatef(0f,0f,-2.0f);
    //t = t+0.5f;
    t=40f;
    mats.glRotatef(t,0f,1f,0f);
    mats.glMatrixMode(GL2.GL_PROJECTION);
    mats.glLoadIdentity();
    mats.glFrustumf(-1f,1f,-1f,1f,1f,100f);
    mats.update();
    gl.glUseProgram(shader.getID());
    gl.glUniformMatrix4fv(uniformMat, 3, false, mats.glGetPMvMviMatrixf());

    obj.display(gl, mats);
    gl.glFlush();
    gl.glUseProgram(0);
  }

  public static void main(String[] args){
    SimpleExample t = new SimpleExample();
    t.start();
  }

}

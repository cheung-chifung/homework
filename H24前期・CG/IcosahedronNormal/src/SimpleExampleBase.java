import javax.swing.*;
import javax.swing.event.*;
import java.lang.reflect.InvocationTargetException;
import javax.media.opengl.*;
import javax.media.opengl.awt.GLCanvas;
import com.jogamp.opengl.util.*;
import java.awt.*;

public abstract class SimpleExampleBase implements GLEventListener{
  private JFrame frame;
  private GLCanvas canvas;

  protected final int SCREENW=640, SCREENH=640;

  public SimpleExampleBase(){
  }

  public void start(){
    try {
      java.awt.EventQueue.invokeAndWait(new Runnable(){
          public void run(){
            initGLUI();
          }
        }
      );
      java.awt.EventQueue.invokeAndWait(new Runnable(){
          public void run(){
            Animator animator = new Animator(canvas);
            animator.start();
          }
        }
      );
    } catch(InvocationTargetException e){
      System.out.println(e);
      e.printStackTrace();
    } catch(InterruptedException e){
      System.out.println(e);
    }
  }

  private void initGLUI(){
    GLProfile profile = GLProfile.get(GLProfile.GL2);
    GLCapabilities capabilities = new GLCapabilities(profile); 
    capabilities.setSampleBuffers(true); //This is effective if supported
    capabilities.setNumSamples(4);
    canvas = new GLCanvas(capabilities);
    canvas.addGLEventListener(this);
    frame = new JFrame("Icosahedron");
    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    Container container = frame.getContentPane();
    //    container.setLayout(new GroupLayout
    frame.add(canvas,BorderLayout.NORTH);
  
    canvas.setPreferredSize(new Dimension(SCREENW, SCREENH));
    frame.pack();
    frame.setResizable(false);
    frame.setVisible(true);
    canvas.setVisible(true);
  }

  public void reshape(GLAutoDrawable drawable,
                        int x, int y, 
                        int width, int height){
  }

  public void dispose(GLAutoDrawable drawable){
  }


}

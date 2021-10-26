import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib import pylab

def draw_line(x0, x1, y0, y1, angle=0.0, xrot=0.0, yrot=0.0,
              linewidth=1.0, color='k'):
    # Draw a line from (x0, y0) to (x1, y1)
    if angle == 0: plt.plot([x0, x1], [y0, y1],
                            linewidth=linewidth, color=color)
    else:
        sin_theta = np.sin(-angle * np.pi / 180)
        cos_theta = np.cos(-angle * np.pi / 180)
        x0_off = x0 - xrot
        x1_off = x1 - xrot
        y0_off = y0 - yrot
        y1_off = y1 - yrot
        x0_new = (x0_off * cos_theta - y0_off * sin_theta) + xrot
        x1_new = (x1_off * cos_theta - y1_off * sin_theta) + xrot
        y0_new = (x0_off * sin_theta + y0_off * cos_theta) + yrot
        y1_new = (x1_off * sin_theta + y1_off * cos_theta) + yrot
        plt.plot([x0_new, x1_new], [y0_new, y1_new],
                 linewidth=linewidth, color=color)

def draw_hriz_line(y, x0, x1, angle=0, xrot=0, yrot=0,
                   linewidth=1.0, color='k'):
    # Draw a horizontal line from (x0, y) to (x1, y)
    draw_line(x0, x1, y, y, angle, xrot, yrot, linewidth, color)

def draw_vert_line(x, y0, y1, angle=0, xrot=0, yrot=0,
                   linewidth=1.0, color='k'):
    # Draw a vertical line from (x, y0) to (x, y1)
    draw_line(x, x, y0, y1, angle, xrot, yrot, linewidth, color)

def draw_rectangle(x0, x1, y0, y1, angle=0, xrot=0, yrot=0,
                   linewidth=1.0, color='k'):
    # Draw a rectangle with opposite corners at (x0, y0) and (x1, y1)
    draw_hriz_line(y0, x0, x1, angle, xrot, yrot, linewidth, color)
    draw_hriz_line(y1, x0, x1, angle, xrot, yrot, linewidth, color)
    draw_vert_line(x0, y0, y1, angle, xrot, yrot, linewidth, color)
    draw_vert_line(x1, y0, y1, angle, xrot, yrot, linewidth, color)

def draw_circle(x, y, r, linewidth=1.0, color='k'):
    # Draw a circle with a center at (x, y) and a radius of r
    circ = pylab.Circle((x, y), radius=r, fill=False,
                        linewidth=linewidth, color=color)
    pylab.gca().add_patch(circ)

def draw_arc(x, y, r, t0, t1, linewidth=1.0, color='k'):
    # Draw an arc with a center at (x, y) and a radius of r and an angular span
    # in degrees from t0 to t1
    if t0 < 0 and t1 > 0:
        arc = patches.Arc((x, y), 2*r, 2*r, angle=0,
                          theta1=0, theta2=t1, fill=False,
                          linewidth=linewidth, color=color)
        pylab.gca().add_patch(arc)
        arc = patches.Arc((x, y), 2*r, 2*r, angle=0,
                          theta1=t0 + 360, theta2=360, fill=False,
                          linewidth=linewidth, color=color)
        pylab.gca().add_patch(arc)
    else:
        arc = patches.Arc((x, y), 2*r, 2*r, angle=0,
                          theta1=t0, theta2=t1, fill=False,
                          linewidth=linewidth, color=color)
        pylab.gca().add_patch(arc)

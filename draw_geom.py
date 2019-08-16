import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib import pylab

def draw_geometry():
    L_mod = 5.0

    # Iron
    draw_vert_line( 6.35, -7.5 ,   7.5)
    draw_vert_line( 6.65, -7.5 ,   7.5)
    draw_vert_line( 5.6 , -7.5 , -20.0)
    draw_vert_line( 5.6 ,  7.5 ,  20.0)
    draw_vert_line( 7.4 , -7.5 , -20.0)
    draw_vert_line( 7.4 ,  7.5 ,  20.0)
    draw_hriz_line(-7.5 ,  6.35,   5.6)
    draw_hriz_line(-7.5 ,  6.65,   7.4)
    draw_hriz_line( 7.5 ,  6.35,   5.6)
    draw_hriz_line( 7.5 ,  6.65,   7.4)

    # Moderator
    draw_hriz_line( 12.5, -20.0, 5.6)
    draw_hriz_line(-12.5, -20.0, 5.6)

    # Target chamber
    draw_rectangle(-20.0, 6.05 - L_mod, -3.0, 3.0)
    draw_rectangle(-20.0, 6.35 - L_mod, -3.3, 3.3)

    # Tally zone
    draw_rectangle(7.4, 7.9, -5.0, 5.0)

# Draw a line from (x0, y0) to (x1, y1)
def draw_line(x0, x1, y0, y1, angle=0, xrot=0, yrot=0):
    if angle == 0: plt.plot([x0, x1], [y0, y1])
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
        plt.plot([x0_new, x1_new], [y0_new, y1_new])

# Draw a horizontal line from (x0, y) to (x1, y)
def draw_hriz_line(y, x0, x1, angle=0, xrot=0, yrot=0):
    draw_line(x0, x1, y, y, angle, xrot, yrot)

# Draw a vertical line from (x, y0) to (x, y1)
def draw_vert_line(x, y0, y1, angle=0, xrot=0, yrot=0):
    draw_line(x, x, y0, y1, angle, xrot, yrot)

# Draw a rectangle with opposite corners at (x0, y0) and (x1, y1)
def draw_rectangle(x0, x1, y0, y1, angle=0, xrot=0, yrot=0):
    draw_hriz_line(y0, x0, x1, angle, xrot, yrot)
    draw_hriz_line(y1, x0, x1, angle, xrot, yrot)
    draw_vert_line(x0, y0, y1, angle, xrot, yrot)
    draw_vert_line(x1, y0, y1, angle, xrot, yrot)

# Draw a circle with a center at (x, y) and a radius of r
def draw_circle(x, y, r):
    circ = pylab.Circle((x, y), radius=r, fill=False)
    pylab.gca().add_patch(circ)

# Draw an arc with a center at (x, y) and a radius of r and an angular span
# in degrees from t0 to t1
def draw_arc(x, y, r, t0, t1):
    if t0 < 0 and t1 > 0:
        arc = patches.Arc((x, y), 2*r, 2*r, angle=0,
                          theta1=0, theta2=t1, fill=False)
        pylab.gca().add_patch(arc)
        arc = patches.Arc((x, y), 2*r, 2*r, angle=0,
                          theta1=t0 + 360, theta2=360, fill=False)
        pylab.gca().add_patch(arc)
    else:
        arc = patches.Arc((x, y), 2*r, 2*r, angle=0,
                          theta1=t0, theta2=t1, fill=False)
        pylab.gca().add_patch(arc)

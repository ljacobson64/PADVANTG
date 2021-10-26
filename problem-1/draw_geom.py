from draw_geom_funcs import *

def draw_geometry():
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
    draw_rectangle(-20.0, 1.05, -3.0, 3.0)
    draw_rectangle(-20.0, 1.35, -3.3, 3.3)

    # Tally zone
    draw_rectangle(7.4, 7.9, -5.0, 5.0)
